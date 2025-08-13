# install.packages(c("shiny", "lavaan", "dplyr", "purrr"))
library(shiny)
library(lavaan)
library(dplyr)
library(purrr)

ui <- fluidPage(
  titlePanel("ShinySEEDMC: Searching for Optimal Design(s) using Monte Carlo Relative Efficiency"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Step 1: Define Design Pool"),
      selectInput("model_type", "Model Type", choices = c("linear"), selected = "linear"),
      sliderInput("time_points", "Number of Time Points", min = 3, max = 10, value = 5),
      radioButtons("sample_input", "Sample Size Input Method",
                   choices = c("Budget and Cost" = "budget", "Fixed Sample Size" = "fixed")),
      conditionalPanel(
        condition = "input.sample_input == 'budget'",
        numericInput("budget", "Total Budget", value = 100000),
        numericInput("unit_cost", "Cost per Measurement per Person", value = 20)
      ),
      conditionalPanel(
        condition = "input.sample_input == 'fixed'",
        numericInput("sample_size", "Sample Size", value = 1000, min = 10)
      ),
      
      # textInput("attrition", "Attrition Rate for Each Time Point (comma-separated, e.g. 0,0.05,0.1...)", 
      #           value = "0, 0.075, 0.15, 0.225, 0.3"),
      
      actionButton("generate_designs", "Generate Designs"),
      br(),
      tableOutput("design_table"),
      
      hr(),
      h4("Step 2: Define Model Parameters"),
      numericInput("int_mean", "Intercept Mean", value = 39.457),
      numericInput("slope_mean", "Slope Mean", value = 8.063),
      numericInput("int_var", "Intercept Variance", value = 28.776),
      numericInput("slope_var", "Slope Variance", value = 8.201),
      numericInput("int_slope_cov", "Intercept-Slope Covariance", value = 1.56),
      numericInput("resid_var", "Residual Variance", value = 30),
      
      # textInput("attrition", "Attrition Rate for Each Time Point (comma-separated, e.g. 0,0.05,0.1...)", 
      #           value = "0, 0.075, 0.15, 0.225, 0.3"),
      
      actionButton("fit_full_model", "Model Check (Full Design and N = 1000)"),
      tableOutput("model_results"),
      
      hr(),
      h4("Step 3: Evaluate Relative Efficiency"),
      numericInput("num_reps", "Number of Replications (large number, e.g., 5000, recommended)", value = 50),
      numericInput("seed", "Simulation Seed",  min = 100, max = 1000000, value = 12345),
      actionButton("evaluate_efficiency", "Run All Designs --> "),
      # tableOutput("efficiency_table")
    ),
    
    mainPanel(
        tableOutput("efficiency_table")
      )
    )
)

server <- function(input, output, session) {
  
  # Step 1: Generate designs
  design_pool <- eventReactive(input$generate_designs, {
    tp <- input$time_points
    # attrition <- as.numeric(unlist(strsplit(input$attrition, ",")))
    # if(length(attrition) != tp) {
    #   showNotification("Attrition vector length must match number of time points", type = "error")
    #   return(NULL)
    # }
    
    ####
    source("source_ShinySEEDMC_design_pool.R")
    if (input$sample_input == "fixed") {
      designs_pool <- designPool(type ="com", traj = "linear", 
                                 time = tp, nfix = input$sample_size)
                                 # attrition = attrition)
    } else {
     designs_pool <- designPool(type ="com", traj = "linear", 
                               time = tp, budget = input$budget, unitcost = input$unit_cost) 
                               # attrition = attrition)
   }
   

  designs_tb <- tibble(
    Design = unlist(designs_pool$Complete$patterns_notation),
    SampleSize = unlist(designs_pool$Complete$N)
  )
    
    designs_tb
  })
  
  output$design_table <- renderTable({
    design_pool()
  })
  
  # Step 2: Lavaan model and estimation with full design
  
  lavaan_model_string <- reactive({
    tp <- input$time_points
    time_vals <- 0:(tp - 1)
    paste0("
      i =~ ", paste(rep("1*", tp), paste0("y", 1:tp), collapse = " + "), "
      s =~ ", paste(time_vals, "*", paste0("y", 1:tp), collapse = " + "), "
      i ~ ", input$int_mean, "*1
      s ~ ", input$slope_mean, "*1
      i ~~ ", input$int_var, "*i
      s ~~ ", input$slope_var, "*s
      i ~~ ", input$int_slope_cov, "*s
      ", paste(paste0("y", 1:tp, " ~~ ", input$resid_var, "*y", 1:tp), collapse = "\n")
    )
  })
  
  lavaan_ana_model_string <- reactive({
    tp <- input$time_points
    time_vals <- 0:(tp - 1)
    paste0("
      i =~ ", paste(rep("1*", tp), paste0("y", 1:tp), collapse = " + "), "
      s =~ ", paste(time_vals, "*", paste0("y", 1:tp), collapse = " + "), "
      i ~ 1
      s ~ 1
      i ~~ i
      s ~~ s
      i ~~ s
      ", paste(paste0("y", 1:tp, " ~~ y", 1:tp), collapse = "\n")
    )
  })

 
  full_model_results <- eventReactive(input$fit_full_model, {
    
    ######
    set.seed(input$seed)
    if (input$sample_input == "fixed") {
      full_n <- input$sample_size
    } else {
      Tc <- input$time_points
      CR <- 1 ### assume cost ratio
      C1 <- input$unit_cost*CR
      full_n <- input$budget/(C1 + (Tc - 1)*input$unit_cost)
    }
    
    data <- simulateData(lavaan_model_string(), sample.nobs = 1000) #### full_n
    ######
      
    fit <- lavaan::growth(lavaan_ana_model_string(), data = data) #, fixed.x = FALSE)
    if (!lavaan::inspect(fit, "converged")) {
      return(data.frame(Parameter = NA, Estimate = NA, SE = NA))
    }
    
    results <- parameterEstimates(fit) %>%
      filter(lhs %in% c("i", "s") & op == "~1") %>%
      select(Parameter = lhs, Estimate = est, SE = se)
    
    results
  })
  
  output$model_results <- renderTable({
    full_model_results()
  })
  
  # Step 3: Run model for all designs
  efficiency_results <- eventReactive(input$evaluate_efficiency, {
    designs <- design_pool()
    if (is.null(designs)) return(NULL)
    
    designs_df <- data.frame(designs)

    full_se <- full_model_results()$SE
    tp <- input$time_points
    
    lavaan_model_string_all_design <- function(times){
      # times <- as.numeric(strsplit(gsub("[{}]", "", designnotation), ",")[[1]])
      time_scores <- times - min(times)  # Align to start at 0 for slope
      paste0("
      i =~ ", paste0("1*y", times, collapse = " + "), "
      s =~ ", paste0(time_scores, "*y", times, collapse = " + "), "
      i ~ ", input$int_mean, "*1
      s ~ ", input$slope_mean, "*1
      i ~~ ", input$int_var, "*i
      s ~~ ", input$slope_var, "*s
      i ~~ ", input$int_slope_cov, "*s
      ", paste(paste0("y", times, " ~~ ", input$resid_var, "*y", times), collapse = "\n"))
      }
    
    lavaan_ana_model_string_all_design <- function(times){
      # times <- as.numeric(strsplit(gsub("[{}]", "", designnotation), ",")[[1]])
      time_scores <- times - min(times)  # Align to start at 0 for slope
      paste0("
      i =~ ", paste0("1*y", times, collapse = " + "), "
      s =~ ", paste0(time_scores, "*y", times, collapse = " + "), "
      i ~ 1
      s ~ 1
      i ~~ i
      s ~~ s
      i ~~ s
      ", paste(paste0("y", times, " ~~ y", times), collapse = "\n")
      )}

    set.seed(12345)
    seed.list <- sample(1:10000000, input$num_reps)
    
    run_one_rep <- function(rep_seed){
      run_one_design <- function(one_design) {
        times <- as.numeric(strsplit(gsub("[{}]", "", unlist(one_design["Design"])), ",")[[1]])
        set.seed(rep_seed)
        data <- simulateData(lavaan_model_string_all_design(times), sample.nobs = as.numeric(one_design["SampleSize"]))
        fit <- tryCatch(growth(lavaan_ana_model_string_all_design(times), data = data), error = function(e) NULL)
        if (is.null(fit) || !lavaan::inspect(fit, "converged"))  Slope_EST <- NA #return(NA)
        
        pe <- parameterEstimates(fit)
        # Intercept_SE <- pe$se[pe$lhs == "i" & pe$op == "~1"]
        Slope_EST <- pe$est[pe$lhs == "s" & pe$op == "~1"]
        # Intercept_RE <- (full_se[1])^2 / (Intercept_SE)^2
        # Slope_RE <- (full_se[2])^2 / (Slope_SE)^2
        # final.results <- data.frame(Design = as.character(unlist(one_design["Design"])), 
        #             SampleSize = as.numeric(one_design["SampleSize"]), 
        #             Int_SE = Intercept_SE, 
        #             Slope_SE = Slope_SE, 
        #             Int_RE = Intercept_RE, 
        #             Slope_RE = Slope_RE)
        final.results <- Slope_EST
        
      final.results
      }
      
      all_models_results <- apply(designs_df, 1, FUN = run_one_design)
      ### all_models_results_tb <- tibble(do.call("rbind", all_models_results))
    }
    
    all_reps_results <- lapply(seed.list, FUN = run_one_rep)
    all_reps_results_v <- do.call("rbind", all_reps_results)
    all_reps_variances <- apply(all_reps_results_v, 2, FUN = var, na.rm = TRUE)
    all_reps_re <- round(all_reps_variances[length(all_reps_variances)]/all_reps_variances, 4)
    out <- data.frame(Design = as.character(unlist(designs_df[,"Design"])), 
                      SampleSize = as.character(designs_df[,"SampleSize"]),  
                      NumReps = as.character(colSums(!is.na(all_reps_results_v))),
                      RE_MeanSlope = all_reps_re)
  
    out
  })

   output$efficiency_table <- renderTable({
    efficiency_results()
  })


 }

shinyApp(ui, server)
