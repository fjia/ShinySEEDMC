# install.packages(c("shiny", "lavaan", "dplyr", "purrr"))
library(shiny)
library(lavaan)
library(dplyr)
library(purrr)

ui <- fluidPage(
  titlePanel("Design-Based Efficiency Evaluation for Linear Growth Model"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Step 1: Define Design Pool"),
      selectInput("model_type", "Model Type", choices = c("linear"), selected = "linear"),
      sliderInput("time_points", "Number of Time Points", min = 3, max = 10, value = 5),
      radioButtons("sample_input", "Sample Size Input Method",
                   choices = c("Fixed Sample Size" = "fixed", "Budget and Cost" = "budget")),
      conditionalPanel(
        condition = "input.sample_input == 'fixed'",
        numericInput("sample_size", "Sample Size", value = 1666, min = 10)
      ),
      conditionalPanel(
        condition = "input.sample_input == 'budget'",
        numericInput("budget", "Total Budget", value = 3000),
        numericInput("unit_cost", "Cost per Measurement per Person", value = 10)
      ),
      textInput("attrition", "Attrition Rate for Each Time Point (comma-separated, e.g. 0,0.05,0.1...)", value = "0,0.05,0.1,0.1,0.15"),
      
      actionButton("generate_designs", "Generate Designs"),
      br(),
      tableOutput("design_table"),
      
      hr(),
      h4("Step 2: Define Model Parameters"),
      numericInput("int_mean", "Intercept Mean", value = 39.457),
      numericInput("slope_mean", "Slope Mean", value = 8.063),
      numericInput("int_var", "Intercept Variance", value = 28.776),
      numericInput("slope_var", "Slope Variance", value = 8.201),
      numericInput("int_slope_cov", "Intercept-Slope Covariance", value = 0),
      numericInput("resid_var", "Residual Variance", value = 30),
      numericInput("seed", "Simulation Seed",  min = 111, max = 1111111, value = 12345),
      
      actionButton("fit_full_model", "Run Model with Full Design"),
      tableOutput("model_results"),
      
      hr(),
      h4("Step 3: Evaluate Efficiency"),
      actionButton("evaluate_efficiency", "Run All Designs"),
      tableOutput("efficiency_table")
    ),
    
    mainPanel()
  )
)

server <- function(input, output, session) {
  
  # Step 1: Generate designs
  design_pool <- eventReactive(input$generate_designs, {
    tp <- input$time_points
    attrition <- as.numeric(unlist(strsplit(input$attrition, ",")))
    if(length(attrition) != tp) {
      showNotification("Attrition vector length must match number of time points", type = "error")
      return(NULL)
    }
    
    ####
    source("Shiny_2_design_pool.R")
    if (input$sample_input == "fixed") {
      designs_pool <- designPool(type ="com", traj = "linear", 
                                 time = tp, nfix = input$sample_size, 
                                 attrition = attrition)
    } else {
     designs_pool <- designPool(type ="com", traj = "linear", 
                               time = tp, budget = input$budget, unitcost = input$unit_cost, 
                               attrition = attrition)
   }
   

  designs_df <- tibble(
    Design = unlist(designs_pool$Complete$patterns_notation),
    SampleSize = unlist(designs_pool$Complete$N)
  )
    
    designs_df
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
      ", paste(paste0("y", 1:tp, " ~~ ", input$resid_var, "*y", 1:tp), collapse = "\n")
    )
  })

 
  full_model_results <- eventReactive(input$fit_full_model, {
    # data <- simulate_data(1666, input$time_points, input$int_mean, input$slope_mean,
    #                       input$int_var, input$slope_var, input$resid_var)
    
    ######
    set.seed(input$seed)
    if (input$sample_input == "fixed") {
      full_n <- input$sample_size
    } else {
      if (input$model_type == "linear") {
        leastT <- 3
      } else if (input$model_type == "quadratic") {
        leastT <- 4
      }
      
      Tc <- leastT:input$time_points
      CR <- 1 ### assume cost ratio
      C1 <- input$unit_cost*CR
      full_n <- floor(input$budget/(C1 + (Tc - 1)*input$unit_cost)) 
    }
    
    data <- simulateData(lavaan_model_string(), sample.nobs = full_n)
    ######
      
    fit <- lavaan::growth(lavaan_ana_model_string(), data = data) #, fixed.x = FALSE)
    if (!lavaan::inspect(fit, "converged")) {
      return(data.frame(Parameter = NA, Estimate = NA, SE = NA))
    }
    
    results <- parameterEstimates(fit) %>%
      filter(lhs %in% c("i", "s") & op == "~1") %>%
      select(Parameter = lhs, Estimate = est, SE = se)
    
    results
    #head(data)
  })
  
  output$model_results <- renderTable({
    full_model_results()
  })
  
  # Step 3: Run model for all designs
  efficiency_results <- eventReactive(input$evaluate_efficiency, {
    designs <- design_pool()
    if (is.null(designs)) return(NULL)

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

  
    run_one_design <- function(one_design) {
      times <- as.numeric(strsplit(gsub("[{}]", "", one_design$Design), ",")[[1]])
      data <- simulateData(lavaan_model_string_all_design(times), sample.nobs = one_design$SampleSize)
      fit <- tryCatch(growth(lavaan_ana_model_string_all_design(times), data = data), error = function(e) NULL)
      if (is.null(fit) || !lavaan::inspect(fit, "converged")) return(NA)
      
      pe <- parameterEstimates(fit)
      Intercept_SE <- pe$se[pe$lhs == "i" & pe$op == "~1"]
      Slope_SE <- pe$se[pe$lhs == "s" & pe$op == "~1"]
      Intercept_Efficiency <- full_se[1] / Intercept_SE
      Slope_Efficiency <- full_se[2] / Slope_SE
      final.results <- data.frame(Design = one_design$Design, 
                  SampleSize = one_design$SampleSize, 
                  Intercept_SE = Intercept_SE, 
                  Slope_SE = Slope_SE, 
                  Intercept_Efficiency = Intercept_Efficiency, 
                  Slope_Efficiency= Slope_Efficiency )
      
    final.results
    }
    
    all_models_results <- apply(designs, 1, FUN = run_one_design)
  
  })

  output$efficiency_table <- renderTable({
    efficiency_results()
  })

 }

shinyApp(ui, server)
