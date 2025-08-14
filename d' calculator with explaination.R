# Load packages
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))

# Load required libraries
library(shiny)
library(shinydashboard)
library(plotly)
library(DT)

# Updated test data with corrected sequences and precomputed chance levels
TEST_DATA <- list(
  "A-not-A" = list(
    test_name = "A-not-A",
    objective = "Is this product the same or different from the reference?",
    reference_product = "Yes",
    samples_per_trial = 1,
    possible_sequences = 2,
    sequence_series = "AA (both samples same), AB (one sample different)",
    variance_sequence = "Low",
    memory_effect = "Medium",
    cognitive_strategy = "In case subjects are sufficiently familiar with the reference product they will use Beta strategy.",
    example = "Panelist receives one sample and compares it with the reference. They decide if its same or different??",
    pc = 0.50,  # Precomputed: 1/2
    model = "SDT"
  ),
  "2AFC" = list(
    test_name = "2AFC", 
    objective = "Which of the products is more intense?",
    reference_product = "No",
    samples_per_trial = 2,
    possible_sequences = 2,
    sequence_series = "AB, BA",
    variance_sequence = "Low",
    memory_effect = "Low",
    cognitive_strategy = "Usually Beta but when subjects are not familiar with reference, subjects can use Tau-COD.",
    example = "Panelist compares two samples and determines which one has stronger sweetness or aroma or any other specified sensory attribute.",
    pc = 0.50,  # Precomputed: 1/2
    model = "Thurstonian"
  ),
  "3AFC" = list(
    test_name = "3AFC",
    objective = "Which of the three products is most intense?",
    reference_product = "No",
    samples_per_trial = 3,
    possible_sequences = 3,
    sequence_series = "AAB, ABA, BAA",
    variance_sequence = "Medium",
    memory_effect = "Low",
    cognitive_strategy = "Usually Beta strategy for directional comparisons.",
    example = "Panelist compares two samples and determines which one has stronger sweetness or aroma or any other specified sensory attribute.",
    pc = 0.333,  # Precomputed: 1/3
    model = "Thurstonian"
  ),
  "Triangle" = list(
    test_name = "Triangle",
    objective = "Which of the three products is different from the other two?",
    reference_product = "No",
    samples_per_trial = 3,
    possible_sequences = 6,
    sequence_series = "AAB, ABA, BAA, BBA, BAB, ABB",
    variance_sequence = "High",
    memory_effect = "High",
    cognitive_strategy = "Tau COD. Re-tasting can cause for some subjects strategy shifts to Beta. This makes it complicated to model with SDT and therefore results in inaccurate d' estimates.",
    example = "Panelist receives three samples (two identical, one different). They identify which sample is the odd one out.",
    pc = 0.333,  # Precomputed: 1/3
    model = "Thurstonian"
  ),
  "Tetrad" = list(
    test_name = "Tetrad",
    objective = "Group the four samples into two pairs of identical products.",
    reference_product = "No",
    samples_per_trial = 4,
    possible_sequences = 6,
    sequence_series = "AABB, ABAB, ABBA, BAAB, BABA, BBAA (Four samples presented: 2 samples of one type and 2 of another)",
    variance_sequence = "High",
    memory_effect = "High",
    cognitive_strategy = "Complex cognitive task involving multiple comparisons and grouping decisions.",
    example = "Four cups of Beer (two from Brand A, two from Brand B) — group the ones that taste the same.",
    pc = 0.167,  # Precomputed: 1/6
    model = "Thurstonian"
  ),
  "Duo-Trio" = list(
    test_name = "Duo-Trio",
    objective = "Which of the two products is the reference?",
    reference_product = "Yes",
    samples_per_trial = 3,
    possible_sequences = 6,
    sequence_series = "RAB, RBA, ARB, ABR, BRA, BAR",
    variance_sequence = "High",
    memory_effect = "Low",
    cognitive_strategy = "Re-tasting can cause strategy shift to more optimal strategies, difficult to model free retesting.",
    example = "Panelist tastes a reference (R), then two samples (A & B). They decide which of the two samples matches the reference they just tasted.",
    pc = 0.50,  # Precomputed: 1/2
    model = "Thurstonian"
  )
)

# Stanislaw & Todorov Standard Error Calculations
calc_se_stanislaw_todorov <- function(n, pa, pc, method) {
  if (method == "Triangle" || method == "Duo-Trio") {
    # For non-directional tests using z(Pa) - z(Pc) formula
    za <- qnorm(pa)
    zc <- qnorm(pc)
    
    # Stanislaw & Todorov formula for non-directional tests
    phi_za <- dnorm(za)
    phi_zc <- dnorm(zc)
    
    se <- sqrt((phi_za^2 * pa * (1 - pa)) / (n * (phi_za^2 + phi_zc^2)))
    
  } else if (method %in% c("2AFC", "3AFC", "Tetrad")) {
    # For directional AFC tests using likelihood-based estimation
    if (method == "2AFC") {
      # Standard error for 2AFC
      se <- sqrt((pa * (1 - pa)) / (n * (dnorm(qnorm(pa)))^2))
    } else if (method == "3AFC") {
      # Standard error for 3AFC
      phi_pa <- dnorm(qnorm(pa))
      se <- sqrt((pa * (1 - pa)) / (n * phi_pa^2))
    } else { # Tetrad
      # Standard error for Tetrad
      phi_pa <- dnorm(qnorm(pa))
      se <- sqrt((pa * (1 - pa)) / (n * phi_pa^2))
    }
  } else {
    # SDT for A-not-A
    phi_pa <- dnorm(qnorm(pa))
    se <- sqrt((pa * (1 - pa)) / (n * phi_pa^2))
  }
  
  return(se)
}

# Calculate d-prime using appropriate method
calc_d_prime <- function(n, x, test_data) {
  if (n <= 0 || x < 0 || x > n) return(list(d_prime = NA, se = NA, pa = NA))
  
  pa <- x / n
  pc <- test_data$pc
  model <- test_data$model
  test_name <- test_data$test_name
  
  # Boundary corrections
  if (pa <= pc) pa <- pc + 0.001
  if (pa >= 1) pa <- 0.999
  
  if (model == "SDT") {
    # Signal Detection Theory for A-not-A
    d_prime <- qnorm(pa) * sqrt(2)
    se <- sqrt(2 * (1 + d_prime^2/2) / n)
    
  } else if (model == "Thurstonian") {
    if (test_name %in% c("Triangle", "Duo-Trio")) {
      # Non-directional tests: d' = z(Pa) - z(Pc)
      d_prime <- qnorm(pa) - qnorm(pc)
      se <- calc_se_stanislaw_todorov(n, pa, pc, test_name)
      
    } else if (test_name %in% c("2AFC", "3AFC", "Tetrad")) {
      # Likelihood-based estimation for directional tests
      if (test_name == "2AFC") {
        d_prime <- qnorm(pa) * sqrt(2)
      } else if (test_name == "3AFC") {
        # For 3AFC, use the relationship for 3-alternative forced choice
        d_prime <- qnorm(pa) * sqrt(2) * 1.128  # Adjustment factor for 3AFC
      } else { # Tetrad
        # For Tetrad test
        d_prime <- qnorm(pa) * sqrt(2) * 1.414  # Adjustment factor for Tetrad
      }
      se <- calc_se_stanislaw_todorov(n, pa, pc, test_name)
    }
  }
  
  list(d_prime = d_prime, se = se, pa = pa)
}

# Calculate confidence intervals
calc_confidence_intervals <- function(d_prime, se, pa, n) {
  alpha <- 0.05
  z_alpha <- qnorm(1 - alpha/2)
  
  # CI for d'
  d_prime_ci_lower <- d_prime - z_alpha * se
  d_prime_ci_upper <- d_prime + z_alpha * se
  
  # CI for proportion correct (using normal approximation)
  se_pa <- sqrt(pa * (1 - pa) / n)
  pa_ci_lower <- max(0, pa - z_alpha * se_pa)
  pa_ci_upper <- min(1, pa + z_alpha * se_pa)
  
  # Probability of discrimination (Pd = Pa - Pc for most tests)
  pd <- pa - 0.5  # Simplified approach
  pd_ci_lower <- pa_ci_lower - 0.5
  pd_ci_upper <- pa_ci_upper - 0.5
  
  list(
    d_prime_ci = c(d_prime_ci_lower, d_prime_ci_upper),
    pa_ci = c(pa_ci_lower, pa_ci_upper),
    pd = pd,
    pd_ci = c(pd_ci_lower, pd_ci_upper)
  )
}

# Calculate minimum proportion for d' = 1 and corresponding number
calc_min_proportion_d1 <- function(test_data, n) {
  model <- test_data$model
  test_name <- test_data$test_name
  pc <- test_data$pc
  
  if (model == "SDT") {
    # For A-not-A: d' = qnorm(pa) * sqrt(2)
    # So pa = pnorm(1 / sqrt(2))
    min_pa <- pnorm(1 / sqrt(2))
    
  } else if (model == "Thurstonian") {
    if (test_name %in% c("Triangle", "Duo-Trio")) {
      # d' = qnorm(pa) - qnorm(pc) = 1
      # So qnorm(pa) = 1 + qnorm(pc)
      min_pa <- pnorm(1 + qnorm(pc))
      
    } else if (test_name == "2AFC") {
      # d' = qnorm(pa) * sqrt(2) = 1
      min_pa <- pnorm(1 / sqrt(2))
      
    } else if (test_name == "3AFC") {
      # d' = qnorm(pa) * sqrt(2) * 1.128 = 1
      min_pa <- pnorm(1 / (sqrt(2) * 1.128))
      
    } else { # Tetrad
      # d' = qnorm(pa) * sqrt(2) * 1.414 = 1
      min_pa <- pnorm(1 / (sqrt(2) * 1.414))
    }
  }
  
  min_number <- ceiling(min_pa * n)  # Minimum number of correct responses
  
  return(list(proportion = min_pa, number = min_number))
}

# Business interpretation
interpret_business <- function(d_prime) {
  if (is.na(d_prime)) return("Invalid calculation")
  
  if (d_prime <= 0.5) {
    return("SIMILARITY: Products perceived as similar (d' ≤ 0.5)")
  } else if (d_prime >= 1.0) {
    return("DIFFERENCE: Products perceived as different (d' ≥ 1.0)")
  } else {
    return("INTERMEDIATE: Can not conclude Products perceived as different or similar")
  }
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "d' Calculator with Proper Statistical Methods"),
  
  dashboardSidebar(
    width = 330,
    h4("Choose Your Test Method", style = "color: white; padding: 10px;"),
    
    selectInput("test_selection", 
                label = NULL,
                choices = list(
                  "A-not-A (SDT)" = "A-not-A",
                  "2AFC (Thurstonian)" = "2AFC", 
                  "3AFC (Thurstonian)" = "3AFC",
                  "Triangle (Thurstonian)" = "Triangle",
                  "Tetrad (Thurstonian)" = "Tetrad",
                  "Duo-Trio (Thurstonian)" = "Duo-Trio"
                ),
                selected = "A-not-A"),
    
    br(),
    h5("Test Objective:", style = "color: white;"),
    div(style = "background: rgba(255,255,255,0.1); padding: 10px; border-radius: 5px; color: white;",
        textOutput("test_objective")
    ),
    
    br(),
    h5("Statistical Method:", style = "color: white;"),
    div(style = "background: rgba(255,255,255,0.1); padding: 10px; border-radius: 5px; color: white;",
        textOutput("statistical_method")
    ),
    
    br(),
    h5("Presentation Sequence Series:", style = "color: white;"),
    div(style = "background: rgba(255,255,255,0.1); padding: 10px; border-radius: 5px; color: white; font-size: 12px;",
        textOutput("sequence_series")
    ),
    
    br(),
    h5("Example:", style = "color: white;"),
    div(style = "background: rgba(255,255,255,0.1); padding: 10px; border-radius: 5px; color: white;",
        textOutput("test_example")
    )
  ),
  
  dashboardBody(
    withMathJax(),
    fluidRow(
      # Test Details and Results
      column(width = 8,
             box(title = "Test Methodology Details", status = "primary", solidHeader = TRUE, width = NULL,
                 tableOutput("methodology_table")
             ),
             
             # Results with Confidence Intervals
             box(title = "Statistical Results", status = "success", solidHeader = TRUE, width = NULL,
                 verbatimTextOutput("detailed_results")
             ),
             
             # SDT Visualization
             box(title = "Statistical Model Visualization", status = "info", solidHeader = TRUE, width = NULL,
                 plotlyOutput("sdt_plot", height = "350px")
             )
      ),
      
      # Calculator
      column(width = 4,
             box(title = "Input & Quick Results", status = "warning", solidHeader = TRUE, width = NULL,
                 numericInput("total_n", "Total Responses (N):", value = 100, min = 1),
                 numericInput("correct_x", "Correct Responses (X):", value = 65, min = 0),
                 
                 hr(),
                 h4("d' Value:"),
                 verbatimTextOutput("d_prime_result"),
                 
                 # Added d' definition
                 div(
                   style = "background: #f9f9f9; padding: 10px; border-radius: 5px; 
           margin: 10px 0; font-size: 14px;",
                   HTML("<b>d' Definition:</b><br>
       d' (sensitivity) is a measure to describe how well someone can tell two things apart—like spotting a difference of taste or texture or smell between two samples.<br>
       $$d' = \\frac{\\text{mean(signal)} - \\text{mean(noise)}}{SE}$$")
                 ),
                 
                 
                 hr(),
                 h4("Business Interpretation:"),
                 div(style = "background: #e8f4fd; padding: 10px; border-radius: 5px; font-weight: bold;",
                     textOutput("business_impact")
                 ),
                 
                 hr(),
                 h4("Min. for d' = 1:"),
                 div(style = "background: #f0f8ff; padding: 10px; border-radius: 5px;",
                     verbatimTextOutput("min_proportion")
                 )
             ),
             
             # Probability of Discrimination Info
             box(title = "Probability of Discrimination", status = "info", solidHeader = TRUE, width = NULL,
                 div(style = "font-size: 13px; line-height: 1.4;",
                     HTML("<b>Definition:</b><br>Probability of Discrimination (Pd) tells you what proportion of panelists can truly sense a difference between two samples, not just guess. It helps Sensory Scientists helps researchers understand how big the sensory difference truly is within the population tested.<br><br>
                     <b>Statistical Parameters (Fixed):</b><br>
                     • Confidence Level: 95%<br>
                     • Statistical Power: 80%<br>
                     • α = 0.05, β = 0.20"))
             )
      )
    ),
    
    fluidRow(
      column(width = 12,
             box(title = "Cognitive Strategy & Method Notes", status = "info", solidHeader = TRUE, width = NULL,
                 textOutput("cognitive_info")
             )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Get selected test data
  current_test <- reactive({
    TEST_DATA[[input$test_selection]]
  })
  
  # Update correct_x max when total_n changes
  observeEvent(input$total_n, {
    updateNumericInput(session, "correct_x", max = input$total_n)
  })
  
  # Test objective and example
  output$test_objective <- renderText({
    current_test()$objective
  })
  
  output$test_example <- renderText({
    current_test()$example
  })
  
  output$statistical_method <- renderText({
    test <- current_test()
    paste("Method:", test$model, 
          if (test$test_name %in% c("Triangle", "Duo-Trio")) "- Non-directional (z(Pa) - z(Pc))" 
          else if (test$test_name %in% c("2AFC", "3AFC", "Tetrad")) "- Likelihood-based estimation"
          else "- Signal Detection Theory")
  })
  
  output$sequence_series <- renderText({
    current_test()$sequence_series
  })
  
  # Methodology table
  output$methodology_table <- renderTable({
    test <- current_test()
    data.frame(
      Characteristic = c("Reference Product", "Samples per Trial", 
                         "Possible Sequences", "Variance", "Memory Effect", 
                         "Chance Probability", "Model"),
      Value = c(test$reference_product, test$samples_per_trial,
                test$possible_sequences, test$variance_sequence, test$memory_effect,
                paste0(test$pc), test$model)
    )
  }, striped = TRUE, hover = TRUE)
  
  # Detailed calculations with confidence intervals
  output$detailed_results <- renderText({
    req(input$total_n, input$correct_x)
    
    if (input$correct_x > input$total_n) {
      return("ERROR: Correct responses cannot exceed total responses")
    }
    
    test <- current_test()
    calc <- calc_d_prime(input$total_n, input$correct_x, test)
    
    if (is.na(calc$d_prime)) {
      return("Unable to calculate d'")
    }
    
    ci <- calc_confidence_intervals(calc$d_prime, calc$se, calc$pa, input$total_n)
    
    paste0("═══ STATISTICAL RESULTS ═══\n",
           "d' = ", round(calc$d_prime, 3), "\n",
           "95% CI for d': [", round(ci$d_prime_ci[1], 3), ", ", round(ci$d_prime_ci[2], 3), "]\n",
           "Standard Error (SE) = ", round(calc$se, 3), "\n\n",
           "Proportion Correct (Pa) = ", round(calc$pa, 3), "\n",
           "95% CI for Pa: [", round(ci$pa_ci[1], 3), ", ", round(ci$pa_ci[2], 3), "]\n",
           "Chance Level (Pc) = ", test$pc, "\n\n",
           "Probability of Discrimination (Pd) = ", round(ci$pd, 3), "\n",
           "95% CI for Pd: [", round(max(0, ci$pd_ci[1]), 3), ", ", round(min(1, ci$pd_ci[2]), 3), "]")
  })
  
  # Quick d' result
  output$d_prime_result <- renderText({
    req(input$total_n, input$correct_x)
    
    if (input$correct_x > input$total_n) {
      return("Check inputs")
    }
    
    test <- current_test()
    calc <- calc_d_prime(input$total_n, input$correct_x, test)
    
    if (is.na(calc$d_prime)) {
      return("Cannot calculate")
    }
    
    round(calc$d_prime, 3)
  })
  
  # Business impact
  output$business_impact <- renderText({
    req(input$total_n, input$correct_x)
    
    if (input$correct_x > input$total_n) return("Check inputs")
    
    test <- current_test()
    calc <- calc_d_prime(input$total_n, input$correct_x, test)
    interpret_business(calc$d_prime)
  })
  
  # Minimum proportion and number for d' = 1
  output$min_proportion <- renderText({
    req(input$total_n)
    
    test <- current_test()
    min_result <- calc_min_proportion_d1(test, input$total_n)
    paste0("Proportion: ", round(min_result$proportion, 3), " (", round(min_result$proportion * 100, 1), "%)\n",
           "Number: ", min_result$number, " out of ", input$total_n, " responses")
  })
  
  # Cognitive strategy
  output$cognitive_info <- renderText({
    test <- current_test()
    paste0("Cognitive Strategy: ", test$cognitive_strategy, "\n\n",
           "Statistical Notes: Using ", test$model, " model. ",
           if (test$test_name %in% c("Triangle", "Duo-Trio")) {
             "Non-directional test using d' = z(Pa) - z(Pc) formula with Stanislaw & Todorov SE."
           } else if (test$test_name %in% c("2AFC", "3AFC", "Tetrad")) {
             "Likelihood-based estimation with Stanislaw & Todorov SE formula."
           } else {
             "Signal Detection Theory with classical d' calculation."
           })
  })
  
  # Statistical Model Visualization
  output$sdt_plot <- renderPlotly({
    req(input$total_n, input$correct_x)
    
    if (input$correct_x > input$total_n) {
      return(plot_ly() %>% add_text(x = 0.5, y = 0.5, text = "Invalid inputs"))
    }
    
    test <- current_test()
    calc <- calc_d_prime(input$total_n, input$correct_x, test)
    
    if (is.na(calc$d_prime)) {
      return(plot_ly() %>% add_text(x = 0.5, y = 0.5, text = "Cannot calculate d'"))
    }
    
    d <- calc$d_prime
    se <- calc$se
    
    # Create distributions
    x <- seq(-3, 3 + d, 0.1)
    noise <- dnorm(x, 0, 1)
    signal <- dnorm(x, d, 1)
    
    # Calculate y_max for positioning annotations
    y_max <- max(c(noise, signal))
    
    plot_ly() %>%
      add_lines(x = x, y = noise, name = "Noise Distribution", 
                line = list(color = "blue", width = 2)) %>%
      add_lines(x = x, y = signal, name = "Signal Distribution", 
                line = list(color = "red", width = 2)) %>%
      add_lines(x = c(0.5, 0.5), y = c(0, y_max*0.3), name = "Similarity Threshold", 
                line = list(color = "orange", dash = "dot", width = 1)) %>%
      add_lines(x = c(1.0, 1.0), y = c(0, y_max*0.3), name = "Difference Threshold", 
                line = list(color = "purple", dash = "dot", width = 1)) %>%
      layout(
        title = paste0(test$model, " Model for ", test$test_name, " (d' = ", round(d, 2), ")"),
        xaxis = list(title = "Intensity"),
        yaxis = list(title = "Probability Density"),
        showlegend = TRUE,
        annotations = list(
          list(x = d + 0.2, y = y_max * 1.05, 
               text = paste("d' =", round(d, 2)), 
               showarrow = FALSE, font = list(color = "black", size = 14)),
          list(x = d + 0.2, y = y_max * 0.95, 
               text = paste("SE =", round(se, 3)), 
               showarrow = FALSE, font = list(color = "darkgreen", size = 12)),
          list(x = 0.5, y = y_max * 0.2, 
               text = "Similarity", 
               showarrow = FALSE, font = list(color = "orange", size = 10)),
          list(x = 1.0, y = y_max * 0.2, 
               text = "Difference", 
               showarrow = FALSE, font = list(color = "purple", size = 10))
        )
      )
  })
}


shinyApp(ui = ui, server = server)
