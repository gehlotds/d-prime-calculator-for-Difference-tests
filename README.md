d' Calculator: An R Shiny App for Sensory Discrimination Analysis

Introduction:
The d' Calculator is an interactive R Shiny application designed for sensory scientists and researchers. This tool goes beyond a simple calculation, providing a comprehensive platform to analyze data from common sensory discrimination tests. It automatically applies the correct statistical models, offers detailed results with confidence intervals, and visualizes the underlying theory to promote a deeper understanding of sensory sensitivity.

Features:
Wide Range of Tests: Supports a variety of sensory tests, including A-not-A, 2AFC, 3AFC, Triangle, Tetrad, and Duo-Trio.

Methodology & Models: Automatically selects the appropriate statistical method (Signal Detection Theory or Thurstonian model) for each test.

Detailed Statistical Output: Provides the d' value, its standard error, 95% confidence intervals, and the probability of discrimination (Pd).

Interactive Visualization: A dynamic plotly graph visualizes the signal and noise distributions, helping you understand the statistical basis of the d' calculation.

Business Interpretation: Offers a simple, clear interpretation of the results to guide decision-making.

Minimum for d' = 1: Calculates the minimum number of correct responses needed to achieve a d' of 1, a key threshold for sensory difference.

Getting Started
To use this application, you will need to have R and RStudio installed.

Install R Packages:
Ensure you have the required R packages by running the following commands in your R console:

install.packages(c("shiny", "shinydashboard", "plotly", "DT"))

Run the App:
Open the d' calculator with explaination.R file in RStudio and click the "Run App" button. Alternatively, you can run the app directly from the R console:

library(shiny)
source("d' calculator with explaination.R")
shinyApp(ui = ui, server = server)

Statistical Methods
The app is built on established sensory science principles:

Signal Detection Theory (SDT) is used for the A-not-A test.

The Thurstonian Model is applied to tests like 2AFC, 3AFC, and Triangle, which are based on comparing stimuli along a single dimension.

Stanislaw & Todorov (1999) formulas are implemented for standard error calculations, ensuring statistical rigor.

This tool aims to provide an accurate, transparent, and easy-to-use solution for analyzing sensory discrimination data.
