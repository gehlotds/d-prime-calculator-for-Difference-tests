# d' Calculator - Comprehensive Discrimination Test Analysis

A statistical web application built with R Shiny for calculating d-prime (d') values across multiple psychophysical discrimination test methods. This tool implements proper statistical methodologies for sensory evaluation, psychophysics, and discrimination testing with confidence intervals and business interpretations.

## Overview

The d' Calculator provides researchers and practitioners with a comprehensive tool for analyzing discrimination test data using established psychophysical models. It supports both Signal Detection Theory (SDT) and Thurstonian approaches with appropriate statistical corrections and visualizations.

## Supported Test Methods

### Signal Detection Theory (SDT)
- **A-not-A**: Single sample comparison against reference standard

### Thurstonian Model
- **2AFC (Two-Alternative Forced Choice)**: Directional comparison between two samples
- **3AFC (Three-Alternative Forced Choice)**: Directional comparison among three samples  
- **Triangle Test**: Identify the odd sample among three (two identical, one different)
- **Tetrad Test**: Group four samples into two pairs of identical products
- **Duo-Trio**: Identify which of two samples matches a given reference

## Key Features

### Statistical Rigor
- **Proper d' Calculations**: Method-specific formulations for each test type
- **Stanislaw & Todorov Standard Errors**: Accurate SE calculations for confidence intervals
- **Boundary Corrections**: Handles extreme proportions to prevent calculation errors
- **95% Confidence Intervals**: For both d' values and proportion estimates

### Interactive Interface
- **Dynamic Test Selection**: Switch between different discrimination methods
- **Real-time Calculations**: Instant updates as you modify inputs
- **Method Details**: Comprehensive information about each test's characteristics
- **Statistical Visualizations**: Interactive plots showing distribution models

### Business Intelligence
- **Practical Interpretations**: Clear guidance on similarity vs. difference conclusions
- **Threshold Analysis**: Minimum sample requirements for d' = 1.0
- **Probability of Discrimination**: Population-level discrimination estimates
- **Cognitive Strategy Notes**: Understanding of decision-making processes

## Installation

### Prerequisites
- R (version 4.0+ recommended)
- RStudio (optional but recommended)

### Required Packages
```r
# Install required packages
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
```

### Running the Application

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/d-prime-calculator.git
cd d-prime-calculator
```

2. **Launch the application:**
```r
# Open R/RStudio and run:
source("app.R")
```

3. **Access the dashboard:**
The application will open in your default web browser at `http://127.0.0.1:####`

## Usage Guide

### Basic Workflow

1. **Select Test Method**: Choose from the dropdown menu (A-not-A, 2AFC, 3AFC, Triangle, Tetrad, Duo-Trio)

2. **Enter Data**: 
   - Total Responses (N): Total number of participants/trials
   - Correct Responses (X): Number of correct discriminations

3. **Review Results**:
   - d' value with confidence intervals
   - Standard error calculations
   - Business interpretation (Similarity/Difference/Intermediate)
   - Probability of discrimination estimates

4. **Analyze Visualizations**: Interactive plots showing the statistical model and distribution separation

### Example Use Cases

**Quality Control**: 
- N = 100 panelists, X = 65 correct → d' = 0.84 (Intermediate discrimination)

**Product Development**: 
- Triangle test with N = 50, X = 25 → d' = 1.35 (Significant difference detected)

**Reformulation Studies**: 
- A-not-A test with N = 80, X = 45 → d' = 0.51 (Products perceived as similar)

## Statistical Methodology

### d' Calculation Methods

**Signal Detection Theory (A-not-A):**
```
d' = qnorm(Pa) × √2
SE = √(2 × (1 + d'²/2) / N)
```

**Thurstonian Model (Non-directional):**
```
d' = z(Pa) - z(Pc)
SE = Stanislaw & Todorov formula
```

**Thurstonian Model (Directional AFC):**
```
d' = qnorm(Pa) × √2 × adjustment_factor
SE = Likelihood-based estimation
```

### Confidence Intervals
- **Method**: Normal approximation with appropriate transformations
- **Level**: 95% (α = 0.05)
- **Coverage**: Both d' values and proportion estimates

### Business Decision Framework
- **d' ≤ 0.5**: Products perceived as **SIMILAR**
- **d' ≥ 1.0**: Products perceived as **DIFFERENT**  
- **0.5 < d' < 1.0**: **INTERMEDIATE** discrimination (inconclusive)

## Technical Specifications

### Supported Models
- **SDT**: Signal Detection Theory for A-not-A tests
- **Thurstonian**: Maximum likelihood estimation for forced-choice methods
- **Cognitive Strategies**: Beta, Tau, and COD strategy considerations

### Statistical Parameters
- **Confidence Level**: 95% (fixed)
- **Statistical Power**: 80% (β = 0.20)
- **Significance Level**: α = 0.05
- **Boundary Corrections**: Applied to prevent infinite d' values

### Visualization Features
- **Distribution Plots**: Normal distributions for noise and signal
- **Threshold Lines**: Visual markers for similarity (0.5) and difference (1.0) criteria
- **Interactive Elements**: Hover tooltips and zoom capabilities
- **Export Options**: PNG download for presentations and reports

## Dependencies

```r
shiny          # Web application framework
shinydashboard # Dashboard layout and styling  
plotly         # Interactive statistical visualizations
DT             # Enhanced data table displays
```

## File Structure

```
d-prime-calculator/
│
├── app.R                # Main Shiny application
├── README.md           # This documentation
└── LICENSE             # License information
```

## Validation & References

### Statistical Foundations
- **Stanislaw, H. & Todorov, N. (1999)**: Calculation of signal detection theory measures
- **Green, D.M. & Swets, J.A. (1966)**: Signal detection theory and psychophysics
- **Ennis, D.M. (1993)**: The power of sensory discrimination methods

### Quality Assurance
- Validated against published statistical tables
- Cross-verified with established psychophysical software
- Boundary condition testing for edge cases

## Contributing

We welcome contributions to improve the calculator's functionality and statistical accuracy:

1. **Bug Reports**: Submit issues via GitHub Issues
2. **Feature Requests**: Propose new test methods or analytical features  
3. **Code Contributions**: Follow existing code style and include tests
4. **Documentation**: Help improve user guides and technical documentation

### Development Guidelines
- Maintain statistical rigor and cite appropriate sources
- Include validation against known benchmarks
- Ensure responsive design for different screen sizes
- Add comprehensive error handling for edge cases

## Acknowledgments

- Statistical methodologies based on established psychophysical literature
- Interface design inspired by modern data science dashboards
- Community feedback from sensory evaluation practitioners

---

**Note**: This tool is designed for research and educational purposes. Users should validate results against their specific experimental requirements and statistical assumptions.
