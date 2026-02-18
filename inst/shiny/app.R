# dosopt Shiny Application
# Full interactive interface for dose optimization design
# Launch via: dosopt::run_dosopt_app()

# Check required packages
required_pkgs <- c("shiny", "shinythemes", "DT", "shinyBS")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                       quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop("Please install required packages: ",
       paste(missing_pkgs, collapse = ", "),
       "\n  install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))",
       call. = FALSE)
}

library(shiny)
library(shinythemes)
library(DT)
library(shinyBS)
library(dosopt)

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- navbarPage(
  title     = "dosopt: Dose Optimization Design",
  theme     = shinythemes::shinytheme("flatly"),
  id        = "navbar",

  # ── Tab 1: General Utility Design ──────────────────────────────────────────
  tabPanel(
    "Utility Design",
    sidebarLayout(
      sidebarPanel(
        h4("Design Parameters"),
        sliderInput("pL",    "Efficacy rate (lower dose, p_L):", 0.05, 0.95, 0.30, 0.01),
        sliderInput("qL",    "Toxicity rate (lower dose, q_L):", 0.05, 0.95, 0.70, 0.01),
        numericInput("delta","Efficacy margin (δ):",  value = 0.10, min = 0.01, max = 0.40, step = 0.01),
        numericInput("d",    "Toxicity margin (d):",  value = 0.15, min = 0.01, max = 0.40, step = 0.01),
        sliderInput("phi",   "Efficacy-safety correlation (φ):", -0.5, 0.5, 0, 0.05),
        numericInput("alpha_L", "Target PCS (lower dose, α_L):", value = 0.80, min = 0.5, max = 0.99, step = 0.01),
        numericInput("alpha_H", "Target PCS (higher dose, α_H):", value = 0.80, min = 0.5, max = 0.99, step = 0.01),
        selectInput("method", "Calculation method:",
                    choices = c("Normal approximation" = "approx", "Exact multinomial" = "exact")),
        actionButton("calc_utility", "Calculate", class = "btn-primary btn-block"),
        hr(),
        helpText("Utility weights are derived from margins: u₂ = 1/(1+r), u₃ = r/(1+r), r = δ/d.")
      ),
      mainPanel(
        h4("Results"),
        verbatimTextOutput("utility_result"),
        hr(),
        h4("Utility Weights"),
        tableOutput("utility_weights"),
        hr(),
        h4("Feasible Correlation Range"),
        tableOutput("phi_range")
      )
    )
  ),

  # ── Tab 2: ROSE Design ─────────────────────────────────────────────────────
  tabPanel(
    "ROSE Design",
    sidebarLayout(
      sidebarPanel(
        h4("ROSE Parameters"),
        sliderInput("rose_pL",    "Efficacy rate (p_L):", 0.05, 0.95, 0.40, 0.01),
        numericInput("rose_delta","Efficacy margin (δ):", value = 0.15, min = 0.01, max = 0.40, step = 0.01),
        numericInput("rose_alpha_L", "α_L:", value = 0.80, min = 0.5, max = 0.99, step = 0.01),
        numericInput("rose_alpha_H", "α_H:", value = 0.80, min = 0.5, max = 0.99, step = 0.01),
        selectInput("rose_method", "Method:",
                    choices = c("Normal approximation" = "approx", "Exact binomial" = "exact")),
        actionButton("calc_rose", "Calculate", class = "btn-primary btn-block"),
        hr(),
        helpText("ROSE is the efficacy-only special case (u = (1,1,0,0)).")
      ),
      mainPanel(
        h4("Results"),
        verbatimTextOutput("rose_result")
      )
    )
  ),

  # ── Tab 3: Bias & Type I Error ─────────────────────────────────────────────
  tabPanel(
    "Bias & Type I Error",
    sidebarLayout(
      sidebarPanel(
        h4("Stage 1 Parameters"),
        sliderInput("b_pL",    "p_L:", 0.05, 0.95, 0.30, 0.01),
        sliderInput("b_qL",    "q_L:", 0.05, 0.95, 0.70, 0.01),
        numericInput("b_delta","δ:", value = 0.10, min = 0.01, max = 0.40, step = 0.01),
        numericInput("b_d",    "d:", value = 0.15, min = 0.01, max = 0.40, step = 0.01),
        numericInput("b_n1",   "n₁ (Stage 1 per arm):", value = 60, min = 10, max = 500),
        numericInput("b_n2",   "n₂ (Stage 2 total):", value = 120, min = 10, max = 1000),
        numericInput("b_alpha","Confirmatory α:", value = 0.05, min = 0.01, max = 0.20, step = 0.005),
        selectInput("b_test",  "Confirmatory test:",
                    choices = c("Z-test" = "z", "Exact binomial" = "binomial")),
        actionButton("calc_bias", "Calculate", class = "btn-primary btn-block")
      ),
      mainPanel(
        h4("Selection Bias"),
        verbatimTextOutput("bias_result"),
        hr(),
        h4("Type I Error Inflation"),
        verbatimTextOutput("t1e_result")
      )
    )
  ),

  # ── Tab 4: Formulas ────────────────────────────────────────────────────────
  tabPanel(
    "Formulas & Methods",
    fluidPage(
      h3("Key Formulas"),
      h4("Utility Score"),
      withMathJax(helpText("$$U = u_1 X(1-Y) + u_2 XY + u_3(1-X)(1-Y) + u_4(1-X)Y$$")),
      h4("Canonical Utility Weights"),
      withMathJax(helpText("$$u_2 = \\frac{1}{1+r},\\quad u_3 = \\frac{r}{1+r},\\quad r = \\frac{\\delta}{d}$$")),
      h4("Sample Size (Normal Approximation)"),
      withMathJax(helpText("$$n = \\left[\\frac{z_{\\alpha_L}\\sqrt{v_L} - z_{1-\\alpha_H}\\sqrt{v_H}}{\\Delta\\mu_H - \\Delta\\mu_L}\\right]^2$$")),
      h4("Selection Bias"),
      withMathJax(helpText("$$\\text{Bias}_1 = \\frac{\\text{Cov}(X,U)}{\\sigma_U\\sqrt{n_1}}\\cdot\\frac{1}{\\sqrt{\\pi}}\\exp\\!\\left(-\\frac{\\lambda_u^2 n_1}{4\\sigma_U^2}\\right)$$")),
      h4("Type I Error Inflation (Z-test)"),
      withMathJax(helpText("$$\\text{TIE} = 1 - \\Phi\\!\\left(z_{1-\\alpha} - \\frac{\\text{Bias}_{\\text{comb}}}{\\text{SE}_0}\\right)$$")),
      hr(),
      h4("References"),
      tags$ul(
        tags$li("Gu X, Xu C, Xu L, Yuan Y (2025). A Utility Score-Based Dose Optimization Design Under FDA's Project Optimus Initiative. Manuscript in preparation."),
        tags$li("Wang et al. (2025). ROSE: Randomized dose Optimization design with a Selection rule based on Efficacy."),
        tags$li("FDA Project Optimus: https://www.fda.gov/patients/drug-development-process/step-3-clinical-research#Optimus")
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Utility Design
  utility_res <- eventReactive(input$calc_utility, {
    req(input$pL, input$qL, input$delta, input$d)
    if (input$method == "approx") {
      calc_sample_size_utility_approx(
        pL = input$pL, qL = input$qL,
        delta = input$delta, d = input$d,
        phi = input$phi,
        alpha_L = input$alpha_L, alpha_H = input$alpha_H
      )
    } else {
      calc_sample_size_utility_exact(
        pL = input$pL, qL = input$qL,
        delta = input$delta, d = input$d,
        phi = input$phi,
        alpha_L = input$alpha_L, alpha_H = input$alpha_H,
        max_n = 500
      )
    }
  })

  output$utility_result <- renderPrint({
    req(utility_res())
    print(utility_res())
  })

  output$utility_weights <- renderTable({
    req(input$delta, input$d)
    u <- calc_utility(r = input$delta / input$d)
    data.frame(
      Outcome = c("Eff+, Tox-", "Eff+, Tox+", "Eff-, Tox-", "Eff-, Tox+"),
      Weight  = round(u, 4)
    )
  })

  output$phi_range <- renderTable({
    req(input$pL, input$qL)
    b <- phi_bounds(input$pL, input$qL)
    data.frame(
      Parameter    = c("φ_min", "φ_max"),
      Value        = round(as.numeric(b), 4)
    )
  })

  # ROSE Design
  rose_res <- eventReactive(input$calc_rose, {
    if (input$rose_method == "approx") {
      calc_sample_size_rose_approx(
        pL = input$rose_pL, delta = input$rose_delta,
        alpha_L = input$rose_alpha_L, alpha_H = input$rose_alpha_H
      )
    } else {
      calc_sample_size_rose_exact(
        pL = input$rose_pL, delta = input$rose_delta,
        alpha_L = input$rose_alpha_L, alpha_H = input$rose_alpha_H,
        max_n = 500
      )
    }
  })

  output$rose_result <- renderPrint({
    req(rose_res())
    print(rose_res())
  })

  # Bias & Type I Error
  bias_res <- eventReactive(input$calc_bias, {
    ss_tmp <- calc_sample_size_utility_approx(
      pL = input$b_pL, qL = input$b_qL,
      delta = input$b_delta, d = input$b_d,
      phi = 0, alpha_L = 0.80, alpha_H = 0.80
    )
    calc_bias(
      pL = input$b_pL, qL = input$b_qL,
      delta = input$b_delta, d = input$b_d,
      phi = 0,
      n1 = input$b_n1, n2 = input$b_n2,
      lambda_u = ss_tmp$lambda_u
    )
  })

  output$bias_result <- renderPrint({
    req(bias_res())
    print(bias_res())
  })

  output$t1e_result <- renderPrint({
    req(bias_res())
    t1e <- calc_type1_error(
      bias_stage1 = bias_res()$bias_stage1,
      n1 = input$b_n1, n2 = input$b_n2,
      p0 = input$b_pL, alpha = input$b_alpha,
      test = input$b_test
    )
    print(t1e)
  })
}

shinyApp(ui = ui, server = server)
