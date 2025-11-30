################################################################################
## SHINY APP — CLAIM SEVERITY DISTRIBUTION FITTING
## Fits: Gamma, Lognormal, Normal, Pareto (Type II / Lomax)
## Includes: AIC/BIC dashboard, QQ plots, PDF plots, styled parameter table
################################################################################

library(shiny)
library(fitdistrplus)
library(MASS)
library(plotly)
library(reactable)

################################################################################
## Pareto (Type II / Lomax) distribution definitions
################################################################################

dparetoL <- function(x, shape, scale) {
  ifelse(x >= 0,
         shape/scale * (1 + x/scale)^(-(shape + 1)),
         0)
}

pparetoL <- function(q, shape, scale) {
  ifelse(q >= 0,
         1 - (1 + q/scale)^(-shape),
         0)
}

qparetoL <- function(p, shape, scale) {
  scale * ((1 - p)^(-1/shape) - 1)
}

rparetoL <- function(n, shape, scale) {
  qparetoL(runif(n), shape, scale)
}

################################################################################
## UI
################################################################################

ui <- fluidPage(
  titlePanel("Insurance Severity Distribution Fitting"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Claim Cost Data (CSV, one column)", accept = ".csv"),
      numericInput("nsim", "Or simulate N lognormal claims:", 2000),
      actionButton("go", "Run Fitting"),
      hr(),
      helpText("Fits Gamma, Lognormal, Normal, and Pareto Type II distributions.")
    ),
    
    mainPanel(
      tabsetPanel(
        
        # Tab 1 — Data summary
        tabPanel("Data Summary",
                 verbatimTextOutput("summary")
        ),
        
        # Tab 2 — Fitted parameters (tidy table)
        tabPanel("Fitted Parameters",
                 h4("Fitted Distribution Parameters"),
                 reactableOutput("params_table")
        ),
        
        # Tab 3 — AIC/BIC dashboard
        tabPanel("Goodness of Fit (AIC/BIC)",
                 plotlyOutput("aic_table"),
                 plotlyOutput("aic_bar"),
                 plotlyOutput("bic_bar")
        ),
        
        # Tab 4 — Interactive QQ plots
        tabPanel("Interactive QQ Plots",
                 h4("Gamma"), plotlyOutput("qq_gamma"),
                 h4("Lognormal"), plotlyOutput("qq_lognorm"),
                 h4("Normal"), plotlyOutput("qq_norm"),
                 h4("Pareto (Type II)"), plotlyOutput("qq_pareto")
        ),
        
        # Tab 5 — PDF comparison
        tabPanel("PDF Plots",
                 h4("PDF Comparison"),
                 plotlyOutput("pdf_plot")
        )
      )
    )
  )
)

################################################################################
## SERVER
################################################################################

server <- function(input, output, session) {
  
  ##############################################################################
  ## CLAIM DATA INPUT
  ##############################################################################
  
  claim_data <- eventReactive(input$go, {
    if (!is.null(input$file)) {
      df <- read.csv(input$file$datapath)
      x <- df[[1]]
      x <- x[x > 0]
    } else {
      x <- abs(rlnorm(input$nsim, meanlog = 9, sdlog = 1))
    }
    x
  })
  
  ##############################################################################
  ## FIT DISTRIBUTIONS
  ##############################################################################
  
  fits <- eventReactive(input$go, {
    x <- claim_data()
    
    fit_gamma <- fitdist(x, "gamma", method = "mle")
    fit_lognormal <- fitdist(x, "lnorm", method = "mle")
    fit_normal <- fitdist(x, "norm", method = "mle")
    fit_pareto <- fitdist(
      x, "paretoL",
      start = list(shape = 1.5, scale = median(x) / 2),
      method = "mle"
    )
    
    list(
      gamma = fit_gamma,
      lognormal = fit_lognormal,
      normal = fit_normal,
      pareto = fit_pareto
    )
  })
  
  ##############################################################################
  ## OUTPUT: DATA SUMMARY
  ##############################################################################
  
  output$summary <- renderPrint({
    summary(claim_data())
  })
  
  ##############################################################################
  ## OUTPUT: FITTED PARAMETERS TABLE (reactable)
  ##############################################################################
  
  output$params_table <- renderReactable({
    
    f <- fits()
    
    params_df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto (Type II)"),
      Parameter_1 = c(
        paste0("shape = ", round(f$gamma$estimate["shape"], 4)),
        paste0("meanlog = ", round(f$lognormal$estimate["meanlog"], 4)),
        paste0("mean = ", round(f$normal$estimate["mean"], 4)),
        paste0("shape = ", round(f$pareto$estimate["shape"], 4))
      ),
      Parameter_2 = c(
        paste0("rate = ", round(f$gamma$estimate["rate"], 4)),
        paste0("sdlog = ", round(f$lognormal$estimate["sdlog"], 4)),
        paste0("sd = ", round(f$normal$estimate["sd"], 4)),
        paste0("scale = ", round(f$pareto$estimate["scale"], 4))
      ),
      stringsAsFactors = FALSE
    )
    
    reactable(
      params_df,
      striped = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      resizable = TRUE,
      defaultColDef = colDef(align = "center"),
      columns = list(
        Distribution = colDef(name = "Distribution", align = "left")
      )
    )
  })
  
  ##############################################################################
  ## AIC/BIC DASHBOARD
  ##############################################################################
  
  output$aic_table <- renderPlotly({
    f <- fits()
    aic_df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      AIC = c(f$gamma$aic, f$lognormal$aic, f$normal$aic, f$pareto$aic),
      BIC = c(f$gamma$bic, f$lognormal$bic, f$normal$bic, f$pareto$bic)
    )
    
    plot_ly(
      type = "table",
      header = list(values = names(aic_df),
                    fill = list(color = "#4a90e2"),
                    font = list(color = "white")),
      cells = list(values = t(aic_df))
    )
  })
  
  output$aic_bar <- renderPlotly({
    f <- fits()
    df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      AIC = c(f$gamma$aic, f$lognormal$aic, f$normal$aic, f$pareto$aic)
    )
    plot_ly(df, x = ~Distribution, y = ~AIC, type = "bar") %>%
      layout(title = "AIC Comparison")
  })
  
  output$bic_bar <- renderPlotly({
    f <- fits()
    df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      BIC = c(f$gamma$bic, f$lognormal$bic, f$normal$bic, f$pareto$bic)
    )
    plot_ly(df, x = ~Distribution, y = ~BIC, type = "bar") %>%
      layout(title = "BIC Comparison")
  })
  
  ##############################################################################
  ## INTERACTIVE QQ PLOTS
  ##############################################################################
  
  make_qq_plot <- function(x, fit, qdist, params) {
    probs <- ppoints(length(x))
    sorted_data <- sort(x)
    q_theo <- do.call(qdist, c(list(probs), params))
    
    plot_ly() %>%
      add_markers(x = q_theo, y = sorted_data) %>%
      add_lines(x = q_theo, y = q_theo, line = list(dash = "dash")) %>%
      layout(
        xaxis = list(title = "Theoretical Quantiles"),
        yaxis = list(title = "Empirical Quantiles")
      )
  }
  
  output$qq_gamma <- renderPlotly({
    x <- claim_data()
    fit <- fits()$gamma
    make_qq_plot(x, fit, qgamma,
                 list(shape = fit$estimate["shape"], rate = fit$estimate["rate"]))
  })
  
  output$qq_lognorm <- renderPlotly({
    x <- claim_data()
    fit <- fits()$lognormal
    make_qq_plot(x, fit, qlnorm,
                 list(meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"]))
  })
  
  output$qq_norm <- renderPlotly({
    x <- claim_data()
    fit <- fits()$normal
    make_qq_plot(x, fit, qnorm,
                 list(mean = fit$estimate["mean"], sd = fit$estimate["sd"]))
  })
  
  output$qq_pareto <- renderPlotly({
    x <- claim_data()
    fit <- fits()$pareto
    make_qq_plot(x, fit, qparetoL,
                 list(shape = fit$estimate["shape"], scale = fit$estimate["scale"]))
  })
  
  ##############################################################################
  ## PDF COMPARISON PLOT
  ##############################################################################
  
  output$pdf_plot <- renderPlotly({
    x <- claim_data()
    f <- fits()
    
    xgrid <- seq(min(x), max(x), length.out = 400)
    
    pdf_gamma <- dgamma(xgrid,
                        shape = f$gamma$estimate["shape"],
                        rate  = f$gamma$estimate["rate"])
    
    pdf_lognorm <- dlnorm(xgrid,
                          meanlog = f$lognormal$estimate["meanlog"],
                          sdlog   = f$lognormal$estimate["sdlog"])
    
    pdf_norm <- dnorm(xgrid,
                      mean = f$normal$estimate["mean"],
                      sd   = f$normal$estimate["sd"])
    
    pdf_pareto <- dparetoL(xgrid,
                           shape = f$pareto$estimate["shape"],
                           scale = f$pareto$estimate["scale"])
    
    plot_ly(x = xgrid, y = pdf_gamma, type = "scatter", mode = "lines", name = "Gamma") %>%
      add_lines(y = pdf_lognorm, name = "Lognormal") %>%
      add_lines(y = pdf_norm, name = "Normal") %>%
      add_lines(y = pdf_pareto, name = "Pareto (Type II)") %>%
      layout(
        title = "PDF Comparison",
        xaxis = list(title = "Claim Cost"),
        yaxis = list(title = "Density")
      )
  })
}

################################################################################
## RUN APP
################################################################################

shinyApp(ui, server)

