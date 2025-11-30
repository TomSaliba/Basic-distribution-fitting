###############################################################################
## SHINY APP — CLAIM SEVERITY DISTRIBUTION FITTING
## Distributions: Gamma, Lognormal, Normal, Pareto Type II (Lomax)
## Includes: AIC/BIC dashboard, QQ plots, PDF plots, parameter table
## FIXED FOR DEPLOYMENT — Includes correct custom distribution registration
###############################################################################

library(shiny)
library(fitdistrplus)
library(MASS)
library(plotly)
library(reactable)

###############################################################################
## CUSTOM PARETO-II (LOMAX) DISTRIBUTION + REGISTRATION FOR FITDIST
###############################################################################

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

## Register with fitdistrplus — **critical fix**
paretoL_dist <- list(
  densfun = dparetoL,
  pfun    = pparetoL,
  qfun    = qparetoL,
  rfun    = rparetoL,
  start   = function(x) {
    list(
      shape = 1.5,
      scale = median(x) / 2
    )
  }
)

###############################################################################
## UI
###############################################################################

ui <- fluidPage(
  titlePanel("Insurance Severity Distribution Fitting"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Claim Cost Data (CSV, one column)", accept = ".csv"),
      numericInput("nsim", "Or simulate N lognormal claims:", 2000),
      actionButton("go", "Run Fitting"),
      hr(),
      helpText("Fits Gamma, Lognormal, Normal, and Pareto Type II (Lomax) distributions.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Summary", verbatimTextOutput("summary")),
        tabPanel("Fitted Parameters", h4("Fitted Distribution Parameters"), reactableOutput("params_table")),
        tabPanel("Goodness of Fit (AIC/BIC)", 
                 plotlyOutput("aic_table"),
                 plotlyOutput("aic_bar"),
                 plotlyOutput("bic_bar")),
        tabPanel("Interactive QQ Plots",
                 h4("Gamma"), plotlyOutput("qq_gamma"),
                 h4("Lognormal"), plotlyOutput("qq_lognorm"),
                 h4("Normal"), plotlyOutput("qq_norm"),
                 h4("Pareto (Type II)"), plotlyOutput("qq_pareto")),
        tabPanel("PDF Plots",
                 h4("PDF Comparison"),
                 plotlyOutput("pdf_plot"))
      )
    )
  )
)

###############################################################################
## SERVER
###############################################################################

server <- function(input, output, session) {
  
  ###########################################################################
  ## CLAIM DATA INPUT
  ###########################################################################
  
  claim_data <- eventReactive(input$go, {
    if (!is.null(input$file)) {
      df <- tryCatch(read.csv(input$file$datapath), error = function(e) NULL)
      validate(need(!is.null(df), "Error reading CSV file."))
      x <- df[[1]]
    } else {
      x <- abs(rlnorm(input$nsim, meanlog = 9, sdlog = 1))
    }
    
    x <- x[x > 0]   # safety filter
    validate(need(length(x) > 10, "Not enough positive values to fit distributions."))
    x
  })
  
  ###########################################################################
  ## FIT DISTRIBUTIONS
  ###########################################################################
  
  fits <- eventReactive(input$go, {
    x <- claim_data()
    
    fit_gamma     <- fitdist(x, "gamma")
    fit_lognorm   <- fitdist(x, "lnorm")
    fit_normal    <- fitdist(x, "norm")
    fit_pareto    <- fitdist(x, "paretoL", custom.dist = paretoL_dist)
    
    list(
      gamma     = fit_gamma,
      lognormal = fit_lognorm,
      normal    = fit_normal,
      pareto    = fit_pareto
    )
  })
  
  ###########################################################################
  ## OUTPUT: DATA SUMMARY
  ###########################################################################
  
  output$summary <- renderPrint({
    summary(claim_data())
  })
  
  ###########################################################################
  ## OUTPUT: PARAMETER TABLE
  ###########################################################################
  
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
      columns = list(Distribution = colDef(align = "left"))
    )
  })
  
  ###########################################################################
  ## AIC/BIC TABLE + BAR CHARTS
  ###########################################################################
  
  output$aic_table <- renderPlotly({
    f <- fits()
    
    df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      AIC = c(f$gamma$aic, f$lognormal$aic, f$normal$aic, f$pareto$aic),
      BIC = c(f$gamma$bic, f$lognormal$bic, f$normal$bic, f$pareto$bic)
    )
    
    plot_ly(type = "table",
            header = list(values = names(df),
                          fill = list(color = "#4a90e2"),
                          font = list(color = "white")),
            cells = list(values = t(df)))
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
  
  ###########################################################################
  ## GENERIC QQ-PLOT BUILDER
  ###########################################################################
  
  make_qq_plot <- function(x, qfun, params) {
    probs <- ppoints(length(x))
    sorted_data <- sort(x)
    q_theo <- do.call(qfun, c(list(probs), params))
    
    plot_ly() %>%
      add_markers(x = q_theo, y = sorted_data) %>%
      add_lines(x = q_theo, y = q_theo, line = list(dash = "dash")) %>%
      layout(
        xaxis = list(title = "Theoretical Quantiles"),
        yaxis = list(title = "Empirical Quantiles")
      )
  }
  
  output$qq_gamma <- renderPlotly({
    f <- fits()$gamma
    make_qq_plot(claim_data(), qgamma,
                 list(shape = f$estimate["shape"], rate = f$estimate["rate"]))
  })
  
  output$qq_lognorm <- renderPlotly({
    f <- fits()$lognormal
    make_qq_plot(claim_data(), qlnorm,
                 list(meanlog = f$estimate["meanlog"], sdlog = f$estimate["sdlog"]))
  })
  
  output$qq_norm <- renderPlotly({
    f <- fits()$normal
    make_qq_plot(claim_data(), qnorm,
                 list(mean = f$estimate["mean"], sd = f$estimate["sd"]))
  })
  
  output$qq_pareto <- renderPlotly({
    f <- fits()$pareto
    make_qq_plot(claim_data(), qparetoL,
                 list(shape = f$estimate["shape"], scale = f$estimate["scale"]))
  })
  
  ###########################################################################
  ## PDF COMPARISON PLOT
  ###########################################################################
  
  output$pdf_plot <- renderPlotly({
    x <- claim_data()
    f <- fits()
    
    xgrid <- seq(min(x), max(x), length.out = 400)
    
    pdf_gamma <- dgamma(xgrid, shape = f$gamma$estimate["shape"], rate = f$gamma$estimate["rate"])
    pdf_lognorm <- dlnorm(xgrid, meanlog = f$lognormal$estimate["meanlog"], sdlog = f$lognormal$estimate["sdlog"])
    pdf_norm <- dnorm(xgrid, mean = f$normal$estimate["mean"], sd = f$normal$estimate["sd"])
    pdf_pareto <- dparetoL(xgrid, shape = f$pareto$estimate["shape"], scale = f$pareto$estimate["scale"])
    
    plot_ly(x = xgrid, y = pdf_gamma, type = "scatter", mode = "lines", name = "Gamma") %>%
      add_lines(y = pdf_lognorm, name = "Lognormal") %>%
      add_lines(y = pdf_norm, name = "Normal") %>%
      add_lines(y = pdf_pareto, name = "Pareto (Type II)") %>%
      layout(title = "PDF Comparison",
             xaxis = list(title = "Claim Cost"),
             yaxis = list(title = "Density"))
  })
}

###############################################################################
## RUN APP
###############################################################################

shinyApp(ui, server)
