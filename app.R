################################################################################
## SHINY APP — CLAIM SEVERITY DISTRIBUTION FITTING (ROBUST VERSION)
## Fixes Pareto MLE failure in R "Great Square Root"
## Stable, gradient-free Pareto-II (Lomax) log-likelihood
################################################################################

library(shiny)
library(fitdistrplus)
library(MASS)
library(plotly)
library(reactable)

################################################################################
## Pareto Type II (Lomax): dpqr
################################################################################

dparetoL <- function(x, shape, scale) {
  ifelse(x >= 0, shape/scale * (1 + x/scale)^(-(shape + 1)), 0)
}

pparetoL <- function(q, shape, scale) {
  ifelse(q >= 0, 1 - (1 + q/scale)^(-shape), 0)
}

qparetoL <- function(p, shape, scale) {
  scale * ((1 - p)^(-1/shape) - 1)
}

rparetoL <- function(n, shape, scale) {
  qparetoL(runif(n), shape, scale)
}

################################################################################
## Robust Pareto-II Log-Likelihood (Fix)
################################################################################

lomax_loglik <- function(param, x) {
  shape <- param[1]
  scale <- param[2]
  
  if (shape <= 0 || scale <= 0) return(-Inf)
  
  n <- length(x)
  
  ll <- n * (log(shape) - log(scale)) -
    (shape + 1) * sum(log1p(x / scale))
  
  return(ll)
}

robust_pareto_fit <- function(x) {
  
  start <- c(shape = 1.5, scale = median(x) / 3)
  
  opt <- optim(
    par = start,
    fn = function(p) -lomax_loglik(p, x),
    method = "L-BFGS-B",
    lower = c(0.05, 0.0001),
    upper = c(50, max(x)),
    control = list(maxit = 500)
  )
  
  if (opt$convergence != 0) return(NULL)
  
  est <- opt$par
  ll  <- -opt$value
  
  fd <- list()
  class(fd) <- "fitdist"
  
  fd$estimate <- c(shape = est[1], scale = est[2])
  fd$loglik   <- ll
  fd$aic      <- -2 * ll + 2 * 2
  fd$bic      <- -2 * ll + log(length(x)) * 2
  fd$distname <- "paretoL"
  fd$data     <- x
  fd$method   <- "mme"
  
  fd
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
      helpText("Fits Gamma, Lognormal, Normal, and Pareto Type II (Lomax).")
    ),
    
    mainPanel(
      tabsetPanel(
        
        tabPanel("Data Summary",
                 verbatimTextOutput("summary")
        ),
        
        tabPanel("Fitted Parameters",
                 h4("Fitted Distribution Parameters"),
                 reactableOutput("params_table")
        ),
        
        tabPanel("Goodness of Fit (AIC/BIC)",
                 plotlyOutput("aic_table"),
                 plotlyOutput("aic_bar"),
                 plotlyOutput("bic_bar")
        ),
        
        tabPanel("Interactive QQ Plots",
                 h4("Gamma"),         plotlyOutput("qq_gamma"),
                 h4("Lognormal"),     plotlyOutput("qq_lognorm"),
                 h4("Normal"),        plotlyOutput("qq_norm"),
                 h4("Pareto (Type II)"), plotlyOutput("qq_pareto")
        ),
        
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
  ## Load or simulate claim data
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
  ## Distribution Fits
  ##############################################################################
  
  fits <- eventReactive(input$go, {
    
    x <- claim_data()
    
    fit_gamma    <- fitdist(x, "gamma", method = "mme")
    fit_lognorm  <- fitdist(x, "lnorm", method = "mle")
    fit_normal   <- fitdist(x, "norm", method = "mle")
    
    pareto_fit <- tryCatch(
      robust_pareto_fit(x),
      error = function(e) { message("Pareto fit failed: ", e$message); NULL }
    )
    
    list(
      gamma = fit_gamma,
      lognormal = fit_lognorm,
      normal = fit_normal,
      pareto = pareto_fit
    )
  })
  
  ##############################################################################
  ## Summary
  ##############################################################################
  
  output$summary <- renderPrint({
    summary(claim_data())
  })
  
  ##############################################################################
  ## Parameter Table
  ##############################################################################
  
  output$params_table <- renderReactable({
    
    f <- fits()
    
    params_df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto (Type II)"),
      Parameter_1 = c(
        paste0("shape = ", round(f$gamma$estimate["shape"], 4)),
        paste0("meanlog = ", round(f$lognormal$estimate["meanlog"], 4)),
        paste0("mean = ", round(f$normal$estimate["mean"], 4)),
        if (!is.null(f$pareto)) paste0("shape = ", round(f$pareto$estimate["shape"], 4)) else "Fit failed"
      ),
      Parameter_2 = c(
        paste0("rate = ", round(f$gamma$estimate["rate"], 4)),
        paste0("sdlog = ", round(f$lognormal$estimate["sdlog"], 4)),
        paste0("sd = ", round(f$normal$estimate["sd"], 4)),
        if (!is.null(f$pareto)) paste0("scale = ", round(f$pareto$estimate["scale"], 4)) else "Fit failed"
      ),
      stringsAsFactors = FALSE
    )
    
    reactable(params_df, striped = TRUE, bordered = TRUE,
              highlight = TRUE, defaultColDef = colDef(align = "center"))
  })
  
  ##############################################################################
  ## AIC/BIC Tables
  ##############################################################################
  
  output$aic_table <- renderPlotly({
    
    f <- fits()
    
    aic_df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      AIC = c(f$gamma$aic, f$lognormal$aic, f$normal$aic,
              ifelse(is.null(f$pareto), NA, f$pareto$aic)),
      BIC = c(f$gamma$bic, f$lognormal$bic, f$normal$bic,
              ifelse(is.null(f$pareto), NA, f$pareto$bic))
    )
    
    plot_ly(type = "table",
            header = list(values = names(aic_df)),
            cells = list(values = t(aic_df)))
  })
  
  output$aic_bar <- renderPlotly({
    
    f <- fits()
    df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      AIC = c(f$gamma$aic, f$lognormal$aic, f$normal$aic,
              ifelse(is.null(f$pareto), NA, f$pareto$aic))
    )
    
    plot_ly(df, x = ~Distribution, y = ~AIC, type = "bar") %>%
      layout(title = "AIC Comparison")
  })
  
  output$bic_bar <- renderPlotly({
    
    f <- fits()
    df <- data.frame(
      Distribution = c("Gamma", "Lognormal", "Normal", "Pareto"),
      BIC = c(f$gamma$bic, f$lognormal$bic, f$normal$bic,
              ifelse(is.null(f$pareto), NA, f$pareto$bic))
    )
    
    plot_ly(df, x = ~Distribution, y = ~BIC, type = "bar") %>%
      layout(title = "BIC Comparison")
  })
  
  ##############################################################################
  ## QQ Plots
  ##############################################################################
  
  make_qq_plot <- function(x, qdist, params) {
    probs <- ppoints(length(x))
    sorted <- sort(x)
    theo   <- do.call(qdist, c(list(probs), params))
    
    plot_ly() %>%
      add_markers(x = theo, y = sorted) %>%
      add_lines(x = theo, y = theo, line = list(dash = "dash"))
  }
  
  output$qq_gamma <- renderPlotly({
    f <- fits()
    make_qq_plot(claim_data(), qgamma,
                 list(shape = f$gamma$estimate["shape"],
                      rate  = f$gamma$estimate["rate"]))
  })
  
  output$qq_lognorm <- renderPlotly({
    f <- fits()
    make_qq_plot(claim_data(), qlnorm,
                 list(meanlog = f$lognormal$estimate["meanlog"],
                      sdlog    = f$lognormal$estimate["sdlog"]))
  })
  
  output$qq_norm <- renderPlotly({
    f <- fits()
    make_qq_plot(claim_data(), qnorm,
                 list(mean = f$normal$estimate["mean"],
                      sd   = f$normal$estimate["sd"]))
  })
  
  output$qq_pareto <- renderPlotly({
    f <- fits()
    if (is.null(f$pareto)) return(NULL)
    
    make_qq_plot(claim_data(), qparetoL,
                 list(shape = f$pareto$estimate["shape"],
                      scale = f$pareto$estimate["scale"]))
  })
  
  ##############################################################################
  ## PDF Comparison Plot
  ##############################################################################
  
  output$pdf_plot <- renderPlotly({
    
    x <- claim_data()
    f <- fits()
    
    xgrid <- seq(min(x), max(x), length.out = 400)
    
    pdf_gamma <- dgamma(xgrid, f$gamma$estimate["shape"], f$gamma$estimate["rate"])
    pdf_logn  <- dlnorm(xgrid, f$lognormal$estimate["meanlog"], f$lognormal$estimate["sdlog"])
    pdf_norm  <- dnorm(xgrid, f$normal$estimate["mean"], f$normal$estimate["sd"])
    
    pdf_par <- if (!is.null(f$pareto))
      dparetoL(xgrid, f$pareto$estimate["shape"], f$pareto$estimate["scale"])
    else
      rep(NA, length(xgrid))
    
    plot_ly(x = xgrid, y = pdf_gamma, type = "scatter", mode = "lines", name = "Gamma") %>%
      add_lines(y = pdf_logn, name = "Lognormal") %>%
      add_lines(y = pdf_norm, name = "Normal") %>%
      add_lines(y = pdf_par, name = "Pareto (Type II)") %>%
      layout(title = "PDF Comparison")
  })
}

################################################################################
## Run App
################################################################################

shinyApp(ui, server)
