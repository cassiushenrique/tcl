# =============================================================================
# UFLA / ICTIN / GCT018 - DISTRIBUICOES AMOSTRAIS
# Media (TCL), proporcao (aproximacao normal) e distribuicao de S^2.
# Para populacao normal, a referencia de S^2 e uma qui-quadrado REESCALADA.
# Autor: material didatico para uso em Estatistica - versao revisada, 2026.
# Dependencias: shiny, ggplot2. Execute: shiny::runApp("app_TCL_unificado.R")
# =============================================================================

library(shiny)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Utilitarios numericos, reprodutibilidade e apresentacao
# -----------------------------------------------------------------------------
fmt <- function(x, digits = 4) {
  if (!is.finite(x)) return(as.character(x))
  formatC(x, format = "f", digits = digits, decimal.mark = ",", big.mark = ".")
}

fmt_pct <- function(x) paste0(fmt(100 * x, 2), "%")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# Nao usa HTML() nas expressoes: tags$span escapa '<' e '>' com seguranca.
math <- function(tex, display = FALSE) {
  tags$span(if (display) paste0("$$", tex, "$$")
            else paste0("\\(", tex, "\\)"))
}

# Quarto momento central obtido a partir dos momentos brutos.
central4 <- function(m1, m2, m3, m4) {
  m4 - 4 * m1 * m3 + 6 * m1^2 * m2 - 3 * m1^4
}

# Familia populacional com parametros, momentos EXATOS e densidade analitica.
make_model <- function(dist, input) {
  if (dist == "normal") {
    mu <- input$mu
    sigma <- input$sigma
    req(mu, sigma)
    variance <- sigma^2
    return(list(
      nome = "Normal", notacao = "X\\sim N(\\mu,\\sigma^2)",
      parametros = "\\mu\\in\\mathbb{R},\\;\\sigma>0;\\quad X\\in\\mathbb{R}",
      f = "f(x)=\\frac{1}{\\sigma\\sqrt{2\\pi}}\\exp\\!\\left[-\\frac{(x-\\mu)^2}{2\\sigma^2}\\right]",
      aplicacao = "Erros de medição e caracteristicas de qualidade aproximadamente simétricas. Se a população e normal, a média amostral e exatamente normal e a variância padronizada e exatamente qui-quadrado.",
      mu = mu, var = variance, mu4 = 3 * variance^2,
      r = function(m) rnorm(m, mean = mu, sd = sigma),
      d = function(x) dnorm(x, mean = mu, sd = sigma),
      q = function(p) qnorm(p, mean = mu, sd = sigma)
    ))
  }
  if (dist == "uniforme") {
    ab <- input$ab
    req(ab)
    validate(need(length(ab) == 2 && ab[2] > ab[1],
                  "Escolha limites com a < b para a distribuição uniforme."))
    a <- ab[1]; b <- ab[2]
    variance <- (b - a)^2 / 12
    return(list(
      nome = "Uniforme contínua", notacao = "X\\sim U(a,b)",
      parametros = "a,b\\in\\mathbb{R},\\;a<b;\\quad X\\in[a,b]",
      f = "f(x)=\\begin{cases}(b-a)^{-1},&a\\le x\\le b,\\\\0,&\\text{caso contrário},\\end{cases}",
      aplicacao = "Modelo de referencia quando valores num intervalo sao equiprovaveis por unidade de comprimento. A distribuição das médias tende a forma normal ao crescer n.",
      mu = (a + b) / 2, var = variance, mu4 = (b - a)^4 / 80,
      r = function(m) runif(m, min = a, max = b),
      d = function(x) dunif(x, min = a, max = b),
      q = function(p) qunif(p, min = a, max = b)
    ))
  }
  if (dist == "lognormal") {
    mul <- input$mul; sigmal <- input$sigmal
    req(mul, sigmal)
    raw <- function(k) exp(k * mul + k^2 * sigmal^2 / 2)
    mu <- raw(1); variance <- raw(2) - mu^2
    return(list(
      nome = "Lognormal (assimetria a direita)",
      notacao = "\\log X\\sim N(\\mu_{\\log},\\sigma_{\\log}^2)",
      parametros = "\\mu_{\\log}\\in\\mathbb{R},\\;\\sigma_{\\log}>0;\\quad X>0",
      f = "f(x)=\\frac{1}{x\\sigma_{\\log}\\sqrt{2\\pi}}\\exp\\!\\left[-\\frac{(\\log x-\\mu_{\\log})^2}{2\\sigma_{\\log}^2}\\right],\\quad x>0",
      aplicacao = "Tempos ou quantidades positivas com assimetria a direita. Com forte assimetria, a convergência visual das médias pode ser lenta mesmo para n relativamente grande.",
      mu = mu, var = variance,
      mu4 = central4(raw(1), raw(2), raw(3), raw(4)),
      r = function(m) rlnorm(m, meanlog = mul, sdlog = sigmal),
      d = function(x) dlnorm(x, meanlog = mul, sdlog = sigmal),
      q = function(p) qlnorm(p, meanlog = mul, sdlog = sigmal)
    ))
  }
  if (dist == "beta") {
    alpha <- input$alpha; beta <- input$beta
    req(alpha, beta)
    raw <- function(k) prod((alpha + 0:(k - 1)) /
                             (alpha + beta + 0:(k - 1)))
    mu <- alpha / (alpha + beta)
    variance <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
    return(list(
      nome = "Beta", notacao = "X\\sim\\operatorname{Beta}(\\alpha,\\beta)",
      parametros = "\\alpha>0,\\;\\beta>0;\\quad 0<X<1",
      f = "f(x)=\\frac{x^{\\alpha-1}(1-x)^{\\beta-1}}{B(\\alpha,\\beta)},\\quad 0<x<1",
      aplicacao = "Modelagem de proporções e indices continuamente mensurados no intervalo (0,1). Diferentes formas da distribuição permitem estudar assimetria e comportamento proximo as fronteiras.",
      mu = mu, var = variance,
      mu4 = central4(raw(1), raw(2), raw(3), raw(4)),
      r = function(m) rbeta(m, shape1 = alpha, shape2 = beta),
      d = function(x) dbeta(x, shape1 = alpha, shape2 = beta),
      q = function(p) qbeta(p, shape1 = alpha, shape2 = beta)
    ))
  }
  if (dist == "exponencial") {
    lambda <- input$lambda
    req(lambda)
    variance <- 1 / lambda^2
    return(list(
      nome = "Exponencial", notacao = "X\\sim\\operatorname{Exp}(\\lambda)",
      parametros = "\\lambda>0;\\quad X\\ge 0",
      f = "f(x)=\\lambda e^{-\\lambda x},\\quad x\\ge 0",
      aplicacao = "Tempos entre ocorrencias em um processo de Poisson homogeneo e tempos até falha sob taxa de falha constante. O TCL se refere aqui as médias, não a propria distribuição exponencial.",
      mu = 1 / lambda, var = variance, mu4 = 9 / lambda^4,
      r = function(m) rexp(m, rate = lambda),
      d = function(x) dexp(x, rate = lambda),
      q = function(p) qexp(p, rate = lambda)
    ))
  }
  stop("Distribuição populacional não reconhecida.")
}

# Densidade para histogramas e curva teorica na mesma escala.
plot_theme <- function(base = 13) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(face = "bold", colour = "#142D50", size = base + 2),
          plot.subtitle = element_text(colour = "#53657B", size = base - 1),
          axis.title = element_text(colour = "#32465E"),
          axis.text = element_text(colour = "#53657B"),
          panel.grid.major.y = element_line(colour = "#E9EEF5", linewidth = .35),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = base - 2),
          strip.background = element_rect(fill = "#EAF1F7", colour = NA),
          strip.text = element_text(face = "bold", colour = "#142D50"),
          plot.margin = margin(14, 16, 12, 12))
}

metric <- function(label, value, detail) {
  div(class = "metric-card",
      div(class = "metric-label", label),
      div(class = "metric-value", value),
      div(class = "metric-detail", detail))
}

# -----------------------------------------------------------------------------
# 2. Interface: um unico app para os tres experimentos
# -----------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("\n      body {background:#F3F6FA;color:#1B304A;font-family:'Segoe UI',Arial,sans-serif;}\n      .container-fluid {max-width:1680px;padding:0 25px 28px;}\n      .app-header {margin:0 -25px 22px;padding:24px 30px;background:#112D52;\n                   border-bottom:5px solid #D5AE48;color:white;}\n      .app-header h1 {font-size:30px;letter-spacing:-.5px;font-weight:750;margin:0 0 6px;}\n      .app-header p {font-size:14px;color:#D5E2F0;margin:0;}\n      .eyebrow {font-size:11px;font-weight:800;letter-spacing:1.7px;color:#F4CC6E;margin-bottom:8px;}\n      .panel-card {border:1px solid #DCE5EF;border-radius:12px;background:#FFF;\n                   padding:21px;margin-bottom:19px;box-shadow:0 2px 12px rgba(16,42,77,.04);}\n      .sidebar-card {padding:23px;}\n      .section-label {font-size:11px;color:#58708B;text-transform:uppercase;\n                      letter-spacing:1.1px;font-weight:800;margin:8px 0 13px;}\n      .sidebar-card label {font-weight:650;color:#213A56;}\n      .sidebar-card .form-group {margin-bottom:20px;}\n      .sidebar-card .irs--shiny .irs-bar,.irs--shiny .irs-single {background:#1A788F;border-color:#1A788F;}\n      .sidebar-card .btn-primary {background:#173E66;border-color:#173E66;}\n      .hero-note {color:#435A71;font-size:13px;line-height:1.65;}\n      .metric-grid {display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-bottom:18px;}\n      .metric-card {background:#FFF;border:1px solid #DDE5EF;border-radius:11px;padding:16px 17px;\n                    border-top:3px solid #247F98;min-width:0;}\n      .metric-label {font-size:11px;letter-spacing:.6px;text-transform:uppercase;font-weight:800;color:#58708B;}\n      .metric-value {font-size:25px;font-weight:800;letter-spacing:-.6px;color:#143658;overflow-wrap:anywhere;}\n      .metric-detail {font-size:12px;color:#6C8093;margin-top:4px;}\n      .nav-tabs {border-bottom:1px solid #DCE5EF;margin-bottom:15px;}\n      .nav-tabs>li>a {font-weight:650;color:#526C85;border-radius:6px 6px 0 0;}\n      .nav-tabs>li.active>a,.nav-tabs>li.active>a:focus {color:#143D63;border-bottom:3px solid #D5AE48;}\n      .math-block {padding:13px 15px;background:#F2F6FB;border-radius:7px;overflow-x:auto;margin:10px 0 17px;}\n      .explain {font-size:14px;line-height:1.8;color:#31465B;}\n      .explain h4 {color:#143658;font-weight:750;margin-top:22px;}\n      .note {background:#EDF6F8;border-left:4px solid #198198;border-radius:5px;padding:11px 14px;\n             margin:15px 0;line-height:1.6;font-size:13px;}\n      .warn {background:#FFF7E7;border-left-color:#C58B24;}\n      .footer {color:#688097;font-size:12px;text-align:right;padding:6px 0;}\n      @media(max-width:1100px){.metric-grid{grid-template-columns:repeat(2,minmax(0,1fr));}}\n      @media(max-width:600px){.container-fluid{padding:0 12px 18px;}\n       .app-header{margin:0 -12px 14px;padding:20px 16px;}\n       .app-header h1{font-size:24px;}.metric-value{font-size:20px;}\n       .metric-card{padding:12px;}.panel-card{padding:14px;}}\n    "))
  ),
  div(class = "app-header",
      div(class = "eyebrow", "UFLA  /  ICTIN  /  GCT018  /  LABORATORIO INTERATIVO"),
      h1("Distribuições amostrais"),
      p("Média e proporção: TCL  |  Variância: distribuição de S², lei exata e aproximação assintótica")),
  fluidRow(
    column(4,
      div(class = "panel-card sidebar-card",
        div(class = "section-label", "01  ·  Escolha o experimento"),
        selectInput("tipo", "Estatística de interesse:",
                    choices = c("Média amostral (TCL)" = "media",
                                "Proporção amostral (TCL)" = "proporcao",
                                "Variância amostral S² (lei exata / aproximação)" = "variancia"),
                    selected = "media"),
        div(class = "section-label", "02  ·  Modelo populacional"),
        uiOutput("escolha_pop"),
        uiOutput("parametros"),
        div(class = "section-label", "03  ·  Plano de simulação"),
        sliderInput("n", "n — observações por amostra:", min = 2, max = 500,
                    value = 30, step = 1, sep = ""),
        sliderInput("k", "B — amostras independentes (replicações):",
                    min = 200, max = 4000, value = 2000, step = 100, sep = ""),
        numericInput("seed", "Semente reprodutível:", value = 2026, min = 1,
                     max = 1000000, step = 1),
        actionButton("nova", "Nova realização aleatória", class = "btn-primary", width = "100%"),
        br(), br(),
        downloadButton("baixar", "Exportar simulação (.csv)", class = "btn-default", width = "100%"),
        div(class = "section-label", "04  ·  Visualização"),
        uiOutput("escala_ui"),
        sliderInput("bins", "Classes do histograma (quando aplicável):",
                    min = 12, max = 80, value = 35, step = 1),
        div(class = "hero-note",
            "n altera a distribuição teórica da estatística. B altera apenas a precisão da aproximação por Monte Carlo. Os gráficos e os indicadores usam as mesmas amostras.")
      )
    ),
    column(8,
      uiOutput("indicadores"),
      div(class = "panel-card",
        tabsetPanel(id = "aba", type = "tabs",
          tabPanel("Distribuição amostral",
            plotOutput("grafico_amostral", height = "435px"),
            uiOutput("legenda_resultado")),
          tabPanel("População e seis amostras",
            plotOutput("grafico_pop", height = "290px"),
            tags$hr(),
            plotOutput("grafico_amostras", height = "440px"),
            div(class = "hero-note",
                "As seis amostras são as primeiras seis replicações usadas no gráfico amostral; seus resumos são calculados com os mesmos dados.")),
          tabPanel("Fundamentação e interpretação",
            uiOutput("fundamentacao"))
        )
      ),
      div(class = "footer", "GCT018 · Estatística: conceitos e experimentações · UFLA / ICTIN")
    )
  )
)

# -----------------------------------------------------------------------------
# 3. Servidor: parametros corretos, simulacao unica, graficos coerentes
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  output$escolha_pop <- renderUI({
    if (input$tipo == "proporcao") {
      tagList(
        div(class = "hero-note", "População de Bernoulli: cada unidade vale 1 (sucesso) ou 0 (fracasso)."),
        sliderInput("p", "p — probabilidade populacional de sucesso:",
                    min = 0, max = 1, value = .5, step = .01)
      )
    } else {
      selectInput("dist", "Distribuição dos dados individuais:",
                  choices = c("Normal" = "normal", "Uniforme contínua" = "uniforme",
                              "Lognormal (cauda a direita)" = "lognormal",
                              "Beta (assimetria configuravel)" = "beta",
                              "Exponencial" = "exponencial"),
                  selected = "normal")
    }
  })

  output$parametros <- renderUI({
    if (input$tipo == "proporcao") return(NULL)
    req(input$dist)
    switch(input$dist,
      normal = tagList(
        sliderInput("mu", "mu — média populacional:", min = -20, max = 20,
                    value = 0, step = .5),
        sliderInput("sigma", "sigma — desvio-padrão populacional:",
                    min = .2, max = 15, value = 2, step = .1)),
      uniforme = sliderInput("ab", "a e b — limites (a < b):",
                            min = -10, max = 20, value = c(0, 10), step = .5),
      lognormal = tagList(
        sliderInput("mul", "mu_log — média de log(X):",
                    min = -1.5, max = 2, value = 0, step = .1),
        sliderInput("sigmal", "sigma_log — DP de log(X):",
                    min = .1, max = 1.5, value = .6, step = .05)),
      beta = tagList(
        sliderInput("alpha", "alpha — parâmetro de forma:",
                    min = .3, max = 12, value = 2, step = .1),
        sliderInput("beta", "beta — parâmetro de forma:",
                    min = .3, max = 12, value = 5, step = .1)),
      exponencial = sliderInput("lambda", "lambda — taxa (> 0):",
                                min = .1, max = 5, value = 1, step = .1)
    )
  })

  output$escala_ui <- renderUI({
    if (input$tipo == "media") {
      radioButtons("escala_media", "Escala da distribuição:",
                   choices = c("Média original (X-barra)" = "natural",
                               "Média padronizada (Z)" = "padrao"), selected = "natural")
    } else if (input$tipo == "proporcao") {
      radioButtons("escala_prop", "Escala da proporção:",
                   choices = c("Proporção original (p-chapeu)" = "natural",
                               "Estatística padronizada (Z)" = "padrao"),
                   selected = "natural")
    } else {
      radioButtons("escala_var", "O que deseja estudar?",
                   choices = c("PRINCIPAL: distribuição de S²" = "natural",
                               "COMPLEMENTAR: Q = (n-1)S²/σ²" = "padrao",
                               "COMPLEMENTAR: estatística centralizada T" = "z"),
                   selected = "natural")
    }
  })

  modelo <- reactive({
    if (input$tipo == "proporcao") {
      req(input$p)
      p <- input$p
      return(list(nome = "Bernoulli", notacao = "X_i\\sim\\operatorname{Bernoulli}(p)",
        parametros = "0\\le p\\le1;\\quad X_i\\in\\{0,1\\}",
        f = "P(X_i=x)=p^x(1-p)^{1-x},\\quad x\\in\\{0,1\\}",
        aplicacao = "Indicadores de conformidade, ocorrencia de defeito por item e respostas binarias. A proporção amostral e a média de variaveis de Bernoulli independentes.",
        mu = p, var = p * (1 - p), mu4 = NA_real_,
        r = function(m) rbinom(m, 1, p)))
    }
    req(input$dist)
    make_model(input$dist, input)
  })

  # Os sliders de aparencia (bins/escala/aba) NAO fazem novas simulacoes.
  # Uma alteracao de n, B, parametros ou botao "Nova realização" faz.
  simulacao <- reactive({
    m <- modelo()
    n <- as.integer(input$n)
    B <- as.integer(input$k)
    validate(need(n >= 2 && B >= 6, "Utilize n >= 2 e B >= 6."))
    req(input$seed)
    seed <- as.integer(input$seed) + as.integer(input$nova %||% 0L)
    set.seed(seed)

    if (input$tipo == "proporcao") {
      # Exatamente Bin(n,p)/n: evita subamostrar uma populacao empirica finita.
      seis <- matrix(rbinom(6L * n, size = 1L, prob = m$mu), nrow = n,
                     ncol = 6L)
      primeiros <- colMeans(seis)
      restantes <- rbinom(B - 6L, size = n, prob = m$mu) / n
      stat <- c(primeiros, restantes)
      return(list(stat = stat, seis = seis, n = n, B = B))
    }

    # B amostras i.i.d. da distribuicao analitica, sem populacao artificial.
    x <- matrix(m$r(n * B), nrow = n, ncol = B)
    medias <- colMeans(x)
    if (input$tipo == "media") {
      stat <- medias
    } else {
      # Estatistica CORRETA: denominador n-1 (e nao n).
      # Centralizacao por coluna para evitar cancelamento numerico.
      residuos <- x - rep(medias, each = n)
      stat <- colSums(residuos^2) / (n - 1)
    }
    list(stat = stat, seis = x[, seq_len(6L), drop = FALSE], n = n, B = B)
  })

  # Media/DP teoricos para a estatistica no painel, nao para X individual.
  referencia <- reactive({
    m <- modelo(); s <- simulacao(); n <- s$n
    tipo <- input$tipo
    if (tipo == "media") {
      escala <- input$escala_media %||% "natural"
      if (escala == "padrao") {
        return(list(x = (s$stat - m$mu) / sqrt(m$var / n),
          mean = 0, sd = 1, modo = "contínua", dens = function(x) dnorm(x),
          quant = function(p) qnorm(p),
          titulo = "Média amostral padronizada", xlabel = "Z = (X-barra - mu)/(sigma/raiz(n))",
          legenda = if (input$dist == "normal") "Normal padrão exata" else "Normal padrão: aproximação TCL"))
      }
      se <- sqrt(m$var / n)
      return(list(x = s$stat, mean = m$mu, sd = se,
        modo = "contínua", dens = function(x) dnorm(x, m$mu, se),
        quant = function(p) qnorm(p, m$mu, se),
        titulo = "Distribuição das médias amostrais",
        xlabel = "Média amostral (X-barra)",
        legenda = if (input$dist == "normal") "Normal exata" else "Normal: aproximação TCL"))
    }
    if (tipo == "proporcao") {
      p <- m$mu; se <- sqrt(m$var / n)
      escala <- input$escala_prop %||% "natural"
      standardized <- escala == "padrao" && se > 0
      x <- if (standardized) (s$stat - p) / se else s$stat
      return(list(x = x, mean = if (standardized) 0 else p,
        sd = if (standardized) 1 else se, modo = "discreta",
        standardized = standardized,
        titulo = "Distribuição das proporções amostrais",
        xlabel = if (standardized) "Z = (p-chapeu - p) / EP(p-chapeu)" else "Proporção amostral (p-chapeu)",
        legenda = "Probabilidade exata: Binomial"))
    }
    sigma2 <- m$var; df <- n - 1
    # A grandeza de interesse didatico e S²; Q e T sao transformacoes opcionais.
    escala <- input$escala_var %||% "natural"
    if (escala == "z") {
      # T transforma uma lei qui-quadrado exata em sua versao centrada.
      # Para nao normais, T usa o EP assintotico dependente de mu_4.
      normal <- identical(input$dist, "normal")
      if (normal) {
        fator <- sqrt(2 * df)
        return(list(x = (df * s$stat / sigma2 - df) / fator,
          mean = 0, sd = 1, modo = "contínua",
          dens = function(t) dchisq(df + fator * t, df) * fator,
          quant = function(p) (qchisq(p, df) - df) / fator,
          dens2 = function(t) dnorm(t),
          legenda2 = "Normal limite (aproximação)",
          titulo = "Qui-quadrado centralizada versus limite normal",
          xlabel = "T = (Q - (n-1)) / raiz(2(n-1))",
          legenda = "Qui-quadrado centralizada: exata"))
      }
      se_as <- sqrt((m$mu4 - sigma2^2) / n)
      se_ex <- sqrt((m$mu4 - ((n - 3)/(n - 1)) * sigma2^2) / n)
      return(list(x = (s$stat - sigma2) / se_as,
        mean = 0, sd = se_ex / se_as, modo = "contínua",
        dens = function(t) dnorm(t), quant = function(p) qnorm(p),
        titulo = "Variância centralizada: limite normal assintótico",
        xlabel = "T = (S^2 - sigma^2) / raiz[(mu_4 - sigma^4)/n]",
        legenda = "Normal limite: aproximação (população não normal)"))
    }
    if (escala == "padrao") {
      normal <- identical(input$dist, "normal")
      return(list(x = df * s$stat / sigma2,
        mean = df, sd = if (normal) sqrt(2 * df) else
          df / sigma2 * sqrt((m$mu4 - ((n - 3) / (n - 1)) * sigma2^2) / n),
        modo = "contínua", dens = if (normal) function(x) dchisq(x, df) else NULL,
        quant = if (normal) function(p) qchisq(p, df) else NULL,
        titulo = "Variância padronizada: Q",
        xlabel = "Q = (n - 1) S^2 / sigma^2",
        legenda = if (normal) "Qui-quadrado exata" else "Sem curva qui-quadrado: população não normal"))
    }
    # Formula exata para Var(S^2), com quarto momento central finito.
    sd_s2 <- sqrt((m$mu4 - ((n - 3) / (n - 1)) * sigma2^2) / n)
    normal <- identical(input$dist, "normal")
    asy_sd <- sqrt(pmax(0, (m$mu4 - sigma2^2) / n))
    return(list(x = s$stat, mean = sigma2, sd = sd_s2,
      modo = "contínua",
      dens = if (normal) function(x) dchisq(x * df / sigma2, df) * df / sigma2
             else function(x) dnorm(x, sigma2, asy_sd),
      quant = if (normal) function(p) qchisq(p, df) * sigma2 / df else NULL,
      titulo = "Distribuição amostral de S²",
      xlabel = "Variância amostral S² (unidades de X ao quadrado)",
      legenda = if (normal) "Densidade de S²: qui-quadrado reescalada (exata)" else "Densidade normal assintótica de S² (aproximação)"))
  })

  output$indicadores <- renderUI({
    s <- simulacao(); ref <- referencia(); m <- modelo()
    raw_mean <- if (input$tipo == "variancia") m$var else m$mu
    raw_sd <- if (input$tipo == "media") sqrt(m$var / s$n) else if (input$tipo == "proporcao")
      sqrt(m$var / s$n) else sqrt((m$mu4 - ((s$n - 3)/(s$n - 1)) * m$var^2) / s$n)
    list_vals <- list(
      metric(if (input$tipo == "variancia") "Variância populacional σ²" else "Centro teórico da estatística",
             fmt(raw_mean), if (input$tipo == "variancia") "E(S²) = σ² (valor de referência)" else "E[estatística] na escala original"),
      metric(if (input$tipo == "variancia") "DP teórico de S²" else "Desvio-padrão teórico",
             fmt(raw_sd), if (input$tipo == "variancia") "Dispersão entre variâncias amostrais" else "Erro-padrão ou DP da estatística"),
      metric("Média das B replicações", fmt(mean(s$stat)), paste0("B = ", s$B)),
      metric("DP das B replicações", fmt(sd(s$stat)), "Variabilidade de Monte Carlo"))
    div(class = "metric-grid", list_vals)
  })

  output$grafico_amostral <- renderPlot({
    m <- modelo(); s <- simulacao(); r <- referencia()
    if (r$modo == "discreta") {
      n <- s$n; p <- m$mu; se <- sqrt(m$var / n)
      ks <- 0:n
      obs <- tabulate(round(s$stat * n) + 1L, nbins = n + 1L) / s$B
      xk <- if (r$standardized) (ks/n - p)/se else ks/n
      dd <- data.frame(x = xk, obs = obs, exact = dbinom(ks, n, p))
      if (se > 0) {
        # Correcao de continuidade nas CELULAS k-0,5 a k+0,5.
        dd$approx <- pnorm((ks + .5)/n, mean = p, sd = se) -
                     pnorm((ks - .5)/n, mean = p, sd = se)
      }
      delta <- if (r$standardized) 1/(n*se) else 1/n
      # Barras e curvas estao todas em PROBABILIDADE POR VALOR, nao em densidade.
      g <- ggplot(dd, aes(x = x)) +
        geom_col(aes(y = obs, fill = "Monte Carlo"), width = .76 * delta,
                 colour = "white", linewidth = .14) +
        geom_line(aes(y = exact, colour = "Binomial exata"), linewidth = 1.1) +
        geom_point(aes(y = exact, colour = "Binomial exata"), size = .9) +
        labs(title = r$titulo,
             subtitle = paste0(s$B, " replicações independentes de tamanho ", n,
                               "  ·  barras e curvas em massa de probabilidade"),
             x = r$xlabel, y = "Probabilidade de cada valor") +
        scale_fill_manual(values = c("Monte Carlo" = "#71BDD0")) +
        scale_colour_manual(values = c("Binomial exata" = "#143B63",
                                       "Normal aprox. (correcao)" = "#BC772D")) +
        plot_theme()
      if (se > 0) {
        g <- g + geom_line(aes(y = approx, colour = "Normal aprox. (correcao)"),
                           linewidth = 1, linetype = "dashed")
      }
      low <- max(0, floor(min(c(qbinom(.001,n,p), round(s$stat*n))) - 1))
      high <- min(n, ceiling(max(c(qbinom(.999,n,p), round(s$stat*n))) + 1))
      if (low == high) {low <- max(0, low - 1); high <- min(n, high + 1)}
      g + coord_cartesian(xlim = range(xk[(low + 1):(high + 1)]),
                          ylim = c(0, max(dd$obs, dd$exact, dd$approx, na.rm = TRUE)*1.12),
                          expand = FALSE) +
        guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2))
    } else {
      vals <- r$x
      rng <- as.numeric(quantile(vals, c(.001, .999), names = FALSE))
      if (!is.null(r$quant)) {
        qtheory <- r$quant(c(.002, .998))
        rng <- range(rng, qtheory[is.finite(qtheory)])
      }
      nonnegative_variance_axis <- input$tipo == "variancia" &&
        !identical(input$escala_var, "z")
      # Inclua sempre o centro populacional no campo visivel do grafico.
      # Na escala natural, este valor e a variancia populacional sigma².
      rng <- range(rng, r$mean)
      if (nonnegative_variance_axis) rng[1] <- max(0, rng[1])
      if (diff(rng) < .Machine$double.eps^(.5)) {
        rng <- c(rng[1] - .5, rng[2] + .5)
      }
      xr <- seq(max(rng[1], if (nonnegative_variance_axis) .Machine$double.eps else -Inf),
                rng[2], length.out = 800)
      curve_data <- if (!is.null(r$dens))
        data.frame(x = xr, y = r$dens(xr)) else NULL
      curve_data2 <- if (!is.null(r$dens2))
        data.frame(x = xr, y = r$dens2(xr)) else NULL
      # Binagem SOMENTE na faixa exibida para que observacoes extremas nao
      # comprimam as barras. Denominador B (e nao numero de observacoes
      # visiveis): a densidade NAO e renormalizada ao recortar as caudas.
      n_classes <- as.integer(input$bins)
      bordas <- seq(rng[1], rng[2], length.out = n_classes + 1L)
      largura <- diff(rng) / n_classes
      classes <- findInterval(vals, bordas, rightmost.closed = TRUE)
      contagens <- tabulate(classes, nbins = n_classes)
      hist_dados <- data.frame(
        centro = head(bordas, -1L) + largura / 2,
        densidade = contagens / (length(vals) * largura))
      ylim <- max(hist_dados$densidade, 0.001) * 1.30
      if (!is.null(curve_data)) {
        finite_y <- curve_data$y[is.finite(curve_data$y)]
        ylim <- max(ylim, as.numeric(quantile(finite_y, .97)) * 1.15)
      }
      if (!is.null(curve_data2)) {
        ylim <- max(ylim, max(curve_data2$y[is.finite(curve_data2$y)]) * 1.15)
      }
      g <- ggplot() +
        geom_col(data = hist_dados,
                 aes(x = centro, y = densidade, fill = "Monte Carlo"),
                 width = largura, colour = "white", linewidth = .28) +
        geom_vline(xintercept = r$mean, colour = "#C38A2E",
                   linetype = "dotted", linewidth = .9) +
        scale_fill_manual(values = c("Monte Carlo" = "#78C0D0")) +
        labs(title = r$titulo,
             subtitle = paste0(s$B, " replicações independentes, n = ", s$n,
               if (input$tipo == "variancia" &&
                   identical(input$escala_var %||% "natural", "natural"))
                 paste0("  ·  linha pontilhada: σ² = ", fmt(m$var))
               else "  ·  linha pontilhada: centro teórico"),
             x = r$xlabel, y = "Densidade") + plot_theme()
      if (!is.null(curve_data)) {
        g <- g + geom_line(data = curve_data, aes(x = x, y = y, colour = r$legenda),
                           inherit.aes = FALSE, linewidth = 1.15)
        cores <- setNames("#163C67", r$legenda)
        if (!is.null(curve_data2)) {
          g <- g + geom_line(data = curve_data2,
                 aes(x = x, y = y, colour = r$legenda2),
                 inherit.aes = FALSE, linewidth = 1, linetype = "dashed")
          cores <- c(cores, setNames("#BC772D", r$legenda2))
        }
        g <- g + scale_colour_manual(values = cores)
      }
      g + coord_cartesian(xlim = rng, ylim = c(0, ylim), expand = FALSE)
    }
  }, res = 110)

  output$legenda_resultado <- renderUI({
    m <- modelo(); s <- simulacao(); r <- referencia()
    if (input$tipo == "media") {
      text <- if (input$dist == "normal")
        "A curva normal e EXATA para qualquer n porque os dados individuais sao normais. Barras irregulares refletem apenas o número finito B de replicações."
      else "A curva normal e ASSINTOTICA (TCL). Aumentar n altera a forma da distribuição das médias; aumentar B reduz o ruído visual do histograma, sem alterar a distribuição teórica."
    } else if (input$tipo == "proporcao") {
      p <- m$mu
      text <- if (p == 0 || p == 1)
        "Caso degenerado: todas as proporções amostrais sao identicas. O erro-padrão e zero e a padronizacao Z não esta definida; a escala original foi utilizada."
      else paste0("A distribuição exata e discreta: n*p = ", fmt(s$n*p, 1),
                  "; n*(1-p) = ", fmt(s$n*(1-p), 1),
                  ". A linha tracejada usa a aproximação normal com correcao de continuidade. ",
                  if (min(s$n*p, s$n*(1-p)) < 10) "Atencao: a condição pedagogica np e n(1-p) >= 10 NÃO e atendida; a aproximação pode ser inadequada."
                  else "As contagens esperadas sao >= 10; compare a aproximação com a curva binomial exata.")
    } else if (input$dist == "normal") {
      text <- if (identical(input$escala_var %||% "natural", "natural"))
        paste0("GRÁFICO PRINCIPAL — As barras mostram as ", s$B,
               " variâncias amostrais S², com denominador n-1. A curva azul é a densidade EXATA de uma qui-quadrado REESCALADA (não a qui-quadrado padrão). A linha pontilhada é σ² = ",
               fmt(m$var), ". Ao aumentar n, S² se concentra em torno de σ²; aumentar B somente suaviza a estimativa por Monte Carlo.")
      else if (identical(input$escala_var, "z"))
        "VISUALIZAÇÃO COMPLEMENTAR — A curva azul é a lei EXATA de (Q-df)/raiz(2df); a curva tracejada é a normal padrão LIMITE. A aproximação melhora com df=n-1 crescentes; nenhuma dessas curvas é ajustada ao histograma."
      else "VISUALIZAÇÃO COMPLEMENTAR — Somente Q=(n-1)S²/σ² tem distribuição qui-quadrado PADRÃO exata sob normalidade. Na visualização principal, S² possui essa lei REESCALADA."
    } else {
      text <- if (identical(input$escala_var %||% "natural", "natural"))
        paste0("GRÁFICO PRINCIPAL — As barras são as ", s$B,
               " variâncias amostrais S². A linha pontilhada é σ² = ", fmt(m$var),
               ". Sem normalidade populacional, não há lei qui-quadrado exata em geral. A curva azul é apenas uma aproximação NORMAL ASSINTÓTICA para S²; para n pequeno ou forte assimetria ela pode ser inadequada e atribuir densidade a valores negativos, impossíveis para S².")
      else if (identical(input$escala_var, "padrao"))
        "VISUALIZAÇÃO COMPLEMENTAR — Para população não normal, Q NÃO segue, em geral, uma qui-quadrado. Nenhuma curva qui-quadrado foi sobreposta; esta transformação é exibida somente para contraste conceitual."
      else "VISUALIZAÇÃO COMPLEMENTAR — T utiliza o quarto momento central da população, NÃO a padronização qui-quadrado. A referência N(0,1) é ASSINTÓTICA e depende de n suficientemente grande."
    }
    div(class = paste("note", if (input$tipo == "variancia" && input$dist != "normal" ||
                                     input$tipo == "proporcao" && min(s$n*m$mu, s$n*(1-m$mu)) < 10)
                       "warn" else ""), text)
  })

  output$grafico_pop <- renderPlot({
    m <- modelo()
    if (input$tipo == "proporcao") {
      dd <- data.frame(resultado = factor(c(0,1)), prob = c(1-m$mu,m$mu))
      ggplot(dd, aes(x = resultado, y = prob)) +
               geom_col(width = .52, fill = "#357F99") +
               geom_text(aes(label = sprintf("%.3f",prob)), vjust = -.6,
                         colour = "#143658", fontface = "bold") +
               scale_y_continuous(limits = c(0, 1.13), expand = c(0,0)) +
               labs(title = "População Bernoulli", subtitle = "Distribuição EXATA dos dados individuais",
                    x = "Resultado X_i", y = "Probabilidade") + plot_theme()
    } else {
    x <- if (identical(input$dist, "beta")) seq(.001, .999, length.out = 1100)
         else seq(m$q(.001), m$q(.999), length.out = 1100)
    y <- m$d(x)
    dd <- data.frame(x = x, dens = y)
    ggplot(dd, aes(x = x, y = dens)) +
      geom_line(colour = "#287F97", linewidth = 1.3) +
      geom_vline(xintercept = m$mu, linetype = "dashed", colour = "#C08A36") +
      labs(title = paste("População:", m$nome),
           subtitle = "Densidade populacional analítica; tracejado = média populacional",
           x = "Observacao individual X", y = "Densidade") + plot_theme()
    }
  }, res = 110)

  output$grafico_amostras <- renderPlot({
    m <- modelo(); s <- simulacao()
    dados <- data.frame(amostra = factor(rep(seq_len(6L), each = s$n),
                                           levels = seq_len(6L)),
                        x = as.vector(s$seis))
    estat <- if (input$tipo == "proporcao") colMeans(s$seis) else
      if (input$tipo == "media") colMeans(s$seis) else apply(s$seis, 2, var)
    nome <- switch(input$tipo, media = "media", proporcao = "proporcao", variancia = "s²")
    levels(dados$amostra) <- paste0("Amostra ", 1:6, "  |  ", nome, " = ",
                                    vapply(estat, fmt, character(1), digits = 3))
    if (input$tipo == "proporcao") {
      d <- aggregate(rep(1, nrow(dados)), by = list(amostra = dados$amostra,
                                                     x = dados$x), FUN = sum)
      names(d)[3] <- "contagem"
      d$freq <- d$contagem / s$n
      ggplot(d, aes(x = factor(x, levels = c(0,1)), y = freq)) +
        geom_col(fill = "#438FA8", width = .5) +
        facet_wrap(~amostra, ncol = 3) +
        scale_y_continuous(limits = c(0,1.05), expand = c(0,0)) +
        labs(title = "Seis amostras Bernoulli", x = "Valor individual", y = "Frequência relativa") +
        plot_theme(10)
    } else {
    means <- aggregate(x ~ amostra, data = dados, FUN = mean)
    ggplot(dados, aes(x = x)) +
      geom_histogram(aes(y = after_stat(density)), bins = 16,
                     fill = "#7ABCCA", colour = "white") +
      geom_vline(data = means, aes(xintercept = x), linewidth = .7,
                 colour = "#BA872B", linetype = "dashed") +
      facet_wrap(~amostra, ncol = 3, scales = "free_y") +
      labs(title = "Seis amostras individuais",
           subtitle = "Cada faceta: n observações; tracejado = média DAQUELA amostra",
           x = "Observacao individual", y = "Densidade") +
      plot_theme(10)
    }
  }, res = 105)

  output$fundamentacao <- renderUI({
    m <- modelo(); s <- simulacao(); n <- s$n
    if (input$tipo == "media") {
      exato <- identical(input$dist, "normal")
      conteudo <- tagList(
        h4("1. Definicao e pressupostos"),
        p("Considere X_1,...,X_n independentes e identicamente distribuidas, com média finita mu e variância populacional sigma² finita e positiva. A média de cada amostra e:"),
        div(class = "math-block", math("\\bar X=\\frac{1}{n}\\sum_{i=1}^n X_i", TRUE)),
        h4("2. Média, variância e erro-padrão — identidades exatas"),
        div(class = "math-block", math("E(\\bar X)=\\mu,\\qquad \\operatorname{Var}(\\bar X)=\\frac{\\sigma^2}{n},\\qquad EP(\\bar X)=\\frac{\\sigma}{\\sqrt n}",TRUE)),
        h4("3. Teorema Central do Limite"),
        div(class = "math-block", math("\\frac{\\sqrt n(\\bar X-\\mu)}{\\sigma}\\xrightarrow{d}N(0,1),\\qquad \\bar X\\;\\dot\\sim\\;N\\!\\left(\\mu,\\frac{\\sigma^2}{n}\\right)",TRUE)),
        div(class = if (exato) "note" else "note warn",
            if (exato) "EXATIDAO: com população normal, a distribuição de X-barra e normal para TODO n; não e necessário invocar o TCL para essa conclusão."
            else "ASSINTOTICA: a aproximação normal depende de n e da forma populacional; n = 30 não e um limiar universal. Populacoes muito assimétricas podem exigir amostras maiores."),
        h4("4. População selecionada"),
        p("Modelo: ", math(m$notacao)),
        div(class = "math-block", math(m$f, TRUE)),
        p("Dominio: ", math(m$parametros)),
        p(m$aplicacao),
        h4("5. Como interpretar os graficos"),
        p("O primeiro grafico mostra os dados individuais; o segundo, seis amostras com n observações cada; o grafico amostral agrega as B médias obtidas. A curva teórica vem dos momentos populacionais EXATOS, nunca do ajuste de uma densidade aos dados simulados."),
        p("Ao aumentar n, o erro-padrão cai com 1/raiz(n). Ao aumentar apenas B, reduz-se o ruído de Monte Carlo, sem alterar a lei teórica de X-barra.")
      )
    } else if (input$tipo == "proporcao") {
      p <- m$mu; se <- sqrt(p * (1-p) / n)
      conteudo <- tagList(
        h4("1. Definicao e distribuição EXATA"),
        p("Considere n respostas independentes de Bernoulli(p), em que 1 = sucesso e 0 = fracasso. K = soma(X_i) conta os sucessos, e p-chapeu = K/n."),
        div(class = "math-block", math("K=\\sum_{i=1}^n X_i\\sim\\operatorname{Bin}(n,p),\\qquad \\hat p=\\frac{K}{n}",TRUE)),
        div(class = "math-block", math("P\\!\\left(\\hat p=\\frac{k}{n}\\right)={n\\choose k}p^k(1-p)^{n-k},\\quad k=0,\\ldots,n",TRUE)),
        h4("2. Momentos e aproximação normal"),
        div(class = "math-block", math("E(\\hat p)=p,\\qquad \\operatorname{Var}(\\hat p)=\\frac{p(1-p)}{n},\\qquad EP(\\hat p)=\\sqrt{\\frac{p(1-p)}{n}}",TRUE)),
        div(class = "math-block", math("Z=\\frac{\\hat p-p}{\\sqrt{p(1-p)/n}}\\xrightarrow{d}N(0,1),\\qquad 0<p<1",TRUE)),
        div(class = "note", paste0("Condição pratica conservadora para inspecao visual: np >= 10 e n(1-p) >= 10. Aqui: np = ",
                                   fmt(n*p,2), "; n(1-p) = ",fmt(n*(1-p),2), ". Trata-se de heuristica, não de teorema nem garantia de erro pequeno.")),
        h4("3. Por que o grafico mostra barras?"),
        p("A proporção amostral possui suporte discreto {0,1/n,...,1}. As barras mostram frequências relativas simuladas, a linha contínua representa a massa BINOMIAL EXATA e a tracejada e a aproximação normal integrada em celulas de largura 1/n (correcao de continuidade). Todas estao na mesma unidade: probabilidade por valor possível."),
        p("Em p = 0 ou p = 1, a distribuição e degenerada, o erro-padrão vale zero e Z não esta definido. A visualizacao permanece na escala original."),
        h4("4. Modelo e aplicação"),
        p("Modelo: ", math(m$notacao)),
        div(class = "math-block", math(m$f,TRUE)),
        p("Dominio: ",math(m$parametros)),
        p(m$aplicacao),
        h4("5. Diferenca entre n e B"),
        p("n define a dispersao de p-chapeu; B e somente o número de repeticoes artificiais. Aumentar B torna as barras mais proximas das probabilidades binomiais exatas, mas não melhora, por si só, a aproximação normal.")
      )
    } else {
      normal <- identical(input$dist,"normal")
      conteudo <- tagList(
        h4("1. Gráfico principal: a distribuição amostral de S²"),
        p("Geramos B amostras independentes, cada uma com n observações. O histograma reúne uma variância S²_j por amostra, usando o denominador n-1. A linha vertical indica a variância populacional σ², que é o valor esperado de S²."),
        div(class = "math-block", math("S_j^2=\\frac{1}{n-1}\\sum_{i=1}^n(X_{ij}-\\bar X_j)^2,\\quad j=1,\\ldots,B;\\qquad E(S^2)=\\sigma^2",TRUE)),
        p("Conforme n aumenta, o desvio-padrão de S² diminui e a distribuição se concentra em torno de σ². Aumentar somente B deixa o histograma menos irregular, mas não altera a distribuição teórica."),
        h4("2. Curva teórica exata de S²: a população deve ser normal"),
        div(class = "math-block", math("X_i\\overset{iid}{\\sim}N(\\mu,\\sigma^2)\\implies Q=\\frac{(n-1)S^2}{\\sigma^2}\\sim\\chi^2_{n-1}", TRUE)),
        p("A relação acima é EXATA para todo n >= 2 sob normalidade. A curva superposta ao histograma principal é a densidade da qui-quadrado REESCALADA, com mudança de variável e respectivo jacobiano:"),
        div(class = "math-block", math("S^2\\overset d=\\frac{\\sigma^2}{n-1}\\chi^2_{n-1},\\qquad f_{S^2}(s)=\\frac{n-1}{\\sigma^2}f_{\\chi^2_{n-1}}\\!\\left(\\frac{(n-1)s}{\\sigma^2}\\right),\\ s>0",TRUE)),
        div(class = "math-block", math("E(S^2)=\\sigma^2,\\qquad \\operatorname{Var}(S^2)=\\frac{2\\sigma^4}{n-1},\\qquad DP(S^2)=\\sigma^2\\sqrt{\\frac{2}{n-1}}\\quad(\\text{normal})",TRUE)),
        p("A transformação Q constitui uma segunda visualização: ela segue a qui-quadrado padrão. S² em escala original não deve ser comparada diretamente com dchisq(s, df=n-1)."),
        div(class = if (normal) "note" else "note warn",
            if (normal) "A população selecionada e normal: a curva qui-quadrado (ou sua versao reescalada) e referencia teórica exata, mesmo para amostras pequenas."
            else "A população selecionada NÃO e normal: a lei qui-quadrado NÃO se aplica em geral. O grafico de Q não apresenta uma curva qui-quadrado; na escala S², apresenta somente uma aproximação normal assintótica."),
        h4("3. Variância de S² fora da normalidade"),
        p("Se o quarto momento central mu_4 = E[(X-mu)^4] for finito, a variância de S² possui a formula EXATA:"),
        div(class = "math-block", math("\\operatorname{Var}(S^2)=\\frac{1}{n}\\left[\\mu_4-\\frac{n-3}{n-1}\\sigma^4\\right]",TRUE)),
        p("Aproximação para n grande, sob condições usuais de independencia e quarto momento finito:"),
        div(class = "math-block", math("\\sqrt n(S^2-\\sigma^2)\\xrightarrow{d}N(0,\\mu_4-\\sigma^4)",TRUE)),
        div(class = "math-block", math("X_i\\sim N(\\mu,\\sigma^2):\\quad\\frac{Q-(n-1)}{\\sqrt{2(n-1)}}\\xrightarrow{d}N(0,1)",TRUE)),
        p("Essa aproximação e em geral NORMAL para S² centralizada e escalonada; NÃO afirma que Q segue qui-quadrado para uma população não normal. No modo T, a qui-quadrado centralizada e comparada a sua normal limite quando a população e normal."),
        h4("4. Modelo populacional selecionado"),
        p("Modelo: ",math(m$notacao)),
        div(class = "math-block",math(m$f,TRUE)),
        p("Dominio: ",math(m$parametros)),
        p(m$aplicacao),
        h4("5. O que muda quando n aumenta?"),
        p("Na escala S², as variâncias amostrais se concentram ao redor de σ²: isso ilustra a consistência do estimador, não a convergência para uma qui-quadrado fixa. Na escala Q, o número de graus de liberdade n-1 muda com n. Na escala T, uma transformação centralizada e padronizada permite comparar a distribuição com seu limite normal. Sem normalidade, o quarto momento populacional altera a dispersão de S²; não se deve impor a qui-quadrado se suas hipóteses falham.")
      )
    }
    withMathJax(div(class = "explain", conteudo,
      tags$hr(), h4("Referencias e reproducibilidade"),
      p("Casella, G.; Berger, R. L. (2010). Inferencia estatística. Cengage Learning."),
      p("Montgomery, D. C.; Runger, G. C. Applied Statistics and Probability for Engineers. Wiley."),
      p("Bartlett, P. Theoretical Statistics, Lecture 7: U-statistics and asymptotic normality of sample variance. University of California, Berkeley. ",
        tags$a(href = "https://www.stat.berkeley.edu/~bartlett/courses/2013spring-stat210b/notes/7notes.pdf",
               target = "_blank", "Notas em PDF")),
      p("NIST/SEMATECH. e-Handbook of Statistical Methods, Section 7.2.4: Normal approximation for proportions. ",
        tags$a(href = "https://itl.nist.gov/div898/handbook/prc/section2/prc24.htm",
               target = "_blank", "Manual NIST")),
      p("R Core Team. R: A Language and Environment for Statistical Computing. Funções stats: rnorm, rbinom, dbinom, dchisq, qchisq e var."),
      p("Todas as estatísticas sao calculadas a partir das amostras geradas; as curvas e os momentos teóricos sao calculados independentemente da simulação. A semente no painel permite repetir uma realização.")))
  })

  output$baixar <- downloadHandler(
    filename = function() paste0("GCT018_simulacao_", input$tipo, "_n",input$n,
                                  "_B", input$k, ".csv"),
    content = function(file) {
      s <- simulacao(); m <- modelo()
      x <- s$stat
      if (input$tipo == "variancia") {
        dd <- data.frame(replicacao = seq_along(x),
                         variancia_amostral = x,
                         Q = (s$n-1)*x/m$var)
      } else if (input$tipo == "media") {
        dd <- data.frame(replicacao = seq_along(x), media_amostral = x,
                         Z = (x-m$mu)/sqrt(m$var/s$n))
      } else {
        dd <- data.frame(replicacao = seq_along(x),
                         sucessos = round(x*s$n), proporcao_amostral = x,
                         Z = if (m$var > 0) (x-m$mu)/sqrt(m$var/s$n) else NA_real_)
      }
      write.csv2(dd, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
}

shinyApp(ui = ui, server = server)
