# Batting Average as the Product of Three Rates — interactive companion
# Illustrates Albert (2026): BA = c * { h + (1 - h) * b }
#   c = contact rate (AB - SO) / AB
#   h = home-run rate HR / (AB - SO)
#   b = BABIP (H - HR) / (AB - SO - HR)
#
# Packages: shiny, bslib, ggplot2, dplyr, MASS, Lahman (version 14.0-0 or later for 2025 data)

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(MASS, include.only = "mvrnorm")
library(Lahman)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

LAST_YEAR <- max(Batting$yearID)
MIN_AB    <- 300

rate_labels <- c(c = "Contact rate", h = "Home-run rate", b = "BABIP", BA = "Batting average")

bat <- Batting |>
  filter(lgID %in% c("AL", "NL", "FL"), !is.na(SO)) |>
  group_by(playerID, yearID) |>
  summarise(AB = sum(AB), H = sum(H), SO = sum(SO), HR = sum(HR), .groups = "drop") |>
  filter(AB >= MIN_AB, AB > SO) |>
  mutate(
    c  = (AB - SO) / AB,
    h  = HR / (AB - SO),
    b  = (H - HR) / (AB - SO - HR),
    f  = H / (AB - SO),
    BA = H / AB
  ) |>
  left_join(People |> dplyr::select(playerID, nameFirst, nameLast), by = "playerID") |>
  mutate(name = paste(nameFirst, nameLast))

cur <- bat |> filter(yearID == LAST_YEAR) |> arrange(nameLast, nameFirst)
hitter_choices <- setNames(cur$playerID, paste0(cur$name, "  (", sprintf("%.3f", cur$BA), ")"))
default_hitter <- if (any(cur$playerID == "judgeaa01")) "judgeaa01" else cur$playerID[which.max(cur$HR)]

# Logits with a 1/2 continuity correction and delta-method sampling variances
logit_obs <- function(y, n) log((y + 0.5) / (n - y + 0.5))
logit_var <- function(y, n) 1 / (y + 0.5) + 1 / (n - y + 0.5)
add_logits <- function(d) {
  d |> mutate(
    zc = logit_obs(AB - SO, AB),            tc = logit_var(AB - SO, AB),
    zh = logit_obs(HR, AB - SO),            th = logit_var(HR, AB - SO),
    zb = logit_obs(H - HR, AB - SO - HR),   tb = logit_var(H - HR, AB - SO - HR)
  )
}
cur <- add_logits(cur)

inv_logit <- function(x) 1 / (1 + exp(-x))

# Table 3 of the paper: maximum likelihood estimates for 2025 (logit scale)
paper_mu    <- c(c = 1.211, h = -3.059, b = -0.853)
paper_sigma <- c(c = 0.369, h = 0.494,  b = 0.078)
paper_rho   <- c(ch = -0.67, cb = -0.22, hb = -0.23)

make_Sigma <- function(sigma, rho) {
  R <- matrix(c(1, rho["ch"], rho["cb"],
                rho["ch"], 1, rho["hb"],
                rho["cb"], rho["hb"], 1), 3, 3)
  D <- diag(sigma)
  D %*% R %*% D
}
is_pd <- function(S) all(eigen(S, symmetric = TRUE, only.values = TRUE)$values > 1e-10)

# Year-to-year correlations for every consecutive pair of seasons, 2004 on
y2y_pairs <- lapply(2004:(LAST_YEAR - 1), function(y) {
  m <- inner_join(bat |> filter(yearID == y), bat |> filter(yearID == y + 1),
                  by = "playerID", suffix = c(".1", ".2"))
  if (nrow(m) < 30) return(NULL)
  data.frame(year = y, n = nrow(m),
             rate = names(rate_labels),
             r = sapply(names(rate_labels), function(v) cor(m[[paste0(v, ".1")]], m[[paste0(v, ".2")]])))
}) |> bind_rows()
pair_years <- sort(unique(y2y_pairs$year))

# Decade summaries, qualifying seasons 1901-2019 (Table 4, Figure 5)
hist_q <- bat |> filter(yearID >= 1901, yearID <= 2019) |>
  mutate(decade = floor(yearID / 10) * 10) |> add_logits()

decade_obs <- hist_q |> group_by(decade) |>
  summarise(n = n(), sdBA = sd(BA), hbar = mean(h),
            `Contact vs home runs` = cor(c, h),
            `Contact vs BABIP`     = cor(c, b),
            `Home runs vs BABIP`   = cor(h, b),
            `Contact vs finishing` = cor(c, f), .groups = "drop")

# Talent parameters by decade. Noise is removed on the rate scale, as in the paper:
# talent variance = observed variance - average binomial variance p(1 - p) / n.
# Covariances need no correction because the three rates use disjoint sets of at-bats.
# Means and SDs are then carried to the logit scale (delta method) for simulation.
decade_talent <- hist_q |> group_by(decade) |> group_modify(function(d, k) {
  R  <- as.matrix(d[, c("c", "h", "b")])
  n  <- cbind(d$AB, d$AB - d$SO, d$AB - d$SO - d$HR)
  noise <- colMeans(R * (1 - R) / n)
  C  <- cov(R)
  tv <- pmax(diag(C) - noise, 1e-8)
  r  <- function(i, j) max(min(C[i, j] / sqrt(tv[i] * tv[j]), 0.9), -0.9)
  m  <- colMeans(R)
  sl <- sqrt(tv) / (m * (1 - m))
  data.frame(mu_c = qlogis(m[1]), mu_h = qlogis(m[2]), mu_b = qlogis(m[3]),
             s_c = sl[1], s_h = sl[2], s_b = sl[3],
             r_ch = r(1, 2), r_cb = r(1, 3), r_hb = r(2, 3))
}) |> ungroup()
decade_choices <- setNames(decade_talent$decade, paste0(decade_talent$decade, "s"))

# ---------------------------------------------------------------------------
# Look
# ---------------------------------------------------------------------------

col_c   <- "#1f4e79"   # contact
col_h   <- "#c0392b"   # home runs
col_b   <- "#d68910"   # BABIP
col_ink <- "#2b2b2b"
col_mut <- "#8a8f98"

theme_app <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(colour = col_mut),
        strip.text = element_text(face = "bold", size = 13))

fmt3 <- function(x) sub("^0", "", sprintf("%.3f", x))

app_theme <- bs_theme(version = 5, primary = col_c, secondary = col_mut,
                      base_font = font_collection("Source Sans Pro", "Helvetica Neue", "Arial", "sans-serif"))

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- page_navbar(
  title = "Three Rates",
  theme = app_theme,
  bg = col_c,
  inverse = TRUE,
  fillable = FALSE,

  # 1 ------------------------------------------------------------------------
  nav_panel("1  Three steps",
    layout_sidebar(
      sidebar = sidebar(width = 320,
        selectizeInput("hitter", paste0("Hitter (", LAST_YEAR, ", ", MIN_AB, "+ AB)"),
                       choices = hitter_choices, selected = default_hitter),
        hr(),
        h6("What if?"),
        p(class = "small text-muted", "Change one rate and see the batting average it builds. The sliders reset when you pick a new hitter."),
        sliderInput("wc", "Contact rate  c", 0.55, 0.95, 0.70, step = 0.005),
        sliderInput("wh", "Home-run rate  h", 0, 0.16, 0.10, step = 0.001),
        sliderInput("wb", "BABIP  b", 0.22, 0.44, 0.30, step = 0.002),
        uiOutput("whatif")
      ),
      p(class = "lead", "BA = c × [ h + (1 − h) × b ]: of every at-bat, how often the hitter avoids a strikeout, how often those at-bats leave the park, and how often the balls that stay in play fall for hits."),
      layout_columns(col_widths = c(3, 3, 3, 3),
        value_box("Batting average", textOutput("vb_ba"),  theme = value_box_theme(bg = "#f4f6f9", fg = col_ink)),
        value_box("Contact rate",    textOutput("vb_c"),   theme = value_box_theme(bg = "#e8eef5", fg = col_c)),
        value_box("Home-run rate",   textOutput("vb_h"),   theme = value_box_theme(bg = "#f7e7e5", fg = col_h)),
        value_box("BABIP",           textOutput("vb_b"),   theme = value_box_theme(bg = "#fbf0dc", fg = col_b))
      ),
      card(card_header("Where the at-bats went"), plotOutput("flow", height = 230)),
      card(card_header(paste0("How the three rates relate across ", nrow(cur), " hitters in ", LAST_YEAR, " (Figure 1)")),
           plotOutput("pairs", height = 340))
    )
  ),

  # 2 ------------------------------------------------------------------------
  nav_panel("2  Skill or luck?",
    layout_sidebar(
      sidebar = sidebar(width = 300,
        radioButtons("y2y_rate", "Rate", choices = setNames(names(rate_labels), rate_labels), selected = "b"),
        selectInput("y2y_year", "First season of the pair", choices = rev(pair_years), selected = max(pair_years)),
        p(class = "small text-muted", "A skill carries over: hitters high on it one year stay high the next. A rate driven by luck scatters. 2020 is skipped because no one reached 300 at-bats.")
      ),
      layout_columns(col_widths = c(7, 5),
        card(card_header(textOutput("y2y_title")), plotOutput("y2y_scatter", height = 430)),
        card(card_header("Every pair of seasons (Figure 3)"), plotOutput("y2y_all", height = 430))
      )
    )
  ),

  # 3 ------------------------------------------------------------------------
  nav_panel("3  Talent vs luck",
    layout_sidebar(
      sidebar = sidebar(width = 320,
        selectizeInput("sh_hitters", "Highlight hitters", choices = hitter_choices, multiple = TRUE,
                       selected = cur$playerID[cur$nameLast %in% c("Judge", "Andujar", "Aranda", "Wilson") &
                                               cur$nameFirst %in% c("Aaron", "Miguel", "Jonathan", "Jacob")]),
        radioButtons("sh_pair", "Plot", choices = c("Home-run rate vs BABIP" = "hb",
                                                    "Contact rate vs BABIP" = "cb",
                                                    "Contact rate vs home-run rate" = "ch")),
        checkboxInput("sh_corr", "Let the three abilities inform each other (use the correlations)", TRUE),
        p(class = "small text-muted", "Estimates use the 2025 model in Table 3 of the paper: each observed rate is pulled toward the league average, hardest where luck is largest.")
      ),
      layout_columns(col_widths = c(7, 5),
        card(card_header("Observed rate (gray) to estimated ability (blue), Figure 4"), plotOutput("shrink", height = 450)),
        card(card_header("How much of each rate is real?"),
             plotOutput("reliab", height = 230),
             p(class = "small text-muted", "Share of the spread between hitters in one season that reflects real differences in ability, averaged over the 2025 hitters. The rest is luck, and the model removes it by pulling each rate toward the league average."))
      ),
      card(card_header("Observed rates and estimated abilities for the highlighted hitters (Table 2)"),
           tableOutput("sh_table"))
    )
  ),

  # 4 ------------------------------------------------------------------------
  nav_panel("4  Where did .400 go?",
    layout_columns(col_widths = c(6, 6),
      card(card_header("Correlations among the rates by decade, qualifying hitters 1900s–2010s (Figure 5)"),
           plotOutput("decades", height = 380),
           p(class = "small text-muted", "The contact–finishing trade-off (dashed) is almost all contact versus home runs. Contact and BABIP never move much together.")),
      card(card_header("Spread of batting averages by decade"),
           plotOutput("sdba", height = 380),
           p(class = "small text-muted", "Points are seasons of .400 or better; the line is the standard deviation of batting average (right scale)."))
    ),
    p(class = "text-muted mt-2", "To see how much the contact–home-run trade-off matters, try the What-if simulator tab.")
  ),

  # 5 ------------------------------------------------------------------------
  nav_panel("5  What-if simulator",
    card(card_header("How the simulator works"),
      layout_columns(col_widths = c(7, 5),
        div(
          p("Real batting averages mix ability and luck, and the three abilities behind them are linked. This simulator builds an imaginary league of 20,000 hitter-seasons so you can change those links and see what happens to the chance of a .400 season."),
          tags$ol(
            tags$li(tags$b("Build a league. "), "Each simulated hitter gets three true abilities: how often he avoids strikeouts (contact), how often a non-strikeout at-bat is a home run (power), and how often a ball in play falls for a hit (BABIP). The league averages and the spread of ability between hitters are copied from the decade you choose."),
            tags$li(tags$b("Decide how the abilities go together. "), "The three sliders set the correlations between abilities. A negative value is a trade-off: hitters who are good at one thing tend to be worse at the other. Zero means the two abilities are unrelated."),
            tags$li(tags$b("Play a season. "), "Each hitter gets the chosen number of at-bats. Every at-bat is decided by chance using his three abilities, the way a weighted coin is flipped. So even a .350 hitter can have a .320 or .380 season."),
            tags$li(tags$b("Read the results. "), "The histogram shows the season batting averages. The gray outline is the same league with that decade's real correlations, for comparison. The numbers on the right count .400 seasons.")
          )
        ),
        div(class = "p-3 rounded", style = "background:#f4f6f9",
          h6("Things to try"),
          tags$ul(class = "small",
            tags$li("Keep the 2010s and press ", tags$b("Remove the contact–home-run trade-off"), ". The spread of ability widens by about a third and .400 seasons appear, though still rarely."),
            tags$li("Switch to the 1920s or 1930s. These leagues had higher averages and a wider spread, so .400 seasons show up even with the trade-off in place."),
            tags$li("Make the contact–home-run correlation strongly negative (−0.8). Fewer hitters are good at both, so the spread narrows and .400 disappears."),
            tags$li("Move the at-bats slider down. Shorter seasons give luck more room, so extreme averages, high and low, become more common.")
          ),
          p(class = "small text-muted mb-0", "The simulator overstates how often .400 happened in early eras (about 160 per 10,000 seasons for the 1920s, versus roughly 56 in the actual record). Use the counts to compare settings, not as forecasts.")
        )
      )
    ),
    layout_sidebar(
      sidebar = sidebar(width = 340,
        selectInput("sim_dec", "1. Choose the league (decade)", choices = decade_choices, selected = 2010),
        uiOutput("sim_profile"),
        hr(),
        h6("2. Set the trade-offs"),
        actionButton("sim_era", "Reset to this decade's real values", class = "btn-sm btn-outline-primary w-100 mb-1"),
        actionButton("sim_noch", "Remove the contact–home-run trade-off", class = "btn-sm btn-outline-danger w-100 mb-3"),
        sliderInput("r_ch", "Contact & home runs", -0.9, 0.5, -0.6, step = 0.01),
        p(class = "small text-muted mt-n2", "Negative: hitters who rarely strike out tend to hit fewer home runs. This is the trade-off the paper studies."),
        sliderInput("r_cb", "Contact & BABIP", -0.6, 0.6, -0.1, step = 0.01),
        p(class = "small text-muted mt-n2", "Negative: hitters who rarely strike out tend to do worse on balls in play."),
        sliderInput("r_hb", "Home runs & BABIP", -0.6, 0.6, 0, step = 0.01),
        p(class = "small text-muted mt-n2", "Negative: power hitters tend to do worse on balls that stay in the park."),
        hr(),
        h6("3. Set the season length"),
        sliderInput("sim_ab", "At-bats per season", 400, 650, 550, step = 10)
      ),
      layout_columns(col_widths = c(8, 4), fill = FALSE,
        card(fill = FALSE, card_header("Simulated season batting averages"),
             plotOutput("sim_hist", height = 400, fill = FALSE),
             p(class = "small text-muted", "Blue bars: your settings. Gray outline: the same league with the decade's real correlations. Dashed red line: .400.")),
        card(fill = FALSE, card_header("Results"), uiOutput("sim_stats"))
      )
    )
  ),

  nav_spacer(),
  nav_item(tags$span(class = "navbar-text small", paste0("Lahman data through ", LAST_YEAR)))
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # 1 Three steps -----------------------------------------------------------
  hit <- reactive(cur |> filter(playerID == input$hitter))

  observeEvent(input$hitter, {
    d <- hit(); req(nrow(d) == 1)
    updateSliderInput(session, "wc", value = round(d$c, 3))
    updateSliderInput(session, "wh", value = round(d$h, 3))
    updateSliderInput(session, "wb", value = round(d$b, 3))
  })

  output$vb_ba <- renderText(fmt3(hit()$BA))
  output$vb_c  <- renderText(fmt3(hit()$c))
  output$vb_h  <- renderText(fmt3(hit()$h))
  output$vb_b  <- renderText(fmt3(hit()$b))

  output$whatif <- renderUI({
    d  <- hit()
    ba <- input$wc * (input$wh + (1 - input$wh) * input$wb)
    diff <- ba - d$BA
    tagList(
      div(class = "mt-2 p-3 rounded", style = "background:#f4f6f9",
        div(class = "small text-muted", "What-if batting average"),
        div(style = "font-size:2rem;font-weight:700", fmt3(ba)),
        div(class = "small", style = paste0("color:", if (diff >= 0) "#1e8449" else col_h),
            sprintf("%+.0f points vs actual %s", 1000 * diff, fmt3(d$BA)))
      )
    )
  })

  output$flow <- renderPlot({
    d <- hit()
    seg <- data.frame(
      stage = factor(rep(c("At-bats", "Not a strikeout", "Ball in play"), c(2, 2, 2)),
                     levels = rev(c("At-bats", "Not a strikeout", "Ball in play"))),
      part  = c("Strikeouts", "Not a strikeout", "Home runs", "Ball in play", "Hits", "Outs"),
      n     = c(d$SO, d$AB - d$SO, d$HR, d$AB - d$SO - d$HR, d$H - d$HR, d$AB - d$SO - d$HR - (d$H - d$HR))
    ) |> group_by(stage) |> mutate(start = cumsum(n) - n, end = cumsum(n), mid = (start + end) / 2,
                                   share = n / sum(n)) |> ungroup()
    seg$fill <- c("#d5d8dc", col_c, col_h, "#e5c79a", col_b, "#ebe3d4")
    seg$lab  <- paste0(seg$part, "\n", seg$n, " (", sprintf("%.1f%%", 100 * seg$share), ")")
    seg$txt  <- ifelse(seg$fill %in% c(col_c, col_h, col_b), "white", col_ink)
    ggplot(seg) +
      geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(stage) - 0.38, ymax = as.numeric(stage) + 0.38),
                fill = seg$fill, colour = "white") +
      geom_text(aes(x = mid, y = as.numeric(stage), label = ifelse(n / max(end) > 0.07, lab, "")),
                colour = seg$txt, size = 3.8, lineheight = 0.9) +
      scale_y_continuous(breaks = 1:3, labels = levels(seg$stage)) +
      labs(x = NULL, y = NULL,
           subtitle = sprintf("%s: %d HR + %d other hits = %d hits in %d at-bats, a %s average",
                              d$name, d$HR, d$H - d$HR, d$H, d$AB, fmt3(d$BA))) +
      theme_app + theme(panel.grid = element_blank(), axis.text.x = element_blank())
  })

  output$pairs <- renderPlot({
    pr <- list(c("c", "h"), c("c", "b"), c("h", "b"))
    long <- bind_rows(lapply(pr, function(p) data.frame(
      panel = paste0(rate_labels[p[1]], " vs ", rate_labels[p[2]],
                     "\nr = ", sprintf("%.2f", cor(cur[[p[1]]], cur[[p[2]]]))),
      x = cur[[p[1]]], y = cur[[p[2]]], me = cur$playerID == input$hitter)))
    long$panel <- factor(long$panel, levels = unique(long$panel))
    ggplot(long, aes(x, y)) +
      geom_point(colour = col_mut, alpha = 0.55, size = 1.8) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = col_ink, linewidth = 0.7) +
      geom_point(data = long[long$me, ], colour = col_h, size = 4) +
      facet_wrap(~panel, scales = "free") +
      labs(x = NULL, y = NULL) + theme_app
  })

  # 2 Skill or luck ---------------------------------------------------------
  y2y_data <- reactive({
    y <- as.integer(input$y2y_year); v <- input$y2y_rate
    inner_join(bat |> filter(yearID == y), bat |> filter(yearID == y + 1),
               by = "playerID", suffix = c(".1", ".2")) |>
      transmute(x = .data[[paste0(v, ".1")]], y = .data[[paste0(v, ".2")]], name = name.1)
  })

  output$y2y_title <- renderText({
    y <- as.integer(input$y2y_year)
    sprintf("%s, %d vs %d (%d hitters with %d+ AB in both)", rate_labels[input$y2y_rate], y, y + 1, nrow(y2y_data()), MIN_AB)
  })

  output$y2y_scatter <- renderPlot({
    d <- y2y_data(); y <- as.integer(input$y2y_year)
    lim <- range(c(d$x, d$y))
    ggplot(d, aes(x, y)) +
      geom_abline(linetype = "dashed", colour = col_mut) +
      geom_point(colour = col_c, alpha = 0.65, size = 2.3) +
      annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 6, fontface = "bold",
               label = sprintf("r = %.2f", cor(d$x, d$y))) +
      coord_equal(xlim = lim, ylim = lim) +
      labs(x = paste(y, "rate"), y = paste(y + 1, "rate")) + theme_app
  })

  output$y2y_all <- renderPlot({
    short <- c(c = "Contact\nrate", h = "Home-run\nrate", b = "BABIP", BA = "Batting\naverage")
    d <- y2y_pairs |> mutate(sel = year == as.integer(input$y2y_year) & rate == input$y2y_rate,
                             rate = factor(short[rate], levels = short))
    med <- d |> group_by(rate) |> summarise(m = median(r), .groups = "drop")
    ggplot(d, aes(rate, r)) +
      geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), colour = col_mut, size = 2.4, alpha = 0.8) +
      geom_point(data = d[d$sel, ], colour = col_h, size = 4.5) +
      geom_errorbar(data = med, aes(y = m, ymin = m, ymax = m), colour = col_h, width = 0.45, linewidth = 1.2) +
      geom_text(data = med, aes(y = m, label = sprintf("%.2f", m)), nudge_x = 0.38, colour = col_h, fontface = "bold") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = NULL, y = "Year-to-year correlation") + theme_app +
      theme(axis.text.x = element_text(size = 11))
  })

  # 3 Talent vs luck --------------------------------------------------------
  shrunk <- reactive({
    S <- make_Sigma(paper_sigma, if (input$sh_corr) paper_rho else c(ch = 0, cb = 0, hb = 0))
    Z <- cbind(cur$zc, cur$zh, cur$zb) - matrix(paper_mu, nrow(cur), 3, byrow = TRUE)
    est <- t(sapply(seq_len(nrow(cur)), function(i) {
      Ti <- diag(c(cur$tc[i], cur$th[i], cur$tb[i]))
      as.vector(S %*% solve(S + Ti, Z[i, ]))
    })) + matrix(paper_mu, nrow(cur), 3, byrow = TRUE)
    cur |> mutate(ec = inv_logit(est[, 1]), eh = inv_logit(est[, 2]), eb = inv_logit(est[, 3]))
  })

  output$shrink <- renderPlot({
    d <- shrunk(); p <- strsplit(input$sh_pair, "")[[1]]
    d$x0 <- d[[p[1]]]; d$y0 <- d[[p[2]]]; d$x1 <- d[[paste0("e", p[1])]]; d$y1 <- d[[paste0("e", p[2])]]
    d$hl <- d$playerID %in% input$sh_hitters
    ggplot(d) +
      geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), colour = "#c5cbd3",
                   arrow = arrow(length = unit(0.12, "cm")), linewidth = 0.4) +
      geom_point(aes(x0, y0), colour = col_mut, size = 1.4) +
      geom_point(aes(x1, y1), colour = col_c, size = 1.6) +
      geom_segment(data = d[d$hl, ], aes(x = x0, y = y0, xend = x1, yend = y1), colour = col_h,
                   arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
      geom_text(data = d[d$hl, ], aes(x0, y0, label = name), size = 3.8, vjust = -0.8, fontface = "bold", colour = col_ink) +
      labs(x = rate_labels[p[1]], y = rate_labels[p[2]]) + theme_app
  })

  output$reliab <- renderPlot({
    rel <- data.frame(
      rate = factor(rate_labels[c("c", "h", "b")], levels = rev(rate_labels[c("c", "h", "b")])),
      val  = c(mean(paper_sigma["c"]^2 / (paper_sigma["c"]^2 + cur$tc)),
               mean(paper_sigma["h"]^2 / (paper_sigma["h"]^2 + cur$th)),
               mean(paper_sigma["b"]^2 / (paper_sigma["b"]^2 + cur$tb))))
    ggplot(rel, aes(val, rate)) +
      geom_col(aes(x = 1), fill = "#eceff3", width = 0.6) +
      geom_col(fill = c(col_c, col_h, col_b), width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%% ability", 100 * val)), hjust = -0.1, size = 4.5, fontface = "bold") +
      scale_x_continuous(limits = c(0, 1.5), breaks = NULL) +
      labs(x = NULL, y = NULL) + theme_app + theme(panel.grid = element_blank())
  })

  output$sh_table <- renderTable({
    d <- shrunk() |> filter(playerID %in% input$sh_hitters)
    req(nrow(d) > 0)
    data.frame(Hitter = d$name,
               `Contact obs` = fmt3(d$c), `Contact est` = fmt3(d$ec),
               `HR obs` = fmt3(d$h), `HR est` = fmt3(d$eh),
               `BABIP obs` = fmt3(d$b), `BABIP est` = fmt3(d$eb), check.names = FALSE)
  }, striped = TRUE, spacing = "s")

  # 4 Where did .400 go -----------------------------------------------------
  output$decades <- renderPlot({
    long <- decade_obs |> tidyr_long()
    ggplot(long, aes(decade + 5, r, colour = pair, linetype = pair)) +
      geom_hline(yintercept = 0, colour = col_mut) +
      geom_line(linewidth = 1.1) + geom_point(size = 2.2) +
      scale_colour_manual(values = c(`Contact vs home runs` = "#1e8449", `Contact vs BABIP` = col_b,
                                     `Home runs vs BABIP` = "#7d3c98", `Contact vs finishing` = col_ink)) +
      scale_linetype_manual(values = c(`Contact vs home runs` = "solid", `Contact vs BABIP` = "solid",
                                       `Home runs vs BABIP` = "solid", `Contact vs finishing` = "dashed")) +
      labs(x = "Decade", y = "Correlation across qualifying hitters", colour = NULL, linetype = NULL) +
      theme_app + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
  })

  output$sdba <- renderPlot({
    four <- hist_q |> filter(BA >= 0.400)
    k <- 0.40 / 0.04
    ggplot() +
      geom_point(data = hist_q, aes(yearID, BA), colour = "#d5d8dc", size = 0.6, alpha = 0.5) +
      geom_point(data = four, aes(yearID, BA), colour = col_h, size = 2.5) +
      geom_hline(yintercept = 0.4, linetype = "dotted", colour = col_h) +
      geom_line(data = decade_obs, aes(decade + 5, sdBA * k), colour = col_c, linewidth = 1.2) +
      geom_point(data = decade_obs, aes(decade + 5, sdBA * k), colour = col_c, size = 2) +
      scale_y_continuous("Batting average", sec.axis = sec_axis(~ . / k, name = "SD of batting average")) +
      labs(x = NULL) + theme_app +
      theme(axis.title.y.right = element_text(colour = col_c), axis.text.y.right = element_text(colour = col_c))
  })

  era <- reactive(decade_talent |> filter(decade == as.integer(input$sim_dec)))

  set_era_rho <- function() {
    e <- era()
    updateSliderInput(session, "r_ch", value = round(e$r_ch, 2))
    updateSliderInput(session, "r_cb", value = round(e$r_cb, 2))
    updateSliderInput(session, "r_hb", value = round(e$r_hb, 2))
  }
  observeEvent(input$sim_dec, set_era_rho())
  observeEvent(input$sim_era, set_era_rho())
  observeEvent(input$sim_noch, updateSliderInput(session, "r_ch", value = 0))

  simulate_league <- function(e, rho, AB, N = 20000) {
    S <- make_Sigma(c(c = e$s_c, h = e$s_h, b = e$s_b), rho)
    if (!is_pd(S)) return(NULL)
    set.seed(2026)
    T  <- mvrnorm(N, c(e$mu_c, e$mu_h, e$mu_b), S)
    cc <- inv_logit(T[, 1]); hh <- inv_logit(T[, 2]); bb <- inv_logit(T[, 3])
    nk  <- rbinom(N, AB, cc)
    hr  <- rbinom(N, nk, hh)
    bip <- rbinom(N, nk - hr, bb)
    list(true = cc * (hh + (1 - hh) * bb), obs = (hr + bip) / AB, N = N)
  }

  sim <- reactive({
    s <- simulate_league(era(), c(ch = input$r_ch, cb = input$r_cb, hb = input$r_hb), input$sim_ab)
    validate(need(!is.null(s), "These three correlations can't all hold at once. Move one of the sliders toward zero."))
    s
  })
  base <- reactive({
    e <- era()
    simulate_league(e, c(ch = e$r_ch, cb = e$r_cb, hb = e$r_hb), input$sim_ab)
  })

  output$sim_profile <- renderUI({
    e <- era()
    m  <- plogis(c(e$mu_c, e$mu_h, e$mu_b))
    sd <- c(e$s_c, e$s_h, e$s_b) * m * (1 - m)
    ba <- m[1] * (m[2] + (1 - m[2]) * m[3])
    row <- function(lbl, mm, ss) tags$tr(tags$td(lbl), tags$td(fmt3(mm)), tags$td(fmt3(ss)))
    tagList(
      p(class = "small text-muted mb-1", "What a typical hitter looked like, and how much true ability varied:"),
      tags$table(class = "table table-sm small mb-1",
        tags$thead(tags$tr(tags$th(""), tags$th("Average"), tags$th("Spread (SD)"))),
        tags$tbody(row("Contact rate", m[1], sd[1]), row("Home-run rate", m[2], sd[2]), row("BABIP", m[3], sd[3]))),
      p(class = "small text-muted", "Typical batting average: ", tags$b(fmt3(ba)))
    )
  })

  output$sim_hist <- renderPlot({
    s <- sim(); b <- base()
    p <- ggplot(data.frame(BA = s$obs), aes(BA)) +
      geom_histogram(binwidth = 0.005, boundary = 0, fill = col_c, colour = "white", alpha = 0.85)
    if (!is.null(b)) p <- p + geom_step(data = data.frame(BA = b$obs) |>
                                          mutate(bin = floor(BA / 0.005) * 0.005) |> count(bin),
                                        aes(x = bin, y = n), colour = col_mut, linewidth = 0.8, direction = "hv")
    p + geom_vline(xintercept = 0.400, colour = col_h, linewidth = 1, linetype = "dashed") +
      annotate("text", x = 0.400, y = Inf, label = ".400", colour = col_h, vjust = 1.5, hjust = -0.2, fontface = "bold") +
      scale_x_continuous(labels = fmt3) +
      coord_cartesian(xlim = c(0.15, 0.45)) +
      labs(x = "Season batting average", y = "Number of hitter-seasons",
           subtitle = sprintf("%s simulated hitter-seasons, %ss league, %d at-bats each",
                              format(s$N, big.mark = ","), input$sim_dec, input$sim_ab)) +
      theme_app
  })

  output$sim_stats <- renderUI({
    s <- sim(); b <- base()
    per10k <- function(x) 1e4 * mean(x >= 0.400)
    stat <- function(lbl, val, ref, help, col = col_ink) div(class = "mb-3",
      div(class = "small fw-semibold", lbl),
      div(style = paste0("font-size:1.7rem;font-weight:700;color:", col), val),
      div(class = "small text-muted", paste0("Decade's real values: ", ref)),
      div(class = "small text-muted fst-italic", help))
    tagList(
      stat(".400 seasons per 10,000", format(round(per10k(s$obs), 1), nsmall = 1),
           format(round(per10k(b$obs), 1), nsmall = 1),
           "How often a hitter would bat .400 or better.", col_h),
      stat("Spread of true ability (SD)", fmt3(sd(s$true)), fmt3(sd(b$true)),
           "How far apart hitters really are, before luck. A wider spread leaves more room for a .400 hitter.", col_c),
      stat("Spread of season averages (SD)", fmt3(sd(s$obs)), fmt3(sd(b$obs)),
           "What you would see in the standings: true spread plus a season's luck.")
    )
  })
}

# Reshape decade correlations for plotting without needing tidyr
tidyr_long <- function(d) {
  pairs <- c("Contact vs home runs", "Contact vs BABIP", "Home runs vs BABIP", "Contact vs finishing")
  bind_rows(lapply(pairs, function(p) data.frame(decade = d$decade, pair = p, r = d[[p]]))) |>
    mutate(pair = factor(pair, levels = pairs))
}

shinyApp(ui, server)
