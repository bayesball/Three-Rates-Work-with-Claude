# analysis.R
# Reproduces every number, table and figure in "Batting Average as the Product
# of Two Rates". Run from this folder:  Rscript analysis.R
# Writes figures to figures/*.pdf and numeric results to numbers.tex, which
# the paper reads with \input{numbers.tex}.

suppressMessages({
  library(Lahman)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(scales)
  library(ggrepel)
  library(purrr)
})

dir.create("figures", showWarnings = FALSE)
set.seed(2025)

# ---- output helpers ---------------------------------------------------------
macros <- character(0)
mac <- function(name, value) {
  macros <<- c(macros, sprintf("\\newcommand{\\%s}{%s}", name, value))
}
f2 <- function(x) sprintf("%.2f", x)
f3 <- function(x) sprintf("%.3f", x)
f4 <- function(x) sprintf("%.4f", x)
p1 <- function(x) sprintf("%.1f\\%%", 100 * x)
p0 <- function(x) sprintf("%.0f\\%%", 100 * x)
tex_rows <- function(df) paste(apply(df, 1, paste, collapse = " & "), collapse = " \\\\\n")

blue <- "#2a78d6"; orange <- "#eb6834"; red <- "#e34948"
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; hair <- "#e1e0d9"
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = hair, linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink),
        strip.text = element_text(face = "bold"), legend.position = "none")
theme_set(theme_paper)
save_fig <- function(p, file, w = 6, h = 4) ggsave(file.path("figures", file), p, width = w, height = h, device = cairo_pdf)

# ---- 2025 sample --------------------------------------------------------------
min_ab <- 300
people_nm <- People |> transmute(playerID, name = paste(nameFirst, nameLast))

b25 <- Batting |> filter(yearID == 2025) |>
  group_by(playerID) |>
  summarise(AB = sum(AB, na.rm = TRUE), SO = sum(SO, na.rm = TRUE), H = sum(H, na.rm = TRUE),
            .groups = "drop") |>
  filter(AB >= min_ab) |>
  mutate(BA = H / AB, contact = (AB - SO) / AB, finish = H / (AB - SO)) |>
  inner_join(people_nm, by = "playerID")

mac("nQual", nrow(b25))
mac("minAB", min_ab)
max_err <- max(abs(b25$BA - b25$contact * b25$finish))
mac("identityErr", sprintf("%.0e", max_err))

# ---- case-study table -------------------------------------------------------------
cases <- c("Aaron Judge", "Jacob Wilson", "Miguel Andujar", "Jonathan Aranda")
case_tbl <- b25 |> filter(name %in% cases) |> arrange(desc(BA))
stopifnot(nrow(case_tbl) == 4)
writeLines(tex_rows(case_tbl |> transmute(name, AB, H, SO, f3(BA), f3(contact), f3(finish))),
           "tab_cases.tex")
judge <- filter(b25, name == "Aaron Judge"); wilson <- filter(b25, name == "Jacob Wilson")
mac("judgeBA", f3(judge$BA)); mac("wilsonBA", f3(wilson$BA))
mac("judgeC", f3(judge$contact)); mac("judgeF", f3(judge$finish))
mac("wilsonC", f3(wilson$contact)); mac("wilsonF", f3(wilson$finish))
ord <- function(n) paste0(n, ifelse(n %% 100 %in% 11:13, "th", c("th","st","nd","rd",rep("th",6))[n %% 10 + 1]))
mac("judgeCrank", ord(sum(b25$contact < judge$contact) + 1))
mac("wilsonCrank", sum(b25$contact > wilson$contact) + 1)

# ---- Figure 1: the 2025 field -----------------------------------------------------
r_cf <- cor(b25$contact, b25$finish)
mac("rcf", f2(r_cf))
r_ci <- tanh(atanh(r_cf) + c(-1, 1) * qnorm(0.975) / sqrt(nrow(b25) - 3))
mac("rcfLo", f2(r_ci[1])); mac("rcfHi", f2(r_ci[2]))
contours <- crossing(ba = c(.200, .250, .300),
                     contact = seq(min(b25$contact), max(b25$contact), length.out = 300)) |>
  mutate(finish = ba / contact) |> filter(finish >= min(b25$finish), finish <= max(b25$finish))
clab <- contours |> group_by(ba) |> slice_max(contact, n = 1)
p1f <- ggplot(b25, aes(contact, finish)) +
  geom_line(data = contours, aes(group = ba), colour = muted, linetype = "22", linewidth = 0.4) +
  geom_text(data = clab, aes(label = sub("^0", "", f3(ba))), hjust = -0.15, size = 3, colour = ink2) +
  geom_point(colour = blue, size = 1.6, alpha = 0.6) +
  geom_point(data = case_tbl, shape = 21, size = 3, stroke = 0.8, colour = ink) +
  geom_text_repel(data = case_tbl, aes(label = name), size = 3, seed = 1, min.segment.length = 0) +
  scale_x_continuous(labels = label_number(accuracy = .01), expand = expansion(mult = c(.02, .08))) +
  scale_y_continuous(labels = label_number(accuracy = .01)) +
  labs(x = "Strikeout-avoidance rate, c", y = "Finishing rate, f")
save_fig(p1f, "fig1_scatter2025.pdf", 6, 4.2)

# ---- variance decomposition -----------------------------------------------------------
lc <- log(b25$contact); lf <- log(b25$finish); lba <- log(b25$BA)
v_c <- var(lc); v_f <- var(lf); cv <- cov(lc, lf); v_ba <- var(lba)
mac("vc", f4(v_c)); mac("vf", f4(v_f)); mac("twocov", f4(2 * cv)); mac("vba", f4(v_ba))
mac("vind", f4(v_c + v_f)); mac("pctred", p0(1 - v_ba / (v_c + v_f)))

# ---- year-to-year stability, 2004-2025 -------------------------------------------------
panel <- Batting |> filter(yearID >= 2004) |>
  group_by(playerID, yearID) |>
  summarise(AB = sum(AB, na.rm = TRUE), SO = sum(SO, na.rm = TRUE), H = sum(H, na.rm = TRUE),
            .groups = "drop") |>
  filter(AB >= min_ab) |>
  mutate(contact = (AB - SO) / AB, finish = H / (AB - SO), BA = H / AB)
yrs <- sort(unique(panel$yearID))
pairs <- tibble(y1 = yrs[-length(yrs)], y2 = yrs[-1]) |> filter(y2 - y1 == 1)
stab <- pairs |> mutate(res = map2(y1, y2, function(a, b) {
  m <- inner_join(filter(panel, yearID == a), filter(panel, yearID == b),
                  by = "playerID", suffix = c("_1", "_2"))
  tibble(n = nrow(m), c = cor(m$contact_1, m$contact_2), f = cor(m$finish_1, m$finish_2),
         ba = cor(m$BA_1, m$BA_2))
})) |> unnest(res)
mac("nPairs", nrow(stab))
mac("medPairN", round(median(stab$n)))
for (v in c("c", "f", "ba")) {
  nm <- c(c = "C", f = "F", ba = "BA")[v]
  mac(paste0("med", nm), f2(median(stab[[v]])))
  mac(paste0("min", nm), f2(min(stab[[v]])))
  mac(paste0("max", nm), f2(max(stab[[v]])))
}

m2425 <- inner_join(filter(panel, yearID == 2024), filter(panel, yearID == 2025),
                    by = "playerID", suffix = c("_1", "_2"))
mac("nBoth", nrow(m2425))
mac("rCtwo", f2(cor(m2425$contact_1, m2425$contact_2)))
mac("rFtwo", f2(cor(m2425$finish_1, m2425$finish_2)))
long2425 <- bind_rows(
  m2425 |> transmute(x = contact_1, y = contact_2, rate = "Strikeout avoidance, c"),
  m2425 |> transmute(x = finish_1, y = finish_2, rate = "Finishing rate, f")) |>
  mutate(rate = factor(rate, levels = c("Strikeout avoidance, c", "Finishing rate, f")))
rl <- long2425 |> group_by(rate) |> summarise(r = cor(x, y))
p2f <- ggplot(long2425, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = muted, linewidth = 0.4) +
  geom_point(aes(colour = rate), size = 1.4, alpha = 0.7) +
  geom_text(data = rl, aes(x = -Inf, y = Inf, label = paste0("r = ", f2(r))),
            hjust = -0.2, vjust = 1.6, size = 3.2, colour = ink) +
  scale_colour_manual(values = c(blue, orange)) +
  scale_x_continuous(labels = label_number(accuracy = .01)) +
  scale_y_continuous(labels = label_number(accuracy = .01)) +
  facet_wrap(~ rate, scales = "free") + labs(x = "2024 rate", y = "2025 rate")
save_fig(p2f, "fig2_repeat.pdf", 6.4, 3.2)

sl <- stab |> select(y1, c, f, ba) |> pivot_longer(-y1, names_to = "rate", values_to = "r") |>
  mutate(rate = factor(rate, levels = c("c", "f", "ba"),
                       labels = c("Strikeout\navoidance, c", "Finishing\nrate, f", "Batting\naverage")))
meds <- sl |> group_by(rate) |> summarise(m = median(r))
p3f <- ggplot(sl, aes(rate, r)) +
  geom_jitter(width = 0.08, height = 0, size = 1.6, alpha = 0.55, colour = ink2) +
  geom_crossbar(data = meds, aes(y = m, ymin = m, ymax = m), width = 0.45, colour = red, linewidth = 0.6) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = NULL, y = "Year-to-year correlation")
save_fig(p3f, "fig3_stability.pdf", 5, 3.4)

# ---- joint bivariate logistic-normal model ------------------------------------------------
long <- bind_rows(
  b25 |> transmute(playerID, type = "c", succ = AB - SO, fail = SO),
  b25 |> transmute(playerID, type = "f", succ = H, fail = (AB - SO) - H)) |>
  mutate(type = factor(type, levels = c("c", "f")))
fit <- glmer(cbind(succ, fail) ~ 0 + type + (0 + type | playerID), data = long,
             family = binomial, control = glmerControl(optimizer = "bobyqa"))
fe <- fixef(fit); se <- sqrt(diag(vcov(fit)))
vc <- as.data.frame(VarCorr(fit))
s_c <- vc$sdcor[vc$var1 == "typec" & is.na(vc$var2)]
s_f <- vc$sdcor[vc$var1 == "typef" & is.na(vc$var2)]
rho <- vc$sdcor[!is.na(vc$var2)]
Sigma <- matrix(c(s_c^2, rho * s_c * s_f, rho * s_c * s_f, s_f^2), 2)

# Extract conditional modes BEFORE profiling: profile()/confint() modify the
# fitted merMod object in place, after which ranef() no longer returns the
# modes at the maximum likelihood estimates.
re <- ranef(fit)$playerID
# profile-likelihood intervals for the variance components
ci <- suppressWarnings(suppressMessages(confint(fit, parm = "theta_", method = "profile", oldNames = FALSE)))
ci_s_c <- ci[grepl("sd_typec", rownames(ci)), ]; ci_s_f <- ci[grepl("sd_typef", rownames(ci)), ]
ci_rho <- ci[grepl("cor_", rownames(ci)), ]

mac("muc", f3(fe[["typec"]])); mac("muf", f3(fe[["typef"]]))
mac("muclo", f3(fe[["typec"]] - 1.96 * se[1])); mac("muchi", f3(fe[["typec"]] + 1.96 * se[1]))
mac("muflo", f3(fe[["typef"]] - 1.96 * se[2])); mac("mufhi", f3(fe[["typef"]] + 1.96 * se[2]))
mac("pmuc", f3(plogis(fe[["typec"]]))); mac("pmuf", f3(plogis(fe[["typef"]])))
mac("sigc", f3(s_c)); mac("sigf", f3(s_f)); mac("rhomc", f2(rho))
mac("sigclo", f3(ci_s_c[1])); mac("sigchi", f3(ci_s_c[2]))
mac("sigflo", f3(ci_s_f[1])); mac("sigfhi", f3(ci_s_f[2]))
mac("rholo", f2(ci_rho[1])); mac("rhohi", f2(ci_rho[2]))

b25 <- b25 |> mutate(tau_c2 = 1 / (AB * contact * (1 - contact)),
                     tau_f2 = 1 / ((AB - SO) * finish * (1 - finish)),
                     w_c = s_c^2 / (s_c^2 + tau_c2), w_f = s_f^2 / (s_f^2 + tau_f2))
mac("relc", f2(mean(b25$w_c))); mac("relf", f2(mean(b25$w_f)))
mac("tauc", f3(mean(sqrt(b25$tau_c2)))); mac("tauf", f3(mean(sqrt(b25$tau_f2))))
mac("relcPct", p0(mean(b25$w_c))); mac("relfPct", p0(mean(b25$w_f)))

b25 <- b25 |> mutate(c_s = plogis(fe[["typec"]] + re[playerID, "typec"]),
                     f_s = plogis(fe[["typef"]] + re[playerID, "typef"]))
# Judge worked example: closed-form bivariate shrinkage vs. lme4 conditional modes
judge <- filter(b25, name == "Aaron Judge")
z <- c(qlogis(judge$contact) - fe[["typec"]], qlogis(judge$finish) - fe[["typef"]])
bj <- Sigma %*% solve(Sigma + diag(c(judge$tau_c2, judge$tau_f2))) %*% z
mac("judgeZc", f3(z[1])); mac("judgeZf", f3(z[2]))
mac("judgeCs", f3(plogis(fe[["typec"]] + bj[1]))); mac("judgeFs", f3(plogis(fe[["typef"]] + bj[2])))
mac("judgeCl", f3(judge$c_s)); mac("judgeFl", f3(judge$f_s))
mac("mnAbsDc", f3(mean(abs(b25$c_s - b25$contact))))
mac("mnAbsDf", f3(mean(abs(b25$f_s - b25$finish))))

cs <- b25 |> filter(name %in% cases) |> arrange(desc(BA))
writeLines(tex_rows(cs |> transmute(name, f3(contact), f3(c_s), f3(finish), f3(f_s))), "tab_shrink.tex")

p4f <- ggplot(b25, aes(contact, finish)) +
  geom_segment(aes(xend = c_s, yend = f_s), colour = muted, alpha = 0.5, linewidth = 0.3,
               arrow = arrow(length = unit(0.04, "cm"))) +
  geom_point(colour = muted, size = 0.9, alpha = 0.6) +
  geom_point(aes(x = c_s, y = f_s), colour = blue, size = 1.1, alpha = 0.8) +
  geom_segment(data = cs, aes(xend = c_s, yend = f_s), colour = red, linewidth = 0.7,
               arrow = arrow(length = unit(0.14, "cm"))) +
  geom_point(data = cs, colour = ink, size = 1.8) +
  geom_text_repel(data = cs, aes(label = name), size = 3, seed = 1, min.segment.length = 0) +
  scale_x_continuous(labels = label_number(accuracy = .01)) +
  scale_y_continuous(labels = label_number(accuracy = .01)) +
  labs(x = "Strikeout-avoidance rate, c", y = "Finishing rate, f")
save_fig(p4f, "fig4_shrinkage.pdf", 6, 4.2)

# ---- historical analysis, 1901-2025 ----------------------------------------------------
# Historical series: AL, NL and Federal League (the Negro Leagues, added to the
# Lahman database in 2024, have incomplete strikeout records).
mlb_lgs <- c("AL", "NL", "FL")
hist <- Batting |> filter(yearID >= 1901, lgID %in% mlb_lgs) |>
  group_by(playerID, yearID) |>
  summarise(AB = sum(AB, na.rm = TRUE), SO = sum(SO, na.rm = TRUE), H = sum(H, na.rm = TRUE),
            nNA = sum(is.na(SO)), .groups = "drop") |>
  filter(AB >= min_ab, nNA == 0) |>   # drop seasons with any missing SO
  mutate(BA = H / AB, contact = (AB - SO) / AB, finish = H / (AB - SO), decade = (yearID %/% 10) * 10)
n_drop <- Batting |> filter(yearID >= 1901, lgID %in% mlb_lgs) |> group_by(playerID, yearID) |>
  summarise(AB = sum(AB, na.rm = TRUE), nNA = sum(is.na(SO)), .groups = "drop") |>
  filter(AB >= min_ab, nNA > 0) |> nrow()
mac("nDropSO", n_drop)
n_negro <- Batting |> filter(yearID >= 1901, !lgID %in% mlb_lgs) |> group_by(playerID, yearID) |>
  summarise(AB = sum(AB, na.rm = TRUE), .groups = "drop") |> filter(AB >= min_ab) |> nrow()
mac("nNegro", n_negro)
q400 <- hist |> filter(BA >= .400) |> inner_join(people_nm, by = "playerID") |> arrange(yearID)
mac("nFourHundred", nrow(q400)); mac("lastFourHundred", max(q400$yearID))
mac("firstHistYear", min(hist$yearID))
dec <- hist |> group_by(decade) |>
  summarise(n = n(), sd_BA = sd(BA), var_lc = var(log(contact)), var_lf = var(log(finish)),
            cov_cf = cov(log(contact), log(finish)), var_lba = var(log(BA)),
            r_cf = cor(contact, finish),
            # expected binomial sampling variance of log(rate) (delta method)
            nv_lc = mean((1 - contact) / (AB * contact)),
            nv_lf = mean((1 - finish) / ((AB - SO) * finish)), .groups = "drop") |>
  mutate(tv_lc = var_lc - nv_lc, tv_lf = var_lf - nv_lf,
         r_talent = cov_cf / sqrt(tv_lc * tv_lf))
d1900 <- filter(dec, decade == 1900); d1940 <- filter(dec, decade == 1940); d2010 <- filter(dec, decade == 2010)
mac("sdBAnineteen", f3(d1900$sd_BA)); mac("sdBAforties", f3(d1940$sd_BA)); mac("sdBAtens", f3(d2010$sd_BA))
mac("sdDrop", p0(1 - d1940$sd_BA / d1900$sd_BA))
mac("vcRatio", sprintf("%.1f", d2010$var_lc / d1900$var_lc))
mac("vfRatio", sprintf("%.2f", d2010$var_lf / d1900$var_lf))
mac("rcfNineteen", f2(d1900$r_cf)); mac("rcfTens", f2(d2010$r_cf))
post60 <- filter(dec, decade >= 1960, decade <= 2010)
mac("rcfPostLo", f2(min(post60$r_cf))); mac("rcfPostHi", f2(max(post60$r_cf)))
cf_var <- d2010$var_lc + d2010$var_lf + 2 * d1900$cov_cf
mac("vlbaTens", f4(d2010$var_lba)); mac("vlbaNineteen", f4(d1900$var_lba)); mac("vlbaCF", f4(cf_var))
mac("vlcTens", f4(d2010$var_lc)); mac("vlfTens", f4(d2010$var_lf))
mac("tvlcNineteen", f4(d1900$tv_lc)); mac("tvlcTens", f4(d2010$tv_lc))
mac("tvlfNineteen", f4(d1900$tv_lf)); mac("tvlfTens", f4(d2010$tv_lf))
mac("tvcRatio", sprintf("%.1f", d2010$tv_lc / d1900$tv_lc)); mac("tvfRatio", sprintf("%.2f", d2010$tv_lf / d1900$tv_lf))
mac("nvlfShareNineteen", p0(d1900$nv_lf / d1900$var_lf)); mac("nvlfShareTens", p0(d2010$nv_lf / d2010$var_lf))
mac("rTalNineteen", f2(d1900$r_talent)); mac("rTalTens", f2(d2010$r_talent))
cf_tal <- d2010$tv_lc + d2010$tv_lf + 2 * d1900$cov_cf
act_tal <- d2010$tv_lc + d2010$tv_lf + 2 * d2010$cov_cf
tal1900 <- d1900$tv_lc + d1900$tv_lf + 2 * d1900$cov_cf
mac("talCF", f4(cf_tal)); mac("talTens", f4(act_tal)); mac("talNineteen", f4(tal1900))
mac("covNineteen", f4(d1900$cov_cf)); mac("covTens", f4(d2010$cov_cf))
writeLines(tex_rows(dec |> filter(decade <= 2010) |>
  transmute(paste0(decade, "s"), n, f3(sd_BA), f4(var_lc), f4(var_lf), f4(2 * cov_cf), f4(var_lba), f2(r_cf),
            f4(tv_lf), f2(r_talent))),
  "tab_decades.tex")
mac("nDecTwenties", filter(dec, decade == 2020)$n)

dec_plot <- dec |> filter(decade <= 2010)
panel_a <- dec_plot |> transmute(decade, value = sd_BA, panel = "(a) SD of batting average", series = "BA")
panel_b <- dec_plot |> transmute(decade, `Var(log c)` = var_lc, `Var(log f)` = var_lf) |>
  pivot_longer(-decade, names_to = "series") |> mutate(panel = "(b) Var(log c) and Var(log f)")
panel_c <- dec_plot |> transmute(decade, value = r_cf, panel = "(c) Correlation of c and f", series = "r")
pd <- bind_rows(panel_a, panel_b, panel_c) |>
  mutate(panel = factor(panel, levels = unique(c(panel_a$panel, panel_b$panel, panel_c$panel))))
lab_b <- panel_b |> group_by(series) |> slice_max(decade, n = 1) |>
  mutate(panel = factor("(b) Var(log c) and Var(log f)", levels = levels(pd$panel)))
rug <- q400 |> transmute(x = yearID, panel = factor("(a) SD of batting average", levels = levels(pd$panel)))
p5f <- ggplot(pd, aes(decade + 5, value, colour = series)) +
  geom_hline(data = tibble(panel = factor("(c) Correlation of c and f", levels = levels(pd$panel)), y = 0),
             aes(yintercept = y), colour = muted, linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "t", colour = red, linewidth = 0.4) +
  geom_text(data = lab_b, aes(label = sub("Var\\(log (.)\\)", "\\1", series)), hjust = -0.6, size = 3.2, fontface = "italic") +
  scale_colour_manual(values = c(BA = blue, `Var(log c)` = blue, `Var(log f)` = orange, r = red)) +
  scale_x_continuous(breaks = seq(1900, 2020, 40), expand = expansion(mult = c(.04, .12))) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) + labs(x = "Decade", y = NULL)
save_fig(p5f, "fig5_history.pdf", 7.5, 3)

# ---- write macros -----------------------------------------------------------------
mac("lahmanVersion", as.character(packageVersion("Lahman")))
mac("lmeVersion", as.character(packageVersion("lme4")))
mac("RVersion", paste(R.version$major, R.version$minor, sep = "."))
writeLines(c("% Generated by analysis.R -- do not edit by hand", macros), "numbers.tex")
cat("Wrote", length(macros), "macros\n")
