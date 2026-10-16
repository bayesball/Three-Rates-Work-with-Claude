# analysis.R
# Reproduces every number, table and figure in "Contact, Power, and Balls in
# Play: A Three-Rate Decomposition of Batting Average". Run from this folder:
#   Rscript analysis.R
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
p0 <- function(x) sprintf("%.0f\\%%", 100 * x)
p1 <- function(x) sprintf("%.1f\\%%", 100 * x)
tex_rows <- function(df) paste(apply(df, 1, paste, collapse = " & "), collapse = " \\\\\n")

blue <- "#2a78d6"; orange <- "#eb6834"; aqua <- "#1baf7a"; violet <- "#4a3aa7"; red <- "#e34948"
ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; hair <- "#e1e0d9"
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = hair, linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink),
        strip.text = element_text(face = "bold"), legend.position = "none")
theme_set(theme_paper)
save_fig <- function(p, file, w = 6, h = 4) ggsave(file.path("figures", file), p, width = w, height = h, device = cairo_pdf)

rate_cols <- c(c = blue, h = aqua, b = orange)

# ---- 2025 sample --------------------------------------------------------------
min_ab <- 300
people_nm <- People |> transmute(playerID, name = paste(nameFirst, nameLast))

b25 <- Batting |> filter(yearID == 2025) |>
  group_by(playerID) |>
  summarise(across(c(AB, SO, H, HR), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
  filter(AB >= min_ab) |>
  mutate(BA = H / AB, c = (AB - SO) / AB, h = HR / (AB - SO), b = (H - HR) / (AB - SO - HR),
         f = H / (AB - SO)) |>
  inner_join(people_nm, by = "playerID")

mac("nQual", nrow(b25)); mac("minAB", min_ab)
mac("identityErr", sprintf("%.0e", max(abs(b25$BA - b25$c * (b25$h + (1 - b25$h) * b25$b)))))
mac("nZeroHR", sum(b25$HR == 0))

# ---- case studies ----------------------------------------------------------------
judge <- filter(b25, name == "Aaron Judge"); wilson <- filter(b25, name == "Jacob Wilson")
band <- b25 |> filter(abs(BA - judge$BA) <= .02, !name %in% c("Aaron Judge", "Jacob Wilson"))
power_pick <- band |> slice_max(h, n = 1, with_ties = FALSE)
babip_pick <- band |> filter(name != power_pick$name) |> slice_max(b, n = 1, with_ties = FALSE)
cases <- bind_rows(judge, wilson, power_pick, babip_pick) |> arrange(desc(BA))
writeLines(tex_rows(cases |> transmute(name, AB, H, SO, HR, f3(BA), f3(c), f3(h), f3(b))), "tab_cases.tex")
mac("powerName", power_pick$name); mac("babipName", babip_pick$name)
mac("powerH", f3(power_pick$h)); mac("babipB", f3(babip_pick$b))
mac("judgeC", f3(judge$c)); mac("judgeBA", f3(judge$BA)); mac("judgeH", f3(judge$h)); mac("judgeB", f3(judge$b))
mac("wilsonBA", f3(wilson$BA)); mac("wilsonH", f3(wilson$h)); mac("wilsonB", f3(wilson$b))

# ---- pairwise correlations, 2025 ----------------------------------------------------
rci <- function(x, y) { r <- cor(x, y); n <- length(x)
  c(r, tanh(atanh(r) + c(-1, 1) * qnorm(.975) / sqrt(n - 3))) }
r_ch <- rci(b25$c, b25$h); r_cb <- rci(b25$c, b25$b); r_hb <- rci(b25$h, b25$b); r_cf <- cor(b25$c, b25$f)
mac("rch", f2(r_ch[1])); mac("rchLo", f2(r_ch[2])); mac("rchHi", f2(r_ch[3]))
mac("rcb", f2(r_cb[1])); mac("rcbLo", f2(r_cb[2])); mac("rcbHi", f2(r_cb[3]))
mac("rhb", f2(r_hb[1])); mac("rhbLo", f2(r_hb[2])); mac("rhbHi", f2(r_hb[3]))
mac("rcf", f2(r_cf))

pairs_df <- bind_rows(
  b25 |> transmute(x = c, y = h, pair = "Contact rate vs home-run rate"),
  b25 |> transmute(x = c, y = b, pair = "Contact rate vs BABIP"),
  b25 |> transmute(x = h, y = b, pair = "Home-run rate vs BABIP")) |>
  mutate(pair = factor(pair, levels = c("Contact rate vs home-run rate", "Contact rate vs BABIP", "Home-run rate vs BABIP")))
rlab <- tibble(pair = factor(levels(pairs_df$pair), levels = levels(pairs_df$pair)),
               r = c(r_ch[1], r_cb[1], r_hb[1]))
p1f <- ggplot(pairs_df, aes(x, y)) +
  geom_point(colour = blue, size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = ink2, linewidth = 0.5) +
  geom_text(data = rlab, aes(x = -Inf, y = Inf, label = paste0("r = ", f2(r))),
            hjust = -0.15, vjust = 1.5, size = 3.2, colour = ink) +
  facet_wrap(~ pair, scales = "free") +
  scale_x_continuous(labels = label_number(accuracy = .01)) +
  scale_y_continuous(labels = label_number(accuracy = .01)) +
  labs(x = NULL, y = NULL)
save_fig(p1f, "fig1_pairs2025.pdf", 7, 2.7)

# ---- variance decomposition, 2025 ----------------------------------------------------
lc <- log(b25$c); lf <- log(b25$f); lba <- log(b25$BA)
mac("vc", f4(var(lc))); mac("vf", f4(var(lf))); mac("twocov", f4(2 * cov(lc, lf))); mac("vba", f4(var(lba)))
# delta-method split of Var(log f) into HR.Rate, BABIP and covariance parts
dm <- function(h, b, f) {
  hb <- mean(h); bb <- mean(b); fb <- mean(f)
  c(hr = (1 - bb)^2 * var(h) / fb^2, babip = (1 - hb)^2 * var(b) / fb^2,
    cov = 2 * (1 - bb) * (1 - hb) * cov(h, b) / fb^2)
}
a25 <- dm(b25$h, b25$b, b25$f)
mac("aHR", f4(a25["hr"])); mac("aB", f4(a25["babip"])); mac("aCov", f4(a25["cov"]))
mac("aTot", f4(sum(a25))); mac("aGap", p1(abs(sum(a25) - var(lf)) / var(lf)))
mac("aHRshare", p0(a25["hr"] / sum(a25))); mac("aBshare", p0(a25["babip"] / sum(a25)))
mac("aCovshare", p0(a25["cov"] / sum(a25)))
# the same split with binomial sampling noise removed from Var(HR.Rate) and Var(BABIP)
tvh25 <- var(b25$h) - mean(b25$h * (1 - b25$h) / (b25$AB - b25$SO))
tvb25 <- var(b25$b) - mean(b25$b * (1 - b25$b) / (b25$AB - b25$SO - b25$HR))
ta25 <- c(hr = (1 - mean(b25$b))^2 * tvh25, babip = (1 - mean(b25$h))^2 * tvb25,
          cov = 2 * (1 - mean(b25$b)) * (1 - mean(b25$h)) * cov(b25$h, b25$b)) / mean(b25$f)^2
mac("taHRshare", p0(ta25["hr"] / sum(ta25))); mac("taBshare", p0(ta25["babip"] / sum(ta25)))
mac("bNoiseCur", p0(1 - tvb25 / var(b25$b))); mac("hNoiseCur", p0(1 - tvh25 / var(b25$h)))

# ---- year-to-year stability, 2004-2025 ------------------------------------------------
panel <- Batting |> filter(yearID >= 2004) |>
  group_by(playerID, yearID) |>
  summarise(across(c(AB, SO, H, HR), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
  filter(AB >= min_ab) |>
  mutate(c = (AB - SO) / AB, h = HR / (AB - SO), b = (H - HR) / (AB - SO - HR), BA = H / AB)
yrs <- sort(unique(panel$yearID))
stab <- tibble(y1 = yrs[-length(yrs)], y2 = yrs[-1]) |> filter(y2 - y1 == 1) |>
  mutate(res = map2(y1, y2, function(a, z) {
    m <- inner_join(filter(panel, yearID == a), filter(panel, yearID == z), by = "playerID", suffix = c("_1", "_2"))
    tibble(n = nrow(m), c = cor(m$c_1, m$c_2), h = cor(m$h_1, m$h_2), b = cor(m$b_1, m$b_2),
           ba = cor(m$BA_1, m$BA_2))
  })) |> unnest(res)
mac("nPairs", nrow(stab)); mac("medPairN", round(median(stab$n)))
for (v in c("c", "h", "b", "ba")) {
  nm <- c(c = "C", h = "H", b = "B", ba = "BA")[v]
  mac(paste0("med", nm), f2(median(stab[[v]])))
  mac(paste0("min", nm), f2(min(stab[[v]]))); mac(paste0("max", nm), f2(max(stab[[v]])))
}
m2425 <- inner_join(filter(panel, yearID == 2024), filter(panel, yearID == 2025), by = "playerID", suffix = c("_1", "_2"))
mac("nBoth", nrow(m2425))
mac("rCtwo", f2(cor(m2425$c_1, m2425$c_2))); mac("rHtwo", f2(cor(m2425$h_1, m2425$h_2)))
mac("rBtwo", f2(cor(m2425$b_1, m2425$b_2)))
rate_lab <- c(c = "Contact rate", h = "Home-run rate", b = "BABIP")
long2425 <- bind_rows(
  m2425 |> transmute(x = c_1, y = c_2, rate = "c"),
  m2425 |> transmute(x = h_1, y = h_2, rate = "h"),
  m2425 |> transmute(x = b_1, y = b_2, rate = "b")) |>
  mutate(rate = factor(rate, levels = c("c", "h", "b"), labels = rate_lab))
rl <- long2425 |> group_by(rate) |> summarise(r = cor(x, y))
p2f <- ggplot(long2425, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = muted, linewidth = 0.4) +
  geom_point(aes(colour = rate), size = 1.2, alpha = 0.7) +
  geom_text(data = rl, aes(x = -Inf, y = Inf, label = paste0("r = ", f2(r))),
            hjust = -0.2, vjust = 1.6, size = 3.2, colour = ink) +
  scale_colour_manual(values = unname(rate_cols)) +
  scale_x_continuous(labels = label_number(accuracy = .01)) +
  scale_y_continuous(labels = label_number(accuracy = .01)) +
  facet_wrap(~ rate, scales = "free") + labs(x = "2024 rate", y = "2025 rate")
save_fig(p2f, "fig2_repeat.pdf", 7, 2.7)

sl <- stab |> select(y1, c, h, b, ba) |> pivot_longer(-y1, names_to = "rate", values_to = "r") |>
  mutate(rate = factor(rate, levels = c("c", "h", "b", "ba"),
                       labels = c("Contact\nrate", "Home-run\nrate", "BABIP", "Batting\naverage")))
meds <- sl |> group_by(rate) |> summarise(m = median(r))
p3f <- ggplot(sl, aes(rate, r)) +
  geom_jitter(width = 0.08, height = 0, size = 1.5, alpha = 0.55, colour = ink2) +
  geom_crossbar(data = meds, aes(y = m, ymin = m, ymax = m), width = 0.45, colour = red, linewidth = 0.6) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = NULL, y = "Year-to-year correlation")
save_fig(p3f, "fig3_stability.pdf", 5.5, 3.4)

# ---- trivariate logistic-normal model, 2025 ------------------------------------------------
long3 <- bind_rows(
  b25 |> transmute(playerID, type = "c", succ = AB - SO, fail = SO),
  b25 |> transmute(playerID, type = "h", succ = HR, fail = AB - SO - HR),
  b25 |> transmute(playerID, type = "b", succ = H - HR, fail = AB - SO - H)) |>
  mutate(type = factor(type, levels = c("c", "h", "b")))
fit <- glmer(cbind(succ, fail) ~ 0 + type + (0 + type | playerID), data = long3,
             family = binomial, control = glmerControl(optimizer = "bobyqa"))
fe <- fixef(fit); se <- sqrt(diag(vcov(fit)))
Sigma <- as.matrix(VarCorr(fit)$playerID)
dimnames(Sigma) <- list(c("c", "h", "b"), c("c", "h", "b"))
sds <- sqrt(diag(Sigma)); Cor <- cov2cor(Sigma)
# conditional modes BEFORE profiling (profile() modifies the merMod in place)
re <- ranef(fit)$playerID
ci <- suppressWarnings(suppressMessages(confint(fit, parm = "theta_", method = "profile", oldNames = FALSE)))
getci <- function(pat) ci[grepl(pat, rownames(ci)), , drop = FALSE][1, ]
ci_s <- list(c = getci("^sd_typec\\|"), h = getci("^sd_typeh\\|"), b = getci("^sd_typeb\\|"))
ci_r <- list(ch = getci("cor_typeh\\.typec|cor_typec\\.typeh"),
             cb = getci("cor_typeb\\.typec|cor_typec\\.typeb"),
             hb = getci("cor_typeb\\.typeh|cor_typeh\\.typeb"))
nmL <- c(c = "C", h = "H", b = "B")
for (k in c("c", "h", "b")) {
  kk <- paste0("type", k)
  mac(paste0("mu", nmL[k]), f3(fe[[kk]]))
  mac(paste0("mu", nmL[k], "lo"), f3(fe[[kk]] - 1.96 * se[kk])); mac(paste0("mu", nmL[k], "hi"), f3(fe[[kk]] + 1.96 * se[kk]))
  mac(paste0("pmu", nmL[k]), f3(plogis(fe[[kk]])))
  mac(paste0("sig", nmL[k]), f3(sds[k]))
  mac(paste0("sig", nmL[k], "lo"), f3(ci_s[[k]][1])); mac(paste0("sig", nmL[k], "hi"), f3(ci_s[[k]][2]))
}
rho <- c(ch = Cor["c", "h"], cb = Cor["c", "b"], hb = Cor["h", "b"])
nmR <- c(ch = "CH", cb = "CB", hb = "HB")
for (k in names(rho)) {
  mac(paste0("rho", nmR[k]), f2(rho[k]))
  mac(paste0("rho", nmR[k], "lo"), f2(ci_r[[k]][1])); mac(paste0("rho", nmR[k], "hi"), f2(ci_r[[k]][2]))
}

# single-season reliabilities (delta-method noise on the logit scale; a
# continuity correction keeps a zero-HR season finite)
b25 <- b25 |> mutate(
  cs = (AB - SO + .5) / (AB + 1), hs = (HR + .5) / (AB - SO + 1), bs = (H - HR + .5) / (AB - SO - HR + 1),
  tau_c2 = 1 / (AB * cs * (1 - cs)), tau_h2 = 1 / ((AB - SO) * hs * (1 - hs)),
  tau_b2 = 1 / ((AB - SO - HR) * bs * (1 - bs)),
  w_c = sds["c"]^2 / (sds["c"]^2 + tau_c2), w_h = sds["h"]^2 / (sds["h"]^2 + tau_h2),
  w_b = sds["b"]^2 / (sds["b"]^2 + tau_b2))
for (k in c("c", "h", "b")) {
  mac(paste0("rel", nmL[k]), f2(mean(b25[[paste0("w_", k)]])))
  mac(paste0("tau", nmL[k]), f3(mean(sqrt(b25[[paste0("tau_", k, "2")]]))))
}

# shrinkage: conditional modes from lme4
b25 <- b25 |> mutate(c_s = plogis(fe[["typec"]] + re[playerID, "typec"]),
                     h_s = plogis(fe[["typeh"]] + re[playerID, "typeh"]),
                     b_s = plogis(fe[["typeb"]] + re[playerID, "typeb"]))
for (k in c("c", "h", "b"))
  mac(paste0("mnAbsD", nmL[k]), f3(mean(abs(b25[[paste0(k, "_s")]] - b25[[k]]))))
cs <- b25 |> filter(name %in% cases$name) |> arrange(desc(BA))
writeLines(tex_rows(cs |> transmute(name, f3(c), f3(c_s), f3(h), f3(h_s), f3(b), f3(b_s))), "tab_shrink.tex")
jd <- filter(b25, name == "Aaron Judge")
mac("judgeHs", f3(jd$h_s)); mac("judgeBs", f3(jd$b_s))
mac("babipShr", f3(filter(b25, name == power_pick$name)$b_s))
# how much of the BABIP shrinkage is due to the correlations: compare with
# univariate shrinkage (diagonal Sigma, normal approximation)
zfun <- function(r) c(qlogis(r$cs) - fe[["typec"]], qlogis(r$hs) - fe[["typeh"]], qlogis(r$bs) - fe[["typeb"]])
multi <- t(sapply(seq_len(nrow(b25)), function(i) {
  r <- b25[i, ]; z <- zfun(r); T <- diag(c(r$tau_c2, r$tau_h2, r$tau_b2))
  c(Sigma %*% solve(Sigma + T) %*% z, diag(diag(Sigma)) %*% solve(diag(diag(Sigma)) + T) %*% z)
}))
b_multi <- plogis(fe[["typeb"]] + multi[, 3]); b_uni <- plogis(fe[["typeb"]] + multi[, 6])
mac("mnAbsBorrow", f3(mean(abs(b_multi - b_uni))))
jdi <- which(b25$name == "Aaron Judge")
mac("judgeBuni", f3(b_uni[jdi])); mac("judgeBmulti", f3(b_multi[jdi]))

p4f <- ggplot(b25, aes(h, b)) +
  geom_segment(aes(xend = h_s, yend = b_s), colour = muted, alpha = 0.5, linewidth = 0.3,
               arrow = arrow(length = unit(0.04, "cm"))) +
  geom_point(colour = muted, size = 0.9, alpha = 0.6) +
  geom_point(aes(x = h_s, y = b_s), colour = blue, size = 1.1, alpha = 0.8) +
  geom_segment(data = cs, aes(xend = h_s, yend = b_s), colour = red, linewidth = 0.7,
               arrow = arrow(length = unit(0.14, "cm"))) +
  geom_point(data = cs, colour = ink, size = 1.8) +
  geom_text_repel(data = cs, aes(label = name), size = 3, seed = 1, min.segment.length = 0) +
  scale_x_continuous(labels = label_number(accuracy = .01)) +
  scale_y_continuous(labels = label_number(accuracy = .01)) +
  labs(x = "Home-run rate", y = "BABIP")
save_fig(p4f, "fig4_shrinkage.pdf", 6, 4.2)

# ---- historical analysis, 1901-2019 ----------------------------------------------------
mlb_lgs <- c("AL", "NL", "FL")
hist <- Batting |> filter(yearID >= 1901, yearID <= 2019, lgID %in% mlb_lgs) |>
  group_by(playerID, yearID) |>
  summarise(AB = sum(AB, na.rm = TRUE), SO = sum(SO, na.rm = TRUE), H = sum(H, na.rm = TRUE),
            HR = sum(HR, na.rm = TRUE), nNA = sum(is.na(SO)), .groups = "drop") |>
  filter(AB >= min_ab, nNA == 0, AB - SO - HR > 0) |>
  mutate(BA = H / AB, c = (AB - SO) / AB, h = HR / (AB - SO), b = (H - HR) / (AB - SO - HR),
         f = H / (AB - SO), decade = (yearID %/% 10) * 10)
mac("nHist", format(nrow(hist), big.mark = ","))
mac("histIdErr", sprintf("%.0e", max(abs(hist$BA - hist$c * (hist$h + (1 - hist$h) * hist$b)))))
q400 <- hist |> filter(BA >= .400)
mac("nFourHundred", nrow(q400)); mac("lastFourHundred", max(q400$yearID))

# sampling-noise variances on the raw scale for a proportion p from n trials
nv <- function(p, n) mean(p * (1 - p) / n)
dec <- hist |> group_by(decade) |>
  summarise(n = n(), sd_BA = sd(BA),
            var_lc = var(log(c)), var_lf = var(log(f)), cov_cf = cov(log(c), log(f)), var_lba = var(log(BA)),
            nv_lc = mean((1 - c) / (AB * c)), nv_lf = mean((1 - f) / ((AB - SO) * f)),
            r_cf = cor(c, f), r_ch = cor(c, h), r_cb = cor(c, b), r_hb = cor(h, b),
            mean_h = mean(h), mean_b = mean(b), mean_f = mean(f),
            var_c = var(c), var_h = var(h), var_b = var(b),
            cov_ch = cov(c, h), cov_cb = cov(c, b), cov_hb = cov(h, b),
            nv_c = nv(c, AB), nv_h = nv(h, AB - SO), nv_b = nv(b, AB - SO - HR),
            .groups = "drop") |>
  mutate(tv_c = var_c - nv_c, tv_h = pmax(var_h - nv_h, 1e-8), tv_b = var_b - nv_b,
         rt_ch = cov_ch / sqrt(tv_c * tv_h), rt_cb = cov_cb / sqrt(tv_c * tv_b), rt_hb = cov_hb / sqrt(tv_h * tv_b),
         tv_lc = var_lc - nv_lc, tv_lf = var_lf - nv_lf,
         a_h = (1 - mean_b)^2 * var_h / mean_f^2, a_b = (1 - mean_h)^2 * var_b / mean_f^2,
         a_cov = 2 * (1 - mean_b) * (1 - mean_h) * cov_hb / mean_f^2, a_tot = a_h + a_b + a_cov,
         gap = abs(a_tot - var_lf) / var_lf,
         ta_h = (1 - mean_b)^2 * tv_h / mean_f^2, ta_b = (1 - mean_h)^2 * tv_b / mean_f^2,
         ta_tot = ta_h + ta_b + a_cov,
         h_share = a_h / a_tot, th_share = ta_h / ta_tot, b_noise_share = nv_b / var_b)
d <- function(y) filter(dec, decade == y)
d1900 <- d(1900); d1920 <- d(1920); d1960 <- d(1960); d2010 <- d(2010)
mac("sdBAnineteen", f3(d1900$sd_BA)); mac("sdBAforties", f3(d(1940)$sd_BA)); mac("sdBAtens", f3(d2010$sd_BA))
mac("meanHnineteen", f3(d1900$mean_h)); mac("meanHtens", f3(d2010$mean_h))
mac("rchNineteen", f2(d1900$r_ch)); mac("rchTwenties", f2(d1920$r_ch)); mac("rchTens", f2(d2010$r_ch))
mac("rcbNineteen", f2(d1900$r_cb)); mac("rcbTens", f2(d2010$r_cb))
mac("rhbNineteen", f2(d1900$r_hb)); mac("rhbTens", f2(d2010$r_hb))
mac("rcfNineteen", f2(d1900$r_cf)); mac("rcfTens", f2(d2010$r_cf))
mac("rtchTens", f2(d2010$rt_ch)); mac("rtcbTens", f2(d2010$rt_cb)); mac("rthbTens", f2(d2010$rt_hb))
mac("rtcbNineteen", f2(d1900$rt_cb))
mac("hShareNineteen", p0(d1900$h_share)); mac("hShareTens", p0(d2010$h_share))
mac("thShareNineteen", p0(d1900$th_share)); mac("thShareTens", p0(d2010$th_share))
mac("bNoiseNineteen", p0(d1900$b_noise_share)); mac("bNoiseTens", p0(d2010$b_noise_share))
mac("gapMed", p1(median(dec$gap))); mac("gapMax", p1(max(dec$gap)))
mac("tvhRatio", sprintf("%.0f", d2010$tv_h / d1900$tv_h)); mac("tvbRatio", f2(d2010$tv_b / d1900$tv_b))
# which covariance channel carries the c-f covariance? cov(c, f) ~ split into HR and BABIP parts
# f = h + (1-h) b  =>  df = (1 - b) dh + (1 - h) db
dec <- dec |> mutate(cfh = (1 - mean_b) * cov_ch, cfb = (1 - mean_h) * cov_cb, cf_share_h = cfh / (cfh + cfb))
mac("cfShareHsixties", p0(d(1960)$cf_share_h))
mac("cfShareHtens", p0(filter(dec, decade == 2010)$cf_share_h))
mac("cfShareHtwenties", p0(filter(dec, decade == 1920)$cf_share_h))
# counterfactual 1 (talent scale): 2010s variances, 1900s c-f covariance
cf1 <- d2010$tv_lc + d2010$tv_lf + 2 * d1900$cov_cf
mac("talTens", f4(d2010$tv_lc + d2010$tv_lf + 2 * d2010$cov_cf))
mac("talNineteen", f4(d1900$tv_lc + d1900$tv_lf + 2 * d1900$cov_cf)); mac("talCF", f4(cf1))
# counterfactual 2: 2010s c-f covariance built with the BABIP channel only (HR channel set to zero)
cov_cf_noH <- (1 - d2010$mean_b) * 0 + (1 - d2010$mean_h) * d2010$cov_cb
cov_log_noH <- cov_cf_noH / (mean(filter(hist, decade == 2010)$c) * d2010$mean_f)
cov_log_full <- ((1 - d2010$mean_b) * d2010$cov_ch + (1 - d2010$mean_h) * d2010$cov_cb) /
  (mean(filter(hist, decade == 2010)$c) * d2010$mean_f)
mac("covLogFull", f4(cov_log_full)); mac("covLogExact", f4(d2010$cov_cf)); mac("covLogNoH", f4(cov_log_noH))
cf2 <- d2010$tv_lc + d2010$tv_lf + 2 * cov_log_noH
mac("talCFnoH", f4(cf2))
act10 <- d2010$tv_lc + d2010$tv_lf + 2 * d2010$cov_cf
mac("hrChannelShare", p0((cf2 - act10) / (cf1 - act10)))
# decade table
dt <- dec |> filter(decade %in% c(1900, 1920, 1940, 1960, 1980, 2000, 2010)) |>
  transmute(paste0(decade, "s"), n, f3(sd_BA), f3(mean_h), f2(r_cf), f2(r_ch), f2(r_cb), f2(r_hb),
            f2(rt_ch), f2(rt_cb), p0(h_share))
writeLines(tex_rows(dt), "tab_decades.tex")

# figure 5: correlations by decade (observed, with talent-scale points for the two contact channels)
cl <- dec |> select(decade, `Contact vs finishing` = r_cf, `Contact vs home runs` = r_ch, `Contact vs BABIP` = r_cb, `Home runs vs BABIP` = r_hb) |>
  pivot_longer(-decade, names_to = "pair", values_to = "r") |>
  mutate(pair = factor(pair, levels = c("Contact vs finishing", "Contact vs home runs", "Contact vs BABIP", "Home runs vs BABIP")))
lab5 <- cl |> filter(decade == 2010) |> mutate(pair_lab = as.character(pair))
p5a <- ggplot(cl, aes(decade + 5, r, colour = pair)) +
  geom_hline(yintercept = 0, colour = muted, linewidth = 0.3) +
  geom_line(aes(linetype = pair == "Contact vs finishing"), linewidth = 0.7) + geom_point(size = 1.3) +
  geom_text_repel(data = lab5, aes(label = pair_lab), direction = "y", nudge_x = 6, hjust = 0, size = 3, segment.colour = NA, seed = 1) +
  scale_colour_manual(values = c(ink2, aqua, orange, violet)) +
  scale_linetype_manual(values = c(`TRUE` = "22", `FALSE` = "solid")) +
  scale_x_continuous(breaks = seq(1900, 2020, 20), expand = expansion(mult = c(.03, .32))) +
  labs(x = "Decade", y = "Correlation across qualifying hitters")
save_fig(p5a, "fig5_tradeoffs.pdf", 6.5, 3.6)

# figure 6: delta-method attribution of Var(log f)
al <- dec |> select(decade, `Home runs` = a_h, `BABIP` = a_b, `Overlap` = a_cov) |>
  pivot_longer(-decade, names_to = "part", values_to = "v") |>
  mutate(part = factor(part, levels = c("Home runs", "BABIP", "Overlap")))
lab6 <- al |> filter(decade == 2010)
p6 <- ggplot(al, aes(decade + 5, v, colour = part)) +
  geom_hline(yintercept = 0, colour = muted, linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  geom_line(data = dec, aes(decade + 5, var_lf), inherit.aes = FALSE, colour = ink, linetype = "22", linewidth = 0.6) +
  geom_text(data = lab6, aes(label = part), hjust = -0.12, size = 3) +
  annotate("text", x = 2015, y = d2010$var_lf, label = "Total", hjust = -0.1, vjust = -0.6, size = 3, colour = ink) +
  scale_colour_manual(values = c(aqua, orange, violet)) +
  scale_x_continuous(breaks = seq(1900, 2020, 20), expand = expansion(mult = c(.03, .22))) +
  labs(x = "Decade", y = "Spread in finishing rate")
save_fig(p6, "fig6_fdecomp.pdf", 6.5, 3.6)


# ---- additions for the self-contained version ---------------------------------------
# two-factor view in 2025: c and f
r_cf_ci <- rci(b25$c, b25$f)
mac("rcfLo", f2(r_cf_ci[2])); mac("rcfHi", f2(r_cf_ci[3]))
mac("vindep", f4(var(lc) + var(lf)))
mac("covShrink", p0(1 - var(lba) / (var(lc) + var(lf))))
mac("pmuF", f3(mean(b25$f)))

# decade decomposition of Var(log BA) on the observed and talent scales
dec <- dec |> mutate(rt_cf = cov_cf / sqrt(tv_lc * tv_lf),
                     tal_lba = tv_lc + tv_lf + 2 * cov_cf,
                     nv_share_lf = nv_lf / var_lf)
mac("tvlcRatio", sprintf("%.1f", d(2010)$tv_lc / d(1900)$tv_lc))
dec10 <- filter(dec, decade == 2010); dec00 <- filter(dec, decade == 1900)
mac("tvlfRatio", f2(dec10$tv_lf / dec00$tv_lf))
mac("rtcfNineteen", f2(dec00$rt_cf)); mac("rtcfTens", f2(dec10$rt_cf))
mac("nvlfNineteen", p0(dec00$nv_share_lf)); mac("nvlfTens", p0(dec10$nv_share_lf))
mac("sdDrop", p0(1 - d(1940)$sd_BA / d(1900)$sd_BA))
rcf60 <- filter(dec, decade >= 1960)$r_cf
mac("rcfMinPost", f2(min(rcf60))); mac("rcfMaxPost", f2(max(rcf60)))
dt2 <- dec |> transmute(paste0(decade, "s"), n, f3(sd_BA), f4(var_lc), f4(var_lf), f4(2 * cov_cf),
                         f4(var_lba), f2(r_cf), f4(tv_lc), f4(tv_lf), f2(rt_cf))
writeLines(tex_rows(dt2), "tab_logdecomp.tex")

# figure 7: (a) SD of BA with .400 seasons, (b) talent variances of log c and log f, (c) talent correlation of c and f
fa <- ggplot(dec, aes(decade + 5, sd_BA)) + geom_line(colour = ink2, linewidth = 0.7) + geom_point(size = 1.3, colour = ink2) +
  geom_rug(data = q400, aes(x = yearID), inherit.aes = FALSE, colour = red, length = unit(0.05, "npc")) +
  labs(x = NULL, y = NULL, title = "(a) SD of batting average")
tvl <- dec |> select(decade, `log c` = tv_lc, `log f` = tv_lf) |> pivot_longer(-decade)
fb <- ggplot(tvl, aes(decade + 5, value, colour = name)) + geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  geom_text(data = filter(tvl, decade == 2010), aes(label = name), hjust = -0.2, size = 3) +
  scale_colour_manual(values = c(blue, orange)) +
  scale_x_continuous(expand = expansion(mult = c(.03, .18))) +
  labs(x = NULL, y = NULL, title = "(b) Talent variance")
fc <- ggplot(dec, aes(decade + 5)) + geom_hline(yintercept = 0, colour = muted, linewidth = 0.3) +
  geom_line(aes(y = r_cf), colour = muted, linetype = "22", linewidth = 0.6) +
  geom_line(aes(y = rt_cf), colour = ink, linewidth = 0.7) + geom_point(aes(y = rt_cf), size = 1.3, colour = ink) +
  labs(x = NULL, y = NULL, title = "(c) Correlation of c and f")
suppressMessages({
  pdf(file.path("figures", "fig7_history.pdf"), width = 7.5, height = 2.6)
  gridExtra::grid.arrange(fa, fb, fc, nrow = 1, bottom = "Decade")
  dev.off()
})

# ---- write macros -----------------------------------------------------------------
mac("lahmanVersion", as.character(packageVersion("Lahman")))
writeLines(c("% Generated by analysis.R -- do not edit by hand", macros), "numbers.tex")
cat("Wrote", length(macros), "macros\n")
print(dec |> select(decade, n, sd_BA, var_lc, var_lf, cov_cf, tv_lc, tv_lf, rt_cf, mean_h, r_cf, r_ch, r_cb, r_hb, rt_ch, rt_cb, rt_hb, h_share, th_share, gap, cf_share_h, b_noise_share), n = 20, width = 200)
print(Sigma); print(Cor); print(ci)
