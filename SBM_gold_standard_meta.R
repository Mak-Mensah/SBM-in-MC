############################################################
## GOLD-STANDARD META-ANALYSIS SCRIPT
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(clubSandwich)
  library(janitor)
  library(readr)
  library(fs)
  library(glue)
  library(flextable)
  library(officer)
  library(stringr)
  library(patchwork)
  library(ggplot2)
  library(grid)
})

## =========================================================
## 1) PATHS
## =========================================================
input_csv  <- "SBMdata.csv"
## input_csv <- "/mnt/data/SBMdata.csv"  # if needed

output_dir <- "meta_outputs_sbm"
dir_create(output_dir)

if (!file.exists(input_csv)) {
  stop("Cannot find '", input_csv, "'. Put it in your working directory or set input_csv to a full path.")
}

## =========================================================
## 2) HELPERS
## =========================================================
pct <- function(x) (exp(x) - 1) * 100

sig_code <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

fmt_ci <- function(est, se) {
  lo <- est - 1.96 * se
  hi <- est + 1.96 * se
  glue("{round(est,3)} [{round(lo,3)}, {round(hi,3)}]")
}

fmt_pct_ci <- function(est_lnrr, se_lnrr) {
  lo <- est_lnrr - 1.96 * se_lnrr
  hi <- est_lnrr + 1.96 * se_lnrr
  glue("{sprintf('%+.1f', pct(est_lnrr))}% [{sprintf('%+.1f', pct(lo))}, {sprintf('%+.1f', pct(hi))}]")
}

z_na <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_file <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\]", "-", x)
  x <- gsub("[:*?\"<>|]", "", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("_+", "_", x)
  trimws(x)
}

theme_pub <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.95)),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

parse_num <- function(x) readr::parse_number(as.character(x), locale = locale(grouping_mark = ","))

fit_mv <- function(df, mods = NULL) {
  if (is.null(mods)) {
    rma.mv(yi, vi,
           random = ~ 1 | study_id/comparison_id,
           data = df, method = "REML")
  } else {
    rma.mv(yi, vi,
           mods = mods,
           random = ~ 1 | study_id/comparison_id,
           data = df, method = "REML")
  }
}

cr2_test <- function(m, cluster) {
  clubSandwich::coef_test(m, cluster = cluster, vcov = "CR2")
}

i2_mv_approx <- function(m) {
  tau2 <- sum(m$sigma2, na.rm = TRUE)
  100 * tau2 / (tau2 + mean(m$vi, na.rm = TRUE))
}

tau2_total <- function(m) sum(m$sigma2, na.rm = TRUE)

R2_tau <- function(m_null, m_mod) {
  t0 <- tau2_total(m_null)
  t1 <- tau2_total(m_mod)
  if (is.na(t0) || t0 <= 0) return(NA_real_)
  pmax(0, (t0 - t1) / t0)
}

save_docx_table <- function(df, title, path) {
  ft <- flextable(df) |>
    autofit() |>
    align(align = "center", part = "all") |>
    add_header_lines(values = title)
  save_as_docx(setNames(list(ft), title), path = path)
}

plot_residual_funnel <- function(resid, sei, title, subtitle = NULL) {
  d <- tibble(resid = resid, sei = sei) |> filter(is.finite(resid), is.finite(sei))
  if (nrow(d) < 5) return(NULL)
  
  ggplot(d, aes(x = resid, y = sei)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6) +
    scale_y_reverse() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Residuals (yi - fitted)",
      y = "SE (sqrt(vi), reversed)"
    ) +
    theme_pub(13)
}

pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf

save_plot_robust <- function(p, filename, width = 9, height = 6, dpi = 450, bg = "white") {
  out <- file.path(output_dir, filename)
  ok <- tryCatch({
    ggsave(out, plot = p, width = width, height = height, dpi = dpi, bg = bg)
    TRUE
  }, error = function(e) {
    message("ggsave failed for ", filename, ": ", conditionMessage(e))
    FALSE
  })
  
  if (!ok) {
    tryCatch({
      png(out, width = width * dpi, height = height * dpi, res = dpi)
      print(p)
      dev.off()
      message("Saved via png() fallback: ", filename)
    }, error = function(e) {
      message("png() fallback failed for ", filename, ": ", conditionMessage(e))
      try(dev.off(), silent = TRUE)
      stop("Failed to save plot: ", filename)
    })
  }
}

fit_mv_stats <- function(d, fml) {
  m0 <- fit_mv(d)
  m1 <- fit_mv(d, mods = fml)
  rob <- cr2_test(m1, d$study_id)
  
  aic0 <- tryCatch(AIC(m0), error = function(e) NA_real_)
  aic1 <- tryCatch(AIC(m1), error = function(e) NA_real_)
  
  list(
    m0 = m0,
    m1 = m1,
    rob = rob,
    AIC_null = aic0,
    AIC_mod  = aic1,
    Delta_AIC = aic1 - aic0,
    R2_tau = R2_tau(m0, m1),
    tau2_null = tau2_total(m0),
    tau2_mod  = tau2_total(m1)
  )
}

ensure_complete <- function(df, mods, min_k = 10, min_studies = 5) {
  mods <- mods[mods %in% names(df)]
  if (length(mods) == 0) stop("No moderators supplied (none found in df).")
  
  while (length(mods) > 0) {
    d <- df |> dplyr::filter(dplyr::if_all(dplyr::all_of(mods), ~ !is.na(.x)))
    if (nrow(d) >= min_k && dplyr::n_distinct(d$study_id) >= min_studies) {
      return(list(mods = mods, data = d))
    }
    na_rate <- vapply(mods, function(m) mean(is.na(df[[m]])), numeric(1))
    if (all(is.na(na_rate)) || length(na_rate) == 0) break
    drop_one <- names(which.max(na_rate))[1]
    miss_pct <- round(100 * max(na_rate, na.rm = TRUE), 1)
    message(sprintf("Dropping '%s' (missingness=%.1f%%) to retain enough complete cases.",
                    drop_one, miss_pct))
    mods <- setdiff(mods, drop_one)
  }
  
  stop("No set of moderators yields enough complete cases. Use fewer moderators or fit separate models.")
}

drop_high_cor <- function(df, vars, cutoff = 0.70) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) < 2) return(vars)
  M <- df |> select(all_of(vars)) |> as.matrix()
  cm <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  
  drop <- character(0)
  while (TRUE) {
    cm2 <- cm; diag(cm2) <- 0
    if (all(is.na(cm2)) || max(abs(cm2), na.rm = TRUE) < cutoff) break
    idx <- which(abs(cm2) == max(abs(cm2), na.rm = TRUE), arr.ind = TRUE)[1, ]
    a <- colnames(cm2)[idx[1]]; b <- colnames(cm2)[idx[2]]
    na_a <- sum(is.na(df[[a]])); na_b <- sum(is.na(df[[b]]))
    drop_one <- if (na_a > na_b) a else if (na_b > na_a) b else b
    drop <- union(drop, drop_one)
    keep <- setdiff(colnames(cm), drop)
    if (length(keep) < 2) break
    cm <- cm[keep, keep, drop = FALSE]
  }
  setdiff(vars, drop)
}

read_plot_rds <- function(path) {
  if (!file.exists(path)) stop("Missing panel RDS: ", path)
  readRDS(path)
}

vif_manual <- function(X) {
  X <- as.matrix(X)
  v <- rep(NA_real_, ncol(X))
  names(v) <- colnames(X)
  for (j in seq_len(ncol(X))) {
    yj <- X[, j]
    xj <- X[, -j, drop = FALSE]
    if (sd(yj, na.rm = TRUE) == 0) { v[j] <- NA_real_; next }
    fit <- lm(yj ~ xj)
    r2 <- summary(fit)$r.squared
    v[j] <- 1 / (1 - r2)
  }
  tibble(Term = names(v), VIF = as.numeric(v)) |> arrange(desc(VIF))
}

collapse_for_mv <- function(df, var, topN = 8, min_n = 2, other = "Other", missing = "Unknown") {
  if (!var %in% names(df)) return(df)
  x <- as.character(df[[var]])
  x[is.na(x) | x == ""] <- missing
  x <- factor(x)
  x <- forcats::fct_lump_n(x, n = topN, other_level = other)
  x <- forcats::fct_lump_min(x, min = min_n, other_level = other)
  df[[var]] <- droplevels(x)
  df
}

fixed_df_burden <- function(df, mods) {
  mods <- mods[mods %in% names(df)]
  if (length(mods) == 0) return(0L)
  sum(vapply(mods, function(v) {
    if (is.factor(df[[v]])) max(0L, nlevels(df[[v]]) - 1L) else 1L
  }, integer(1)))
}

fit_mv_safe <- function(df, mods, random = ~ 1 | study_id/comparison_id, method = "REML") {
  tries <- list(
    list(control = list(optimizer = "nlminb", eval.max = 20000, iter.max = 20000)),
    list(control = list(optimizer = "optim", optmethod = "BFGS", maxit = 200000)),
    list(control = list(optimizer = "optim", optmethod = "Nelder-Mead", maxit = 200000))
  )
  
  for (i in seq_along(tries)) {
    m <- tryCatch(
      metafor::rma.mv(yi, vi, mods = mods, random = random, data = df, method = method,
                      control = tries[[i]]$control),
      error = function(e) NULL
    )
    if (!is.null(m)) return(m)
  }
  
  for (i in seq_along(tries)) {
    m <- tryCatch(
      metafor::rma.mv(yi, vi, mods = mods, random = ~ 1 | study_id, data = df, method = method,
                      control = tries[[i]]$control),
      error = function(e) NULL
    )
    if (!is.null(m)) return(m)
  }
  
  NULL
}

fit_mv_safe2 <- function(df, fml, method = "REML") {
  m <- fit_mv_safe(df, mods = fml, method = method)
  if (!is.null(m)) return(m)
  tryCatch(fit_mv_safe(df, mods = fml, method = "ML"), error = function(e) NULL)
}

## ---- NULL-safe patchwork wrapper ----
wrap_plots_safe <- function(..., ncol = 1, tag_levels = "A") {
  plots <- list(...)
  plots <- purrr::compact(plots)  # drop NULLs
  if (length(plots) == 0) return(NULL)
  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(tag_levels = tag_levels)
}

save_rds_png <- function(p, stub, width = 9.5, height = 6.0) {
  if (is.null(p)) return(invisible(FALSE))
  saveRDS(p, file.path(output_dir, paste0(stub, ".rds")))
  ggsave(file.path(output_dir, paste0(stub, ".png")),
         p, width = width, height = height, dpi = 600, bg = "white")
  invisible(TRUE)
}

read_rds_stub_safe <- function(stub) {
  path <- fs::path(output_dir, paste0(stub, ".rds"))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

save_composite <- function(p, out_stub, width, height) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(fs::path(output_dir, paste0(out_stub, ".pdf")),
         plot = p, width = width, height = height, units = "in",
         device = pdf_device)
  ggsave(fs::path(output_dir, paste0(out_stub, ".png")),
         plot = p, width = width, height = height, units = "in",
         dpi = 600, bg = "white")
}

## =========================================================
## 3) READ + CLEAN (SBM) + MAP TO REQUIRED NAMES
## =========================================================
dat_raw <- read_csv(input_csv, locale = locale(encoding = "latin1"),
                    show_col_types = FALSE) |>
  clean_names() |>
  mutate(across(where(is.character), ~na_if(.x, "")))

dat_raw <- dat_raw |>
  rename(
    outcome_variable = outcome,
    treatment_mean   = exp_mean,
    treatment_sd     = exp_sd,
    treatment_n      = exp_n,
    control_mean     = ctrl_mean,
    control_sd       = ctrl_sd,
    control_n        = ctrl_n
  )

req <- c("study_id","outcome_variable",
         "treatment_mean","treatment_sd","treatment_n",
         "control_mean","control_sd","control_n")
miss <- setdiff(req, names(dat_raw))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

num_cols <- intersect(
  c("initial_ec_d_sm","initial_p_h","duration_years",
    "treatment_mean","treatment_sd","treatment_n",
    "control_mean","control_sd","control_n"),
  names(dat_raw)
)

dat <- dat_raw |>
  mutate(across(all_of(num_cols), parse_num)) |>
  mutate(
    study_id = factor(study_id),
    
    outcome_type = factor(ifelse(is.na(outcome_variable) | outcome_variable == "",
                                 "Unknown", outcome_variable)),
    
    ## REQUESTED: group equals outcome_type (Yield/SOC/Enzyme)
    group = droplevels(outcome_type),
    
    study_type = if ("experiment_type" %in% names(dat_raw)) factor(experiment_type) else factor("Unknown"),
    comparison_id = row_number()
  ) |>
  filter(
    treatment_mean > 0, control_mean > 0,
    treatment_n >= 2, control_n >= 2,
    treatment_sd > 0, control_sd > 0
  )

write_csv(dat, file.path(output_dir, "data_after_cleaning.csv"))

## =========================================================
## 4) EFFECT SIZES (lnRR via ROM) + sei + z-scores
## =========================================================
dat_es <- escalc(
  measure = "ROM",
  m1i = treatment_mean, sd1i = treatment_sd, n1i = treatment_n,
  m2i = control_mean,   sd2i = control_sd,   n2i = control_n,
  data = dat
) |>
  as_tibble() |>
  mutate(sei = sqrt(vi))

## Ensure key categorical moderators are factors
key_fac_mods <- c("climate", "soil_texture", "crop", "cropping_system", "irrigation_type", "fertilizer")
for (v in intersect(key_fac_mods, names(dat_es))) {
  dat_es[[v]] <- dat_es[[v]] |>
    as.character() |>
    na_if("") |>
    replace_na("Unknown") |>
    factor()
}

## Climate × Fertilizer interaction for requested Fig *A
if (all(c("climate","fertilizer") %in% names(dat_es))) {
  dat_es <- dat_es |>
    mutate(
      climate_x_fertilizer = interaction(climate, fertilizer, sep = " × ", drop = TRUE) |> droplevels()
    )
}

## z-scores for numeric moderators (incl initial EC/pH)
exclude_from_z <- c(
  "yi","vi","sei","zi","pval","ci_lb","ci_ub",
  "treatment_mean","treatment_sd","treatment_n",
  "control_mean","control_sd","control_n"
)
num_mods <- names(dat_es)[sapply(dat_es, is.numeric)]
num_mods <- setdiff(num_mods, exclude_from_z)
for (nm in num_mods) dat_es[[paste0(nm, "_z")]] <- z_na(dat_es[[nm]])

write_csv(dat_es, file.path(output_dir, "effect_sizes_lnRR.csv"))

## =========================================================
## 4C) REQUESTED PANELS: helper plot builders
## =========================================================
make_factor_plot_outcome <- function(d, fvar, outcome_label,
                                     min_k_total = 10, min_studies = 5, min_k_per_level = 3,
                                     topN = 10) {
  if (!fvar %in% names(d)) return(NULL)
  d <- d |> filter(!is.na(.data[[fvar]]), is.finite(yi), is.finite(vi))
  if (nrow(d) < min_k_total || n_distinct(d$study_id) < min_studies) return(NULL)
  
  lvl_counts <- d |>
    count(.data[[fvar]], name = "k_level") |>
    mutate(level_chr = as.character(.data[[fvar]])) |>
    filter(!is.na(level_chr), level_chr != "") |>
    filter(k_level >= min_k_per_level) |>
    arrange(desc(k_level)) |>
    slice_head(n = topN)
  
  keep_lvls <- lvl_counts |> pull(level_chr)
  d <- d |> filter(as.character(.data[[fvar]]) %in% keep_lvls)
  if (n_distinct(d[[fvar]]) < 2) return(NULL)
  
  fml <- as.formula(paste0("~ 0 + ", fvar))
  m <- tryCatch(
    rma.mv(yi, vi, mods = fml,
           random = ~ 1 | study_id/comparison_id,
           data = d, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  rob <- tryCatch(cr2_test(m, d$study_id), error = function(e) NULL)
  if (is.null(rob)) return(NULL)
  
  tt <- tibble(
    term = rownames(rob),
    estimate = rob$beta,
    se = rob$SE,
    p = rob$p_Satt
  ) |>
    mutate(
      level = str_replace(term, paste0("^", fvar), ""),
      level = str_replace_all(level, "^\\s*|\\s*$", ""),
      lo = estimate - 1.96*se,
      hi = estimate + 1.96*se,
      pct_est = pct(estimate),
      pct_lo  = pct(lo),
      pct_hi  = pct(hi),
      sig = sig_code(p)
    ) |>
    arrange(pct_est)
  
  ggplot(tt, aes(x = pct_est, y = reorder(level, pct_est))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.8) +
    geom_errorbarh(aes(xmin = pct_lo, xmax = pct_hi), height = 0.22, linewidth = 0.6) +
    labs(
      title = glue("{outcome_label}: {fvar} (multilevel; CR2)"),
      subtitle = glue("Top-{nrow(lvl_counts)} levels (k per level ≥ {min_k_per_level}); k={nrow(d)}, studies={n_distinct(d$study_id)}."),
      x = "Percent change vs control",
      y = NULL
    ) +
    theme_pub(12)
}

make_bubble_plot_outcome <- function(d, xvar, outcome_label,
                                     min_k_total = 10, min_studies = 5) {
  if (!xvar %in% names(d)) return(NULL)
  d <- d |> filter(is.finite(yi), is.finite(vi), is.finite(.data[[xvar]]))
  if (nrow(d) < min_k_total || n_distinct(d$study_id) < min_studies) return(NULL)
  
  fml <- as.formula(paste0("~ ", xvar))
  stats <- tryCatch(fit_mv_stats(d, fml), error = function(e) NULL)
  if (is.null(stats)) return(NULL)
  rob <- stats$rob
  
  idx <- match(xvar, rownames(rob))
  if (is.na(idx)) return(NULL)
  
  slope_est <- rob$beta[idx]
  slope_se  <- rob$SE[idx]
  slope_p   <- rob$p_Satt[idx]
  
  lab <- glue(
    "CR2 slope: β = {round(slope_est, 4)} (SE = {round(slope_se, 4)}), p = {signif(slope_p, 3)}\n",
    "R²(τ²) = {ifelse(is.na(stats$R2_tau), 'NA', sprintf('%.3f', stats$R2_tau))}; ΔAIC = {ifelse(is.na(stats$Delta_AIC), 'NA', sprintf('%.2f', stats$Delta_AIC))}\n",
    "k = {nrow(d)}, studies = {n_distinct(d$study_id)}"
  )
  
  ggplot(d, aes(x = .data[[xvar]], y = yi, size = 1/vi)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.65) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
    annotate("text", x = Inf, y = Inf, label = lab,
             hjust = 1.05, vjust = 1.1, size = 3.4, fontface = "italic") +
    labs(
      title = glue("{outcome_label}: lnRR vs {xvar}"),
      subtitle = "Point size = inverse-variance weight (1/vi). Line = fitted trend with 95% CI band.",
      x = xvar,
      y = "ln Response Ratio (lnRR)",
      size = "Weight (1/vi)"
    ) +
    theme_pub(13)
}

## =========================================================
## 4D) BUILD REQUESTED FIGS + PANELS (NULL-safe)
## =========================================================
outcomes_requested <- c("Yield","SOC","Enzyme")
outcomes_present <- intersect(outcomes_requested, levels(dat_es$outcome_type))

for (ot in outcomes_present) {
  d_ot <- dat_es |> filter(outcome_type == ot)
  
  ## A: Climate × Fertilizer
  pA <- make_factor_plot_outcome(d_ot, "climate_x_fertilizer", outcome_label = ot, topN = 12)
  save_rds_png(pA, glue("REQ_{ot}_A_climate_x_fertilizer"), width = 10.5, height = 6.2)
  
  ## B: EC/pH vs outcome (Yield/SOC: both EC + pH; Enzyme: pH only)
  if (ot %in% c("Yield","SOC")) {
    pB1 <- make_bubble_plot_outcome(d_ot, "initial_ec_d_sm_z", outcome_label = ot)
    pB2 <- make_bubble_plot_outcome(d_ot, "initial_p_h_z",     outcome_label = ot)
    save_rds_png(pB1, glue("REQ_{ot}_B1_ec"), width = 9.2, height = 5.6)
    save_rds_png(pB2, glue("REQ_{ot}_B2_ph"), width = 9.2, height = 5.6)
    
    pB <- wrap_plots_safe(
      read_rds_stub_safe(glue("REQ_{ot}_B1_ec")),
      read_rds_stub_safe(glue("REQ_{ot}_B2_ph")),
      ncol = 1,
      tag_levels = "a"
    )
    if (!is.null(pB)) {
      saveRDS(pB, file.path(output_dir, glue("REQ_{ot}_B_EC_pH.rds")))
      save_composite(pB, glue("REQ_{ot}_B_EC_pH"), width = 10.5, height = 12)
    }
  } else if (ot == "Enzyme") {
    pB <- make_bubble_plot_outcome(d_ot, "initial_p_h_z", outcome_label = ot)
    save_rds_png(pB, glue("REQ_{ot}_B_ph"), width = 10.5, height = 6.2)
  }
  
  ## C: outcome-specific
  if (ot == "Yield") {
    pC1 <- make_factor_plot_outcome(d_ot, "cropping_system", outcome_label = ot, topN = 12)
    pC2 <- make_factor_plot_outcome(d_ot, "irrigation_type", outcome_label = ot, topN = 12)
    save_rds_png(pC1, glue("REQ_{ot}_C1_cropping"), width = 10.5, height = 6.2)
    save_rds_png(pC2, glue("REQ_{ot}_C2_irrigation"), width = 10.5, height = 6.2)
    
    pC <- wrap_plots_safe(
      read_rds_stub_safe(glue("REQ_{ot}_C1_cropping")),
      read_rds_stub_safe(glue("REQ_{ot}_C2_irrigation")),
      ncol = 1,
      tag_levels = "a"
    )
    if (!is.null(pC)) {
      saveRDS(pC, file.path(output_dir, glue("REQ_{ot}_C_Cropping_Irrigation.rds")))
      save_composite(pC, glue("REQ_{ot}_C_Cropping_Irrigation"), width = 10.5, height = 12)
    }
  }
  
  if (ot == "SOC") {
    pC <- make_factor_plot_outcome(d_ot, "cropping_system", outcome_label = ot, topN = 12)
    save_rds_png(pC, glue("REQ_{ot}_C_cropping_system"), width = 10.5, height = 6.2)
  }
  
  if (ot == "Enzyme") {
    pC <- make_factor_plot_outcome(d_ot, "fertilizer", outcome_label = ot, topN = 12)
    save_rds_png(pC, glue("REQ_{ot}_C_fertilizer"), width = 10.5, height = 6.2)
  }
}

## ---- Assemble Panels (NULL-safe; never crashes) ----
if ("Yield" %in% outcomes_present) {
  p1A <- read_rds_stub_safe("REQ_Yield_A_climate_x_fertilizer")
  p1B <- read_rds_stub_safe("REQ_Yield_B_EC_pH")
  p1C <- read_rds_stub_safe("REQ_Yield_C_Cropping_Irrigation")
  
  Panel_1_Yield <- wrap_plots_safe(p1A, p1B, p1C, ncol = 1, tag_levels = "A")
  save_composite(Panel_1_Yield, "Panel_1_Yield", width = 11, height = 19)
}

if ("SOC" %in% outcomes_present) {
  p2A <- read_rds_stub_safe("REQ_SOC_A_climate_x_fertilizer")
  p2B <- read_rds_stub_safe("REQ_SOC_B_EC_pH")
  p2C <- read_rds_stub_safe("REQ_SOC_C_cropping_system")
  
  Panel_2_SOC <- wrap_plots_safe(p2A, p2B, p2C, ncol = 1, tag_levels = "A")
  save_composite(Panel_2_SOC, "Panel_2_SOC", width = 11, height = 19)
}

if ("Enzyme" %in% outcomes_present) {
  p3A <- read_rds_stub_safe("REQ_Enzyme_A_climate_x_fertilizer")
  p3B <- read_rds_stub_safe("REQ_Enzyme_B_ph")
  p3C <- read_rds_stub_safe("REQ_Enzyme_C_fertilizer")
  
  Panel_3_Enzyme <- wrap_plots_safe(p3A, p3B, p3C, ncol = 1, tag_levels = "A")
  save_composite(Panel_3_Enzyme, "Panel_3_Enzyme", width = 11, height = 19)
}

## =========================================================
## 5) OVERALL MULTILEVEL MODELS (by outcome_type)
## =========================================================
overall_outcome_tbl <- tibble()
for (ot in levels(dat_es$outcome_type)) {
  d <- dat_es |> filter(outcome_type == ot)
  if (nrow(d) < 3) next
  m <- fit_mv(d)
  rob <- cr2_test(m, d$study_id)
  overall_outcome_tbl <- bind_rows(
    overall_outcome_tbl,
    tibble(
      Outcome_Type = as.character(ot),
      k = nrow(d),
      Studies = n_distinct(d$study_id),
      lnRR = as.numeric(coef(m)),
      CI_low = m$ci.lb,
      CI_high = m$ci.ub,
      Percent = pct(as.numeric(coef(m))),
      p_model = m$pval,
      p_CR2 = rob$p_Satt[1],
      tau2_total = sum(m$sigma2),
      I2_approx = i2_mv_approx(m)
    )
  )
}
write_csv(overall_outcome_tbl, file.path(output_dir, "table_overall_by_outcome_type.csv"))

## =========================================================
## 6) PUBLICATION BIAS (study-aggregated by outcome)
## =========================================================
agg <- dat_es |>
  group_by(group, study_id) |>
  summarise(
    yi = mean(yi, na.rm = TRUE),
    vi = mean(vi, na.rm = TRUE),
    sei = sqrt(mean(vi, na.rm = TRUE)),
    .groups = "drop"
  )

egger_tbl    <- tibble()
begg_tbl     <- tibble()
trimfill_tbl <- tibble()

for (g in unique(agg$group)) {
  dg <- agg |> filter(group == g)
  if (nrow(dg) < 5) next
  
  dg2 <- dg |> arrange(yi)
  m2  <- rma(yi, vi, data = dg2, method = "REML")
  
  tau2 <- as.numeric(m2$tau2)
  I2   <- as.numeric(m2$I2)
  QE_p <- as.numeric(m2$QEp)
  
  eg <- NULL
  if (nrow(dg2) >= 10) {
    eg <- tryCatch(regtest(m2, model = "rma", predictor = "sei"),
                   error = function(e) NULL)
  }
  eg_z <- if (!is.null(eg)) as.numeric(eg$zval) else NA_real_
  eg_p <- if (!is.null(eg)) as.numeric(eg$pval) else NA_real_
  eg_note <- if (nrow(dg2) < 10) "Egger not run (k<10)" else if (is.na(eg_p)) "Egger failed" else glue("Egger: z={round(eg_z,2)}, p={signif(eg_p,3)}")
  
  egger_tbl <- bind_rows(
    egger_tbl,
    tibble(Outcome = as.character(g), k = nrow(dg2), tau2 = tau2, I2 = I2, QEp = QE_p,
           Egger_z = eg_z, Egger_p = eg_p, Egger_note = eg_note)
  )
  
  bg <- tryCatch(ranktest(m2), error = function(e) NULL)
  bg_tau <- if (!is.null(bg)) as.numeric(bg$tau) else NA_real_
  bg_p   <- if (!is.null(bg)) as.numeric(bg$pval) else NA_real_
  bg_note <- if (is.na(bg_p)) "Begg failed" else glue("Begg: tau={round(bg_tau,2)}, p={signif(bg_p,3)}")
  
  begg_tbl <- bind_rows(
    begg_tbl,
    tibble(Outcome = as.character(g), k = nrow(dg2), Begg_tau = bg_tau, Begg_p = bg_p, Begg_note = bg_note)
  )
  
  taf <- tryCatch(trimfill(m2), error = function(e) NULL)
  if (!is.null(taf)) {
    k0 <- tryCatch(as.numeric(taf$k0), error = function(e) NA_real_)
    est0 <- as.numeric(coef(m2))
    est1 <- as.numeric(coef(taf))
    
    trimfill_tbl <- bind_rows(
      trimfill_tbl,
      tibble(
        Outcome = as.character(g),
        k_observed = nrow(dg2),
        k0_imputed = k0,
        est_lnRR_original = est0,
        est_lnRR_trimfill = est1,
        delta_lnRR = est1 - est0
      )
    )
  }
}

write_csv(egger_tbl,    file.path(output_dir, "table_egger_tests_by_outcome.csv"))
write_csv(begg_tbl,     file.path(output_dir, "table_begg_ranktest_by_outcome.csv"))
write_csv(trimfill_tbl, file.path(output_dir, "table_trimfill_by_outcome.csv"))

## =========================================================
## 6E) ORCHARD PLOTS (WITH intercept; recommended)
## =========================================================
orch_log <- file.path(output_dir, "orchard_log.txt")
cat("ORCHARD LOG\n", file = orch_log)

if (!requireNamespace("orchaRd", quietly = TRUE)) {
  cat("Package orchaRd NOT installed; orchard plots skipped.\n", file = orch_log, append = TRUE)
} else {
  cat("Package orchaRd detected.\n", file = orch_log, append = TRUE)
  
  dat_es <- dat_es |>
    mutate(
      study_id = droplevels(factor(study_id)),
      group = droplevels(factor(group))
    )
  
  m_out <- tryCatch(
    metafor::rma.mv(
      yi, vi,
      mods   = ~ group,  # intercept included
      random = ~ 1 | study_id/comparison_id,
      data   = dat_es,
      method = "REML"
    ),
    error = function(e) {
      cat("m_out FAILED: ", conditionMessage(e), "\n", file = orch_log, append = TRUE)
      NULL
    }
  )
  
  if (!is.null(m_out)) {
    p_orch <- tryCatch(
      orchaRd::orchard_plot(
        object = m_out,
        mod    = "group",
        group  = "study_id",
        xlab   = "ln Response Ratio (lnRR)"
      ) + theme_pub(13) + labs(title = "Orchard plot – Outcome effects"),
      error = function(e) {
        cat("orchard_plot FAILED: ", conditionMessage(e), "\n", file = orch_log, append = TRUE)
        NULL
      }
    )
    
    if (!is.null(p_orch)) {
      save_plot_robust(p_orch, "figure_orchard_outcomes.png", width = 9.2, height = 6.0, dpi = 450)
      cat("Saved: figure_orchard_outcomes.png\n", file = orch_log, append = TRUE)
    }
  }
  
  cat("Done.\n", file = orch_log, append = TRUE)
}

## =========================================================
## 10) SOFTWARE VERSIONS + SESSION INFO
## =========================================================
pkgs_used <- c("tidyverse", "metafor", "clubSandwich", "janitor", "readr", "fs", "glue",
               "flextable", "officer", "stringr", "patchwork", "ggplot2")
pkg_tbl <- tibble(
  R = R.version.string,
  Package = pkgs_used,
  Version = map_chr(pkgs_used, ~ as.character(packageVersion(.x)))
)
write_csv(pkg_tbl, file.path(output_dir, "table_software_versions.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))

message("✅ Meta-analysis completed successfully.")
message("✅ Outputs saved in: ", output_dir)
message("✅ Panels saved (when possible): Panel_1_Yield, Panel_2_SOC, Panel_3_Enzyme (PDF + PNG).")
message("✅ Orchard plot attempted; see orchard_log.txt for details.")
