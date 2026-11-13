# ============================================================
# BIOCHAR META-ANALYSIS — INTEGRATED ANALYSIS + FIGURES
# BiocharDS V1.0 | Maize and wheat
#
# Outcomes: N2O and GWP only
# SOC_initial remains a continuous moderator.
#
# Primary model: multilevel REML
# Primary inference: article-clustered CR2-Satterthwaite
# Moderator models: one moderator at a time
# Nested fixed-effect comparisons: ML
#
# Figures: PNG, PDF, TIFF (600 dpi, LZW)
# Supplementary tables: S1-S8 (S8a/S8b panels) + Excel workbook
# ============================================================

# -------------------- 1. Packages ----------------------------

packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "tibble",
  "purrr", "metafor", "clubSandwich", "ggplot2", "patchwork",
  "scales", "openxlsx"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Install required packages first:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

# -------------------- 2. Configuration ----------------------

input_xlsx <- if (file.exists("BiocharDS_V1.0 Maize and wheat.xlsx")) {
  "BiocharDS_V1.0 Maize and wheat.xlsx"
} else if (file.exists("BiocharDS_V1.0 Maize and wheat(2).xlsx")) {
  "BiocharDS_V1.0 Maize and wheat(2).xlsx"
} else {
  "BiocharDS_V1.0 Maize and wheat.xlsx"
}
input_sheet <- "Sheet1"
output_root <- "biochar_meta_outputs"

# Deliberately exclude SOC as an outcome.
outcomes_requested <- c("N2O", "GWP")

outcome_labels <- c(
  N2O = "N2O emissions",
  GWP = "Global warming potential (GWP)"
)

confidence_level <- 0.95
alpha <- 1 - confidence_level
low_cr2_df_threshold <- 4

min_effects_moderator <- 10L
min_articles_moderator <- 5L
min_unique_continuous_values <- 3L
min_effects_per_category <- 5L
min_articles_per_category <- 3L

min_effects_egger <- 10L
min_articles_egger <- 10L
min_articles_trimfill <- 10L

expected_feedstocks <- c("Wheat straw", "Maize straw")

main_text_continuous <- c(
  "pyrolysis_temp",
  "biochar_rate_log10",
  "trial_duration_log10",
  "ph_initial"
)

main_text_categorical <- c(
  "feedstock", "texture", "climate_zone", "treatment_type"
)

main_text_moderators <- c(
  main_text_continuous, main_text_categorical
)

# -------------------- 3. Helpers -----------------------------

lnrr_to_pct <- function(x) (exp(x) - 1) * 100

safe_num <- function(x) suppressWarnings(as.numeric(x))

clean_chr <- function(x) {
  y <- stringr::str_squish(as.character(x))
  y[y %in% c("", "NA", "N/A", "na", "n/a")] <- NA_character_
  y
}

clean_texture <- function(x) {
  clean_chr(x) |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("-", "_") |>
    stringr::str_replace_all("\\s+", "_")
}

clean_feedstock <- function(x) {
  y <- clean_chr(x)
  unexpected <- setdiff(unique(na.omit(y)), expected_feedstocks)
  
  if (length(unexpected)) {
    stop(
      "Unexpected Type_biochar value(s): ",
      paste(unexpected, collapse = ", "),
      ". Expected: ",
      paste(expected_feedstocks, collapse = "; ")
    )
  }
  
  y
}

parse_num_mid <- function(x) {
  x <- clean_chr(x)
  
  vapply(x, function(s) {
    if (is.na(s)) return(NA_real_)
    z <- stringr::str_extract_all(s, "[0-9]+\\.?[0-9]*")[[1]]
    z <- suppressWarnings(as.numeric(z))
    if (!length(z)) NA_real_ else mean(z, na.rm = TRUE)
  }, numeric(1))
}

recode_treatment <- function(x) {
  lookup <- c(
    B = "B (vs matched control)",
    BM = "B+M (manure)",
    BN = "B+N",
    BNK = "B+N+K",
    BNP = "B+N+P",
    BNPK = "B+NPK",
    BNPKM = "B+NPK+M",
    BP = "B+P",
    BPK = "B+P+K"
  )
  
  y <- clean_chr(x)
  out <- y
  matched <- !is.na(y) & y %in% names(lookup)
  out[matched] <- unname(lookup[y[matched]])
  out
}

col_or_na <- function(data, nm) {
  if (nm %in% names(data)) safe_num(data[[nm]]) else
    rep(NA_real_, nrow(data))
}

total_tau2 <- function(fit) sum(fit$sigma2, na.rm = TRUE)

pseudo_r2 <- function(null_fit, mod_fit) {
  t0 <- total_tau2(null_fit)
  t1 <- total_tau2(mod_fit)
  if (!is.finite(t0) || t0 <= 0) return(NA_real_)
  100 * (t0 - t1) / t0
}

safe_cr2 <- function(fit, cluster) {
  tryCatch(
    clubSandwich::coef_test(
      fit,
      vcov = "CR2",
      cluster = cluster,
      test = "Satterthwaite"
    ),
    error = function(e) structure(
      list(error = conditionMessage(e)),
      class = "cr2_error"
    )
  )
}

cr2_tidy <- function(x, model_name) {
  if (inherits(x, "cr2_error")) {
    return(tibble(
      model = model_name, term = NA_character_,
      estimate = NA_real_, se_cr2 = NA_real_,
      df_satt = NA_real_, p_cr2 = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_,
      inference = "CR2 failed", note = x$error
    ))
  }
  
  z <- as.data.frame(x)
  crit <- stats::qt(
    1 - alpha / 2,
    df = z$df_Satt
  )
  
  tibble(
    model = model_name,
    term = rownames(z),
    estimate = z$beta,
    se_cr2 = z$SE,
    df_satt = z$df_Satt,
    p_cr2 = z$p_Satt,
    ci_lb = z$beta - crit * z$SE,
    ci_ub = z$beta + crit * z$SE,
    inference = "CR2-Satterthwaite",
    low_df_warning = is.finite(z$df_Satt) &
      z$df_Satt < low_cr2_df_threshold,
    note = ifelse(
      is.finite(z$df_Satt) &
        z$df_Satt < low_cr2_df_threshold,
      "Very low CR2 degrees of freedom; interpret cautiously",
      NA_character_
    )
  )
}

safe_wald <- function(fit, constraints, cluster) {
  tryCatch(
    clubSandwich::Wald_test(
      fit,
      constraints = constraints,
      vcov = "CR2",
      cluster = cluster,
      test = "HTZ"
    ),
    error = function(e) structure(
      list(error = conditionMessage(e)),
      class = "wald_error"
    )
  )
}

wald_tidy <- function(x, moderator) {
  if (inherits(x, "wald_error")) {
    return(tibble(
      moderator = moderator, test = "CR2 HTZ omnibus",
      Fstat = NA_real_, df_num = NA_real_,
      df_denom = NA_real_, p_omnibus = NA_real_,
      note = x$error
    ))
  }
  
  z <- as.data.frame(x)
  
  tibble(
    moderator = moderator,
    test = "CR2 HTZ omnibus",
    Fstat = z$Fstat[1],
    df_num = z$df_num[1],
    df_denom = z$df_denom[1],
    p_omnibus = z$p_val[1],
    note = ifelse(
      z$df_denom[1] < low_cr2_df_threshold,
      "Very low denominator df; interpret cautiously",
      NA_character_
    )
  )
}

extract_heterogeneity <- function(fit, data, model_name) {
  s <- fit$sigma2
  tau_article <- if (length(s) >= 1) s[1] else NA_real_
  tau_effect <- if (length(s) >= 2) s[2] else NA_real_
  tau_total <- sum(s, na.rm = TRUE)
  mean_vi <- mean(data$vi, na.rm = TRUE)
  
  tibble(
    model = model_name,
    k = fit$k,
    m = dplyr::n_distinct(data$article_id),
    QE = as.numeric(fit$QE),
    QE_df = if (!is.null(fit$QEdf)) {
      as.numeric(fit$QEdf)
    } else {
      fit$k - fit$p
    },
    QE_p = as.numeric(fit$QEp),
    tau2_between_article = tau_article,
    tau2_within_article = tau_effect,
    tau2_total = tau_total,
    mean_sampling_variance = mean_vi,
    I2_approx_pct = if (
      is.finite(tau_total + mean_vi) && tau_total + mean_vi > 0
    ) 100 * tau_total / (tau_total + mean_vi) else NA_real_
  )
}

compare_models_ml <- function(data, moderator_formula, moderator) {
  null_fit <- tryCatch(
    metafor::rma.mv(
      yi, vi, random = ~ 1 | article_id/effect_id,
      method = "ML", data = data
    ),
    error = function(e) NULL
  )
  
  full_fit <- tryCatch(
    metafor::rma.mv(
      yi, vi, mods = moderator_formula,
      random = ~ 1 | article_id/effect_id,
      method = "ML", data = data
    ),
    error = function(e) NULL
  )
  
  if (is.null(null_fit) || is.null(full_fit)) {
    return(tibble(
      moderator = moderator, k = nrow(data),
      m = n_distinct(data$article_id),
      AIC_null = NA_real_, AIC_full = NA_real_,
      delta_AIC = NA_real_, BIC_null = NA_real_,
      BIC_full = NA_real_, delta_BIC = NA_real_,
      LRT = NA_real_, df_LRT = NA_real_, p_LRT = NA_real_,
      status = "ML model fitting failed"
    ))
  }
  
  ll0 <- logLik(null_fit)
  ll1 <- logLik(full_fit)
  ddf <- attr(ll1, "df") - attr(ll0, "df")
  lrt0 <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  valid <- is.finite(lrt0) && is.finite(ddf) &&
    ddf > 0 && lrt0 >= -1e-7
  lrt <- if (valid) max(0, lrt0) else NA_real_
  
  tibble(
    moderator = moderator,
    k = nrow(data),
    m = n_distinct(data$article_id),
    AIC_null = AIC(null_fit),
    AIC_full = AIC(full_fit),
    delta_AIC = AIC(full_fit) - AIC(null_fit),
    BIC_null = BIC(null_fit),
    BIC_full = BIC(full_fit),
    delta_BIC = BIC(full_fit) - BIC(null_fit),
    LRT = lrt,
    df_LRT = if (valid) ddf else NA_real_,
    p_LRT = if (valid) pchisq(lrt, ddf, lower.tail = FALSE) else NA_real_,
    status = if (valid) {
      "ML comparison of nested fixed-effect models"
    } else {
      "AIC/BIC available; LRT unavailable or invalid"
    }
  )
}

category_thresholds <- function(var) {
  if (var == "texture") return(list(k = 3L, m = 2L))
  if (var == "country") return(list(k = 4L, m = 2L))
  list(k = min_effects_per_category, m = min_articles_per_category)
}

make_working_V <- function(vi, cluster, rho) {
  V <- diag(vi)
  for (ii in split(seq_along(vi), cluster)) {
    if (length(ii) > 1L) {
      block <- rho * outer(sqrt(vi[ii]), sqrt(vi[ii]))
      diag(block) <- vi[ii]
      V[ii, ii] <- block
    }
  }
  V
}

write_csv_safe <- function(x, path) {
  readr::write_csv(x, path, na = "")
}

# -------------------- 4. Moderator specifications ------------

continuous_spec <- tibble(
  var = c(
    "pyrolysis_temp", "biochar_rate_log10", "c_biochar",
    "n_biochar", "ph_biochar", "trial_duration_log10",
    "soil_depth", "soc_initial", "ph_initial",
    "temperature", "precipitation"
  ),
  label = c(
    "Pyrolysis temperature", "Biochar application rate",
    "Biochar C content", "Biochar N content", "Biochar pH",
    "Trial duration", "Sampling soil depth", "Initial SOC",
    "Initial soil pH", "Mean annual temperature",
    "Mean annual precipitation"
  ),
  unit = c(
    "degC", "t/ha", "g/kg", "g/kg", "", "days",
    "cm", "g/kg", "", "degC", "mm"
  ),
  transform = c(
    "none", "log10", "none", "none", "none", "log10",
    "none", "none", "none", "none", "none"
  )
)

categorical_spec <- tibble(
  var = c(
    "feedstock", "treatment_type", "crop",
    "climate_zone", "country", "texture"
  ),
  source = c(
    "Type_biochar", "Treatment", "CropType",
    "ClimateZone", "Country", "Texture_initial"
  )
)

supplement_continuous <- setdiff(
  continuous_spec$var, main_text_continuous
)
supplement_categorical <- setdiff(
  categorical_spec$var, main_text_categorical
)
supplement_moderators <- c(
  supplement_continuous, supplement_categorical
)

# -------------------- 5. Load and validate data --------------

if (!file.exists(input_xlsx)) {
  stop("Input workbook not found: ", input_xlsx)
}

dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_root, "templates"),
           showWarnings = FALSE, recursive = TRUE)

raw <- readxl::read_excel(
  input_xlsx,
  sheet = input_sheet,
  .name_repair = "minimal"
)

names(raw) <- stringr::str_remove_all(names(raw), "\\s+")

required <- c(
  "Study_ID", "Country", "ClimateZone", "SOC_initial",
  "pH_initial", "Texture_initial", "Type_biochar",
  "PyrolysisTemperature_biochar", "C_biochar", "N_biochar",
  "pH_biochar", "Treatment", "BiocharAddition", "CropType",
  "TrialDuration", "SoilDepth"
)

missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

outcome_cols <- unlist(lapply(
  outcomes_requested,
  function(o) paste0(
    o, c("_T", "_T_SD", "_T_N", "_CK", "_CK_SD", "_CK_N")
  )
))

missing_outcome_cols <- setdiff(outcome_cols, names(raw))
if (length(missing_outcome_cols)) {
  stop(
    "Missing required outcome columns: ",
    paste(missing_outcome_cols, collapse = ", ")
  )
}

raw <- raw |>
  mutate(Study_ID = clean_chr(Study_ID))

if (anyNA(raw$Study_ID)) {
  stop(
    "Missing Study_ID in Excel data row(s): ",
    paste(which(is.na(raw$Study_ID)) + 1L, collapse = ", ")
  )
}

base <- raw |>
  mutate(sheet_row = row_number() + 1L) |>
  transmute(
    sheet_row,
    Study_ID,
    article_id = Study_ID,
    comparison_id = paste0(Study_ID, " [row ", sheet_row, "]"),
    country = clean_chr(Country),
    climate_zone = clean_chr(ClimateZone),
    texture = clean_texture(Texture_initial),
    feedstock_raw = clean_feedstock(Type_biochar),
    feedstock = factor(
      clean_feedstock(Type_biochar),
      levels = expected_feedstocks
    ),
    treatment_type = recode_treatment(Treatment),
    crop = clean_chr(CropType),
    pyrolysis_temp = parse_num_mid(PyrolysisTemperature_biochar),
    biochar_rate_log10 = ifelse(
      safe_num(BiocharAddition) > 0,
      log10(safe_num(BiocharAddition)), NA_real_
    ),
    c_biochar = safe_num(C_biochar),
    n_biochar = safe_num(N_biochar),
    ph_biochar = safe_num(pH_biochar),
    trial_duration_log10 = ifelse(
      safe_num(TrialDuration) > 0,
      log10(safe_num(TrialDuration)), NA_real_
    ),
    soil_depth = safe_num(SoilDepth),
    soc_initial = safe_num(SOC_initial),
    ph_initial = safe_num(pH_initial),
    temperature = col_or_na(raw, "Temperature"),
    precipitation = col_or_na(raw, "Precipitation")
  )

# Audits
write_csv_safe(
  base |> count(Study_ID, name = "n_comparisons") |>
    arrange(desc(n_comparisons), Study_ID),
  file.path(output_root, "study_id_cluster_audit.csv")
)

write_csv_safe(
  base |> select(
    sheet_row, Study_ID, comparison_id, feedstock_raw, feedstock
  ),
  file.path(output_root, "feedstock_row_audit.csv")
)

write_csv_safe(
  base |> count(feedstock, .drop = FALSE, name = "n_rows"),
  file.path(output_root, "feedstock_value_audit.csv")
)

write_csv_safe(
  bind_rows(
    continuous_spec |> transmute(
      requested_moderator = label,
      analysis_variable = var, type = "continuous"
    ),
    categorical_spec |> transmute(
      requested_moderator = source,
      analysis_variable = var, type = "categorical"
    )
  ) |>
    mutate(
      n_nonmissing_all_rows = vapply(
        analysis_variable,
        function(v) sum(!is.na(base[[v]])),
        numeric(1)
      ),
      analysis_priority = ifelse(
        analysis_variable %in% main_text_moderators,
        "Main text", "Supplement"
      )
    ),
  file.path(output_root, "moderator_availability_audit.csv")
)

write_csv_safe(
  tibble(
    item = c(
      "Study_ID/article_id", "comparison_id",
      "Pyrolysis temperature", "Biochar application rate",
      "Trial duration", "Feedstock", "SOC_initial"
    ),
    treatment = c(
      "Study_ID identifies the independent publication cluster",
      "Study_ID plus Excel row number identifies a comparison",
      "Range values parsed as midpoints",
      "Log10 transformed before standardisation",
      "Log10 transformed before standardisation",
      "Original Wheat straw and Maize straw categories retained",
      "Retained as a continuous moderator; not analysed as an outcome"
    )
  ),
  file.path(output_root, "data_coding_audit.csv")
)

write_csv_safe(
  tibble(
    Study_ID = sort(unique(base$Study_ID)),
    reviewer_1 = NA_character_,
    reviewer_2 = NA_character_,
    consensus = NA_character_,
    sampling_design = NA_character_,
    comparability_of_control = NA_character_,
    outcome_measurement = NA_character_,
    completeness_of_reporting = NA_character_,
    selective_reporting = NA_character_,
    overall_judgment = NA_character_,
    rationale = NA_character_
  ),
  file.path(output_root, "templates", "risk_of_bias_appraisal_template.csv")
)

# -------------------- 6. Outcome analysis --------------------

run_outcome <- function(outcome) {
  message("\nAnalysing outcome: ", outcome)
  
  out <- file.path(output_root, outcome)
  fig_dir <- file.path(out, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "diagnostics"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "main_text"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "supplement"),
             recursive = TRUE, showWarnings = FALSE)
  
  cols <- paste0(
    outcome, c("_T", "_T_SD", "_T_N", "_CK", "_CK_SD", "_CK_N")
  )
  
  d0 <- base |>
    mutate(
      m1i = safe_num(raw[[cols[1]]]),
      sd1i = safe_num(raw[[cols[2]]]),
      n1i = safe_num(raw[[cols[3]]]),
      m2i = safe_num(raw[[cols[4]]]),
      sd2i = safe_num(raw[[cols[5]]]),
      n2i = safe_num(raw[[cols[6]]])
    )
  
  reported <- d0 |>
    filter(if_any(
      c(m1i, sd1i, n1i, m2i, sd2i, n2i), ~ !is.na(.x)
    ))
  
  excluded <- reported |>
    mutate(
      reason = case_when(
        !is.finite(m1i) | !is.finite(m2i) ~ "Non-finite mean",
        !is.finite(sd1i) | !is.finite(sd2i) ~ "Non-finite SD",
        !is.finite(n1i) | !is.finite(n2i) ~ "Non-finite sample size",
        m1i <= 0 | m2i <= 0 ~ "Non-positive mean; lnRR undefined",
        sd1i < 0 | sd2i < 0 ~ "Negative SD",
        n1i < 2 | n2i < 2 ~ "Sample size below 2",
        sd1i == 0 | sd2i == 0 ~
          "Zero SD; review before calculating sampling variance",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(reason))
  
  write_csv_safe(excluded, file.path(out, "excluded_rows.csv"))
  
  # Strictly exclude zero SDs: escalc ROM sampling variance is not
  # safely handled by silently retaining zero-SD comparisons.
  dat <- reported |>
    filter(
      is.finite(m1i), is.finite(m2i),
      is.finite(sd1i), is.finite(sd2i),
      is.finite(n1i), is.finite(n2i),
      m1i > 0, m2i > 0,
      sd1i > 0, sd2i > 0,
      n1i >= 2, n2i >= 2
    ) |>
    mutate(
      effect_id = row_number(),
      article_id = factor(article_id),
      country = factor(country),
      climate_zone = factor(climate_zone),
      texture = factor(texture),
      feedstock = factor(
        feedstock, levels = expected_feedstocks
      ),
      treatment_type = factor(treatment_type),
      crop = factor(crop)
    )
  
  if (nrow(dat) < min_effects_moderator ||
      n_distinct(dat$article_id) < 2L) {
    stop(
      outcome, ": insufficient usable comparisons or independent articles."
    )
  }
  
  es <- metafor::escalc(
    measure = "ROM",
    m1i = m1i, sd1i = sd1i, n1i = n1i,
    m2i = m2i, sd2i = sd2i, n2i = n2i,
    data = dat
  ) |>
    as_tibble() |>
    mutate(
      sei = sqrt(vi),
      yi_pct = lnrr_to_pct(yi)
    )
  
  if (any(!is.finite(es$vi)) || any(es$vi <= 0)) {
    stop(outcome, ": invalid or non-positive sampling variances remain.")
  }
  
  counts <- tibble(
    outcome = outcome,
    rows_in_sheet = nrow(d0),
    rows_reporting_outcome = nrow(reported),
    rows_excluded_unusable = nrow(excluded),
    k_comparisons = nrow(es),
    m_articles = n_distinct(es$article_id),
    countries = n_distinct(es$country, na.rm = TRUE),
    climate_zones = n_distinct(es$climate_zone, na.rm = TRUE)
  )
  
  write_csv_safe(counts, file.path(out, "evidence_base_counts.csv"))
  write_csv_safe(es, file.path(out, "effect_sizes_lnRR.csv"))
  
  write_csv_safe(
    es |> count(article_id, name = "k_effects"),
    file.path(out, "article_cluster_audit.csv")
  )
  
  # Primary intercept-only multilevel REML model.
  main_model <- metafor::rma.mv(
    yi, vi,
    random = ~ 1 | article_id/effect_id,
    method = "REML",
    data = es
  )
  
  main_cr2 <- safe_cr2(main_model, es$article_id)
  main_cr2_tab <- cr2_tidy(main_cr2, "Overall multilevel model") |>
    mutate(
      estimate_pct = lnrr_to_pct(estimate),
      pct_lb = lnrr_to_pct(ci_lb),
      pct_ub = lnrr_to_pct(ci_ub)
    )
  
  write_csv_safe(
    main_cr2_tab,
    file.path(out, "overall_effect_CR2.csv")
  )
  
  overall_heterogeneity <- extract_heterogeneity(
    main_model, es, "Overall intercept-only model"
  )
  
  write_csv_safe(
    overall_heterogeneity,
    file.path(out, "overall_heterogeneity_assessment.csv")
  )
  
  main_summary <- tibble(
    k = main_model$k,
    m = n_distinct(es$article_id),
    lnRR = as.numeric(main_model$b[1]),
    model_se = main_model$se[1],
    model_ci_lb = main_model$ci.lb[1],
    model_ci_ub = main_model$ci.ub[1],
    model_p = main_model$pval[1],
    pct = lnrr_to_pct(main_model$b[1]),
    pct_lb = lnrr_to_pct(main_model$ci.lb[1]),
    pct_ub = lnrr_to_pct(main_model$ci.ub[1]),
    tau2_total = total_tau2(main_model)
  )
  
  write_csv_safe(
    main_summary,
    file.path(out, "overall_model_and_heterogeneity.csv")
  )
  
  pred <- predict(main_model)
  write_csv_safe(
    tibble(
      estimate_lnRR = pred$pred,
      ci_lb_lnRR = pred$ci.lb,
      ci_ub_lnRR = pred$ci.ub,
      pi_lb_lnRR = pred$pi.lb,
      pi_ub_lnRR = pred$pi.ub,
      estimate_pct = lnrr_to_pct(pred$pred),
      ci_lb_pct = lnrr_to_pct(pred$ci.lb),
      ci_ub_pct = lnrr_to_pct(pred$ci.ub),
      pi_lb_pct = lnrr_to_pct(pred$pi.lb),
      pi_ub_pct = lnrr_to_pct(pred$pi.ub)
    ),
    file.path(out, "prediction_interval.csv")
  )
  
  # ---------------- Continuous moderators ------------------
  
  cont_terms <- list()
  cont_tests <- list()
  cont_models <- list()
  cont_support <- list()
  cont_comparisons <- list()
  cont_heterogeneity <- list()
  
  for (v in continuous_spec$var) {
    d <- es |>
      filter(is.finite(.data[[v]])) |>
      droplevels()
    
    k <- nrow(d)
    m <- n_distinct(d$article_id)
    u <- n_distinct(d[[v]])
    
    eligible <- k >= min_effects_moderator &&
      m >= min_articles_moderator &&
      u >= min_unique_continuous_values &&
      is.finite(sd(d[[v]])) && sd(d[[v]]) > 0
    
    cont_support[[v]] <- tibble(
      moderator = v, k = k, m = m, unique_values = u,
      support_status = ifelse(eligible, "Eligible", "Not eligible")
    )
    
    if (!eligible) {
      cont_tests[[v]] <- tibble(
        moderator = v, test = "Not fitted",
        Fstat = NA_real_, df_num = NA_real_,
        df_denom = NA_real_, p_omnibus = NA_real_,
        pseudo_R2_pct = NA_real_, k = k, m = m,
        note = "Insufficient support or moderator variation"
      )
      next
    }
    
    d <- d |> mutate(mod_z = as.numeric(scale(.data[[v]])))
    
    null_fit <- metafor::rma.mv(
      yi, vi, random = ~ 1 | article_id/effect_id,
      method = "REML", data = d
    )
    
    fit <- metafor::rma.mv(
      yi, vi, mods = ~ mod_z,
      random = ~ 1 | article_id/effect_id,
      method = "REML", data = d
    )
    
    cont_models[[v]] <- fit
    
    cont_comparisons[[v]] <- compare_models_ml(
      d, ~ mod_z, v
    )
    
    cont_heterogeneity[[v]] <- extract_heterogeneity(
      fit, d, paste0("Continuous: ", v)
    )
    
    cont_terms[[v]] <- cr2_tidy(
      safe_cr2(fit, d$article_id),
      paste0("Continuous: ", v)
    ) |>
      mutate(
        moderator = v, k = k, m = m,
        scale_note = paste0(
          "Coefficient per SD; mean=", round(mean(d[[v]]), 3),
          "; SD=", round(sd(d[[v]]), 3)
        )
      )
    
    cont_tests[[v]] <- wald_tidy(
      safe_wald(fit, clubSandwich::constrain_zero(2), d$article_id),
      v
    ) |>
      mutate(
        pseudo_R2_pct = pseudo_r2(null_fit, fit),
        k = k, m = m
      )
  }
  
  cont_terms_df <- bind_rows(cont_terms)
  cont_tests_df <- bind_rows(cont_tests)
  cont_support_df <- bind_rows(cont_support)
  cont_comparisons_df <- bind_rows(cont_comparisons)
  cont_heterogeneity_df <- bind_rows(cont_heterogeneity)
  
  # ---------------- Categorical moderators -----------------
  
  cat_levels <- list()
  cat_tests <- list()
  cat_models <- list()
  cat_support <- list()
  cat_level_support <- list()
  cat_comparisons <- list()
  cat_heterogeneity <- list()
  
  for (v in categorical_spec$var) {
    d0v <- es |>
      filter(!is.na(.data[[v]])) |>
      droplevels()
    
    thr <- category_thresholds(v)
    
    level_support <- d0v |>
      mutate(level = as.character(.data[[v]])) |>
      group_by(level) |>
      summarise(
        k = n(),
        m = n_distinct(article_id),
        .groups = "drop"
      ) |>
      mutate(
        moderator = v,
        eligible_level = k >= thr$k & m >= thr$m
      )
    
    cat_level_support[[v]] <- level_support
    
    keep <- level_support |>
      filter(eligible_level) |>
      pull(level)
    
    d <- d0v |>
      filter(as.character(.data[[v]]) %in% keep) |>
      droplevels()
    
    lvls <- levels(d[[v]])
    k <- nrow(d)
    m <- n_distinct(d$article_id)
    eligible <- length(lvls) >= 2L &&
      k >= min_effects_moderator &&
      m >= min_articles_moderator
    
    cat_support[[v]] <- tibble(
      moderator = v, eligible_levels = length(lvls),
      k = k, m = m,
      support_status = ifelse(eligible, "Eligible", "Not eligible")
    )
    
    if (!eligible) {
      cat_tests[[v]] <- tibble(
        moderator = v, test = "Not fitted",
        Fstat = NA_real_, df_num = NA_real_,
        df_denom = NA_real_, p_omnibus = NA_real_,
        pseudo_R2_pct = NA_real_, k = k, m = m,
        note = "Insufficient supported levels or article support"
      )
      next
    }
    
    # No-intercept parameterisation estimates one mean per level.
    form <- stats::as.formula(paste0("~ 0 + ", v))
    
    null_fit <- metafor::rma.mv(
      yi, vi, random = ~ 1 | article_id/effect_id,
      method = "REML", data = d
    )
    
    fit <- metafor::rma.mv(
      yi, vi, mods = form,
      random = ~ 1 | article_id/effect_id,
      method = "REML", data = d
    )
    
    cat_models[[v]] <- fit
    cat_comparisons[[v]] <- compare_models_ml(d, form, v)
    cat_heterogeneity[[v]] <- extract_heterogeneity(
      fit, d, paste0("Categorical: ", v)
    )
    
    # Robust coefficient-to-level mapping:
    # use the model matrix column names, not the order of factor levels.
    mm <- model.matrix(form, data = d)
    coefficient_names <- colnames(mm)
    coefficient_names <- make.names(coefficient_names, unique = FALSE)
    
    cr <- safe_cr2(fit, d$article_id)
    crtab <- cr2_tidy(cr, paste0("Categorical: ", v))
    
    if (!inherits(cr, "cr2_error")) {
      z <- as.data.frame(cr)
      cr_names <- rownames(z)
      
      level_map <- tibble(
        term = cr_names,
        coefficient_index = seq_along(cr_names)
      )
      
      # For ~ 0 + factor, model-matrix columns correspond to factor
      # levels. Match the actual fitted coefficient names.
      mm_names <- colnames(model.matrix(form, data = d))
      level_map$level <- NA_character_
      
      for (j in seq_along(mm_names)) {
        candidate <- mm_names[j]
        matched <- which(
          cr_names == candidate |
            make.names(cr_names) == make.names(candidate)
        )
        if (length(matched) == 1L) {
          level_map$level[matched] <- lvls[j]
        }
      }
      
      # A failed exact mapping is reported, never silently position-matched.
      if (anyNA(level_map$level) || anyDuplicated(level_map$level)) {
        stop(
          "Could not unambiguously map CR2 coefficients to levels for ",
          v, ". Coefficients: ",
          paste(cr_names, collapse = ", "),
          "; model matrix: ",
          paste(mm_names, collapse = ", ")
        )
      }
      
      ct_levels <- crtab |>
        left_join(level_map |> select(term, level), by = "term")
    } else {
      ct_levels <- tibble(
        model = paste0("Categorical: ", v),
        term = NA_character_,
        estimate = NA_real_,
        se_cr2 = NA_real_,
        df_satt = NA_real_,
        p_cr2 = NA_real_,
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        inference = "CR2 failed",
        level = lvls,
        note = cr$error
      )
    }
    
    level_counts <- d |>
      mutate(level = as.character(.data[[v]])) |>
      group_by(level) |>
      summarise(
        k = n(), m = n_distinct(article_id),
        .groups = "drop"
      )
    
    cat_levels[[v]] <- ct_levels |>
      left_join(level_counts, by = "level") |>
      mutate(
        moderator = v,
        estimate_pct = lnrr_to_pct(estimate),
        pct_lb = lnrr_to_pct(ci_lb),
        pct_ub = lnrr_to_pct(ci_ub)
      )
    
    # Omnibus equality of all level means.
    cat_tests[[v]] <- wald_tidy(
      safe_wald(
        fit,
        clubSandwich::constrain_equal(seq_along(lvls)),
        d$article_id
      ),
      v
    ) |>
      mutate(
        pseudo_R2_pct = pseudo_r2(null_fit, fit),
        k = k, m = m
      )
  }
  
  cat_levels_df <- bind_rows(cat_levels)
  cat_tests_df <- bind_rows(cat_tests)
  cat_support_df <- bind_rows(cat_support)
  cat_level_support_df <- bind_rows(cat_level_support)
  cat_comparisons_df <- bind_rows(cat_comparisons)
  cat_heterogeneity_df <- bind_rows(cat_heterogeneity)
  
  # ---------------- Combined moderator tables ---------------
  
  moderator_omnibus <- bind_rows(
    cont_tests_df |> mutate(model_type = "Continuous"),
    cat_tests_df |> mutate(model_type = "Categorical")
  ) |>
    mutate(
      analysis_priority = ifelse(
        moderator %in% main_text_moderators,
        "Main text", "Supplement"
      )
    ) |>
    group_by(analysis_priority) |>
    mutate(p_holm = p.adjust(p_omnibus, method = "holm")) |>
    ungroup()
  
  model_comparisons <- bind_rows(
    cont_comparisons_df |> mutate(model_type = "Continuous"),
    cat_comparisons_df |> mutate(model_type = "Categorical")
  ) |>
    mutate(
      analysis_priority = ifelse(
        moderator %in% main_text_moderators,
        "Main text", "Supplement"
      )
    ) |>
    group_by(analysis_priority) |>
    mutate(p_LRT_holm = p.adjust(p_LRT, method = "holm")) |>
    ungroup()
  
  moderator_heterogeneity <- bind_rows(
    cont_heterogeneity_df,
    cat_heterogeneity_df
  ) |>
    mutate(
      moderator = stringr::str_remove(
        model, "^(Continuous|Categorical): "
      ),
      analysis_priority = ifelse(
        moderator %in% main_text_moderators,
        "Main text", "Supplement"
      )
    )
  
  # Save full tables.
  tables <- list(
    continuous_moderator_support = cont_support_df,
    continuous_moderator_coefficients_CR2 = cont_terms_df,
    continuous_moderator_omnibus_R2 = cont_tests_df,
    continuous_model_comparisons_ML = cont_comparisons_df,
    continuous_moderator_heterogeneity = cont_heterogeneity_df,
    categorical_moderator_support = cat_support_df,
    categorical_moderator_level_support = cat_level_support_df,
    categorical_subgroup_complete_results_CR2 = cat_levels_df,
    categorical_moderator_omnibus_R2 = cat_tests_df,
    categorical_model_comparisons_ML = cat_comparisons_df,
    categorical_moderator_heterogeneity = cat_heterogeneity_df,
    all_moderators_omnibus = moderator_omnibus,
    model_comparisons_ML_all = model_comparisons,
    moderator_heterogeneity_all = moderator_heterogeneity
  )
  
  purrr::iwalk(tables, function(x, nm) {
    write_csv_safe(x, file.path(out, paste0(nm, ".csv")))
    write_csv_safe(
      filter(x, if ("moderator" %in% names(x)) {
        moderator %in% main_text_moderators
      } else {
        rep(TRUE, nrow(x))
      }),
      file.path(out, "main_text", paste0(nm, ".csv"))
    )
    write_csv_safe(
      filter(x, if ("moderator" %in% names(x)) {
        moderator %in% supplement_moderators
      } else {
        rep(FALSE, nrow(x))
      }),
      file.path(out, "supplement", paste0(nm, ".csv"))
    )
  })
  
  # ---------------- Leave-one-article-out --------------------
  
  loo <- purrr::map_dfr(
    levels(es$article_id),
    function(s) {
      d <- es |> filter(article_id != s) |> droplevels()
      if (n_distinct(d$article_id) < 3L) return(tibble())
      
      fit <- tryCatch(
        metafor::rma.mv(
          yi, vi, random = ~ 1 | article_id/effect_id,
          method = "REML", data = d
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(tibble())
      
      cr <- cr2_tidy(
        safe_cr2(fit, d$article_id),
        paste0("Omit ", s)
      ) |> slice(1)
      
      tibble(
        omitted_article = s,
        k = nrow(d),
        m = n_distinct(d$article_id),
        lnRR = as.numeric(fit$b[1]),
        pct = lnrr_to_pct(fit$b[1]),
        model_ci_lb = fit$ci.lb[1],
        model_ci_ub = fit$ci.ub[1],
        pct_lb = lnrr_to_pct(cr$ci_lb[1]),
        pct_ub = lnrr_to_pct(cr$ci_ub[1]),
        cr2_se = cr$se_cr2[1],
        cr2_df = cr$df_satt[1],
        cr2_p = cr$p_cr2[1]
      )
    }
  )
  
  write_csv_safe(loo, file.path(out, "leave_one_article_out.csv"))
  
  # ---------------- Assumed dependence sensitivity ---------
  
  rho_values <- c(0, 0.3, 0.6, 0.9)
  rho_models <- list()
  
  rho_results <- purrr::map_dfr(rho_values, function(rho) {
    V <- make_working_V(es$vi, es$article_id, rho)
    
    fit <- metafor::rma.mv(
      yi, V = V,
      random = ~ 1 | article_id/effect_id,
      method = "REML", data = es
    )
    
    rho_models[[as.character(rho)]] <<- fit
    
    cr <- cr2_tidy(
      safe_cr2(fit, es$article_id),
      paste0("rho=", rho)
    ) |> slice(1)
    
    tibble(
      assumed_rho = rho,
      lnRR = as.numeric(fit$b[1]),
      pct = lnrr_to_pct(fit$b[1]),
      model_ci_lb = fit$ci.lb[1],
      model_ci_ub = fit$ci.ub[1],
      pct_lb = lnrr_to_pct(cr$ci_lb[1]),
      pct_ub = lnrr_to_pct(cr$ci_ub[1]),
      cr2_se = cr$se_cr2[1],
      cr2_df = cr$df_satt[1],
      cr2_p = cr$p_cr2[1]
    )
  })
  
  write_csv_safe(
    rho_results,
    file.path(out, "within_article_correlation_sensitivity.csv")
  )
  
  writeLines(
    c(
      "Interpretation note",
      "The rho values are assumed working within-article sampling correlations.",
      "They are not estimates of the true within-article correlation.",
      "Results assess sensitivity to alternative dependence assumptions."
    ),
    file.path(out, "diagnostics", "dependence_sensitivity_note.txt")
  )
  
  # ---------------- Small-study effects --------------------
  
  article_es <- es |>
    mutate(w = 1 / vi) |>
    group_by(article_id) |>
    summarise(
      k_within = n(),
      yi = sum(w * yi) / sum(w),
      vi = 1 / sum(w),
      sei = sqrt(vi),
      .groups = "drop"
    )
  
  write_csv_safe(
    article_es,
    file.path(out, "article_level_effects_for_funnel.csv")
  )
  
  egger_fit <- NULL
  egger_result <- tibble(
    outcome = outcome,
    k = nrow(es),
    m = n_distinct(es$article_id),
    status = "Not tested: insufficient independent article support",
    term = NA_character_,
    estimate = NA_real_,
    se_cr2 = NA_real_,
    df_satt = NA_real_,
    p_cr2 = NA_real_,
    ci_lb = NA_real_,
    ci_ub = NA_real_
  )
  
  if (nrow(es) >= min_effects_egger &&
      n_distinct(es$article_id) >= min_articles_egger) {
    egger_fit <- tryCatch(
      metafor::rma.mv(
        yi, vi, mods = ~ sei,
        random = ~ 1 | article_id/effect_id,
        method = "REML", data = es
      ),
      error = function(e) NULL
    )
    
    if (!is.null(egger_fit)) {
      egger_result <- cr2_tidy(
        safe_cr2(egger_fit, es$article_id),
        "Multilevel Egger-type regression"
      ) |>
        filter(term == "sei") |>
        mutate(
          outcome = outcome,
          k = nrow(es),
          m = n_distinct(es$article_id),
          status = "Fitted; asymmetry is not proof of publication bias"
        )
    }
  }
  
  write_csv_safe(
    egger_result,
    file.path(out, "small_study_effects_multilevel_CR2.csv")
  )
  
  article_uni <- tryCatch(
    metafor::rma.uni(yi, vi, method = "REML", data = article_es),
    error = function(e) NULL
  )
  
  rank_res <- if (!is.null(article_uni)) {
    tryCatch(metafor::ranktest(article_uni), error = function(e) NULL)
  } else NULL
  
  write_csv_safe(
    tibble(
      diagnostic = "Exploratory article-level Begg-Mazumdar test",
      m_articles = nrow(article_es),
      statistic = if (is.null(rank_res)) NA_real_ else unname(rank_res$tau),
      p_value = if (is.null(rank_res)) NA_real_ else rank_res$pval,
      note = "Exploratory only; low power with few independent articles"
    ),
    file.path(out, "rank_correlation_exploratory.csv")
  )
  
  tf <- if (!is.null(article_uni) &&
            nrow(article_es) >= min_articles_trimfill) {
    tryCatch(metafor::trimfill(article_uni), error = function(e) NULL)
  } else NULL
  
  write_csv_safe(
    tibble(
      status = if (is.null(tf)) {
        "Not performed or model failed; exploratory diagnostic only"
      } else {
        "Exploratory article-level trim-and-fill"
      },
      m_articles = nrow(article_es),
      imputed_articles = if (is.null(tf)) NA_real_ else tf$k0,
      adjusted_lnRR = if (is.null(tf)) NA_real_ else as.numeric(tf$b[1]),
      adjusted_ci_lb = if (is.null(tf)) NA_real_ else tf$ci.lb[1],
      adjusted_ci_ub = if (is.null(tf)) NA_real_ else tf$ci.ub[1],
      adjusted_pct = if (is.null(tf)) NA_real_ else lnrr_to_pct(tf$b[1])
    ),
    file.path(out, "trim_and_fill_status.csv")
  )
  
  # ---------------- Reproducibility bundle -----------------
  
  bundle <- list(
    outcome = outcome,
    outcome_label = unname(outcome_labels[outcome]),
    data = es,
    counts = counts,
    main_model = main_model,
    main_cr2 = main_cr2_tab,
    main_cr2_raw = main_cr2,
    overall_heterogeneity = overall_heterogeneity,
    continuous_spec = continuous_spec,
    categorical_spec = categorical_spec,
    continuous_models = cont_models,
    categorical_models = cat_models,
    continuous_support = cont_support_df,
    categorical_support = cat_support_df,
    categorical_level_support = cat_level_support_df,
    continuous_terms = cont_terms_df,
    categorical_terms = cat_levels_df,
    continuous_tests = cont_tests_df,
    categorical_tests = cat_tests_df,
    moderator_omnibus = moderator_omnibus,
    model_comparisons = model_comparisons,
    moderator_heterogeneity = moderator_heterogeneity,
    loo = loo,
    rho_results = rho_results,
    rho_models = rho_models,
    article_level_effects = article_es,
    article_aggregated_model = article_uni,
    egger_model = egger_fit,
    main_text_continuous = main_text_continuous,
    main_text_categorical = main_text_categorical,
    supplement_continuous = supplement_continuous,
    supplement_categorical = supplement_categorical,
    session_info = sessionInfo()
  )
  
  saveRDS(bundle, file.path(out, "biochar_meta_models.rds"))
  
  message(
    outcome, ": k=", nrow(es),
    "; m=", n_distinct(es$article_id),
    ". Results saved in ", normalizePath(out)
  )
  
  invisible(bundle)
}

# Run N2O and GWP only.
bundles <- setNames(
  lapply(outcomes_requested, run_outcome),
  outcomes_requested
)

# -------------------- 7. Supplementary tables S1-S8 ---------

# S1 is a transparent template: actual database names, dates,
# search strings and retrieved counts must be filled from records.
S1 <- tibble(
  database = NA_character_,
  platform = NA_character_,
  search_date = NA_character_,
  complete_search_string = NA_character_,
  date_limits = NA_character_,
  language_limits = NA_character_,
  records_retrieved = NA_integer_,
  note = "Template: populate from the actual search records; do not infer."
)

S2 <- base |>
  count(Study_ID, country, climate_zone, crop, name = "n_dataset_rows") |>
  arrange(Study_ID)

S3 <- bind_rows(lapply(bundles, function(b) {
  b$counts
}))

S4 <- bind_rows(lapply(names(bundles), function(o) {
  b <- bundles[[o]]
  cr2 <- b$main_cr2 |>
    dplyr::slice(1) |>
    dplyr::select(
      estimate, se_cr2, df_satt, p_cr2, ci_lb, ci_ub,
      estimate_pct, pct_lb, pct_ub
    )
  het <- b$overall_heterogeneity |>
    dplyr::select(
      k, m, QE, QE_df, QE_p, tau2_between_article,
      tau2_within_article, tau2_total, I2_approx_pct
    )
  tibble::tibble(
    outcome = o,
    model_based_lnRR = as.numeric(b$main_model$b[1]),
    model_based_se = as.numeric(b$main_model$se[1]),
    model_based_ci_lb = as.numeric(b$main_model$ci.lb[1]),
    model_based_ci_ub = as.numeric(b$main_model$ci.ub[1])
  ) |>
    dplyr::bind_cols(cr2, het)
}))

S5 <- bind_rows(lapply(names(bundles), function(o) {
  b <- bundles[[o]]
  bind_rows(
    b$continuous_support,
    b$categorical_support
  ) |>
    mutate(outcome = o)
}))

S6 <- bind_rows(lapply(names(bundles), function(o) {
  b <- bundles[[o]]
  bind_rows(
    b$continuous_terms,
    b$categorical_terms
  ) |>
    mutate(outcome = o)
}))

S7 <- bind_rows(lapply(names(bundles), function(o) {
  b <- bundles[[o]]
  b$model_comparisons |>
    mutate(outcome = o)
}))

S8 <- bind_rows(lapply(names(bundles), function(o) {
  b <- bundles[[o]]
  bind_rows(
    b$loo |> mutate(sensitivity = "Leave-one-article-out"),
    b$rho_results |> mutate(sensitivity = "Assumed within-article correlation")
  ) |>
    mutate(outcome = o)
}))

# Exploratory small-study-effect diagnostics form a second panel of Table S8.
S8_small_study <- bind_rows(lapply(names(bundles), function(o) {
  out_dir <- file.path(output_root, o)
  egger <- readr::read_csv(
    file.path(out_dir, "small_study_effects_multilevel_CR2.csv"),
    show_col_types = FALSE
  ) |> mutate(outcome = o, diagnostic = "Multilevel Egger-type regression")
  begg <- readr::read_csv(
    file.path(out_dir, "rank_correlation_exploratory.csv"),
    show_col_types = FALSE
  ) |> mutate(outcome = o, diagnostic = "Exploratory Begg-Mazumdar test")
  trimfill <- readr::read_csv(
    file.path(out_dir, "trim_and_fill_status.csv"),
    show_col_types = FALSE
  ) |> mutate(outcome = o, diagnostic = "Exploratory article-level trim-and-fill")
  bind_rows(
    tibble::as_tibble(egger),
    tibble::as_tibble(begg),
    tibble::as_tibble(trimfill)
  )
}))

supplementary_tables <- list(
  S1_Search_strategy = S1,
  S2_Study_characteristics = S2,
  S3_Evidence_base = S3,
  S4_Overall_effects_CR2 = S4,
  S5_Moderator_support = S5,
  S6_Moderator_estimates = S6,
  S7_Model_comparisons = S7,
  S8a_Sensitivity_analyses = S8,
  S8b_Small_study_diagnostics = S8_small_study
)

manuscript_crosswalk <- tibble::tribble(
  ~item, ~content, ~methods_location, ~results_location,
  "Figure 1", "Overall N2O and GWP effects", "2.2 Meta-analysis", "3.1 Overall effects",
  "Figure 2", "Feedstock and treatment-type subgroup estimates", "2.2 Moderator analyses", "3.3 Categorical moderators",
  "Figure 3", "Climate-zone and soil-texture subgroup estimates", "2.2 Moderator analyses", "3.3 Categorical moderators",
  "Figure 4", "Continuous moderator meta-regressions", "2.2 Moderator analyses", "3.4 Continuous moderators",
  "Figure S1", "PRISMA flow diagram; supplied from screening records", "2.1 Database compilation", "Study selection",
  "Figure S2", "Leave-one-article-out sensitivity", "2.2 Sensitivity analyses", "3.2 Robustness",
  "Figure S3", "Assumed within-article correlation sensitivity", "2.2 Sensitivity analyses", "3.2 Robustness",
  "Figure S4", "Article-level funnel plots", "2.2 Small-study-effect diagnostics", "3.2 Robustness",
  "Figure S5", "Evidence-base counts", "2.1-2.2 Data extraction and effect sizes", "Opening Results paragraph",
  "Figure S6", "Moderator support diagnostics", "2.2 Moderator support", "3.3-3.4 Moderator analyses",
  "Table S1", "Database search strategy", "2.1 Database compilation", "Study identification",
  "Table S2", "Included-study characteristics", "2.1 Database compilation", "Opening Results paragraph",
  "Table S3", "Effect-size and evidence-base counts", "2.1-2.2 Data extraction and effect sizes", "Opening Results paragraph",
  "Table S4", "Overall CR2 effects and heterogeneity", "2.2 Meta-analysis", "3.1 Overall effects",
  "Table S5", "Moderator support and omnibus tests", "2.2 Moderator analyses", "3.3-3.4 Moderator analyses",
  "Table S6", "Moderator coefficient and subgroup estimates", "2.2 Moderator analyses", "3.3-3.4 Moderator analyses",
  "Table S7", "Nested fixed-effect model comparisons using ML", "2.2 Model comparisons", "3.3-3.4 where reported",
  "Table S8a", "Leave-one-article-out and assumed-rho sensitivity", "2.2 Sensitivity analyses", "3.2 Robustness",
  "Table S8b", "Exploratory Egger, rank-correlation and trim-and-fill diagnostics", "2.2 Small-study-effect diagnostics", "3.2 Robustness"
)

write_csv_safe(
  manuscript_crosswalk,
  file.path(output_root, "Manuscript_Figure_Table_Crosswalk.csv")
)

for (nm in names(supplementary_tables)) {
  write_csv_safe(
    supplementary_tables[[nm]],
    file.path(output_root, paste0(nm, ".csv"))
  )
}

openxlsx::write.xlsx(
  supplementary_tables,
  file = file.path(output_root, "Supplementary_Tables_S1-S8.xlsx"),
  overwrite = TRUE,
  asTable = TRUE
)

# -------------------- 8. PRISMA template --------------------

# These counts are the previously supplied counts, not verified
# against screening records. Confirm them before publication.
prisma <- tibble(
  stage = c(
    "Records identified",
    "Duplicates removed",
    "Records screened",
    "Records excluded",
    "Reports sought for retrieval",
    "Reports not retrieved",
    "Reports assessed for eligibility",
    "Full-text reports excluded",
    "Studies included"
  ),
  n = c(2000, 300, 1700, 1250, 450, 20, 430, 363, 67),
  verification = "User-supplied count; verify against screening log"
)

write_csv_safe(prisma, file.path(output_root, "PRISMA_counts_check.csv"))

# -------------------- 9. Figure helpers ---------------------

theme_pub <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(
        fill = "grey92", color = "grey35"
      ),
      strip.text = element_text(face = "bold")
    )
}

save_plot_all <- function(plot, filename, width, height, dpi = 600) {
  # IMPORTANT: filename may contain a directory. Use basename() here;
  # otherwise file.path(folder, stem) duplicates the directory and causes
  # ggsave() to fail because the resulting nested path does not exist.
  folder <- dirname(filename)
  stem <- tools::file_path_sans_ext(basename(filename))

  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

  output_files <- c(
    png = file.path(folder, paste0(stem, ".png")),
    pdf = file.path(folder, paste0(stem, ".pdf")),
    tiff = file.path(folder, paste0(stem, ".tiff"))
  )

  ggplot2::ggsave(
    filename = output_files[["png"]],
    plot = plot, width = width, height = height, dpi = dpi,
    limitsize = FALSE
  )

  pdf_device <- if (isTRUE(base::capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  
  ggplot2::ggsave(
    filename = output_files[["pdf"]],
    plot = plot, width = width, height = height,
    device = pdf_device,
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename = output_files[["tiff"]],
    plot = plot, width = width, height = height, dpi = dpi,
    compression = "lzw", limitsize = FALSE
  )

  missing_files <- output_files[!file.exists(output_files)]
  if (length(missing_files)) {
    stop(
      "Figure export failed. Missing file(s): ",
      paste(missing_files, collapse = ", ")
    )
  }

  invisible(unname(output_files))
}

save_base_plot_all <- function(filename_stem, plot_fun,
                               width = 8, height = 7,
                               res = 600) {
  grDevices::png(
    paste0(filename_stem, ".png"),
    width = width * res, height = height * res, res = res
  )
  plot_fun()
  grDevices::dev.off()
  
  grDevices::tiff(
    paste0(filename_stem, ".tiff"),
    width = width * res, height = height * res,
    res = res, compression = "lzw"
  )
  plot_fun()
  grDevices::dev.off()
  
  grDevices::cairo_pdf(
    paste0(filename_stem, ".pdf"),
    width = width, height = height
  )
  plot_fun()
  grDevices::dev.off()
}

# -------------------- 10. Figures per outcome ---------------

make_figures <- function(outcome) {
  b <- readRDS(file.path(
    output_root, outcome, "biochar_meta_models.rds"
  ))
  
  out <- file.path(output_root, outcome)
  fig <- file.path(out, "figures")
  main_fig <- file.path(fig, "main_text")
  supp_fig <- file.path(fig, "supplement")
  
  dir.create(main_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(supp_fig, recursive = TRUE, showWarnings = FALSE)
  
  es <- b$data
  label <- b$outcome_label
  
  # Forest plot: comparison-level REML intervals and CR2 pooled CI.
  cr2 <- readr::read_csv(
    file.path(out, "overall_effect_CR2.csv"),
    show_col_types = FALSE
  ) |> slice(1)
  
  forest_fun <- function() {
    k <- nrow(es)
    rows <- (k + 4):5
    metafor::forest(
      x = es$yi,
      vi = es$vi,
      slab = es$comparison_id,
      rows = rows,
      ylim = c(0, k + 8),
      refline = 0,
      xlab = "Log response ratio (lnRR)",
      header = c("Comparison", "lnRR [95% CI]"),
      cex = 0.55
    )
    metafor::addpoly.default(
      x = as.numeric(cr2$estimate),
      ci.lb = as.numeric(cr2$ci_lb),
      ci.ub = as.numeric(cr2$ci_ub),
      rows = 3,
      mlab = "Overall estimate (article-clustered CR2 95% CI)",
      annotate = TRUE
    )
  }
  
  save_base_plot_all(
    file.path(fig, "forest_multilevel_lnRR_CR2"),
    forest_fun,
    width = 11,
    height = max(7, 0.18 * nrow(es) + 3)
  )
  
  # Article-level funnel plot.
  article_es <- b$article_level_effects
  article_uni <- b$article_aggregated_model
  
  if (!is.null(article_uni) && nrow(article_es) >= 2L) {
    funnel_fun <- function() {
      metafor::funnel(
        article_uni,
        yaxis = "sei",
        level = c(90, 95, 99),
        shade = c("white", "grey85", "grey70"),
        xlab = "Article-level log response ratio",
        ylab = "Standard error"
      )
    }
    
    save_base_plot_all(
      file.path(fig, "funnel_article_level"),
      funnel_fun, width = 7, height = 6
    )
  }
  
  # Feedstock composition.
  feedstock_counts <- es |>
    filter(!is.na(feedstock)) |>
    mutate(feedstock = factor(
      as.character(feedstock), levels = expected_feedstocks
    )) |>
    group_by(feedstock) |>
    summarise(
      k = n(), m = n_distinct(article_id), .groups = "drop"
    )
  
  p_feedstock <- ggplot(
    feedstock_counts,
    aes(x = feedstock, y = k)
  ) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = paste0("k = ", k, "\nm = ", m)),
      hjust = -0.1
    ) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = paste("Feedstock composition:", label),
      x = "Biochar feedstock type",
      y = "Number of effect sizes"
    ) +
    theme_pub()
  
  save_plot_all(
    p_feedstock,
    file.path(main_fig, "feedstock_type_distribution.png"),
    width = 8, height = 5
  )
  
  # CR2 categorical subgroup forest-style plot.
  sub <- b$categorical_terms |>
    filter(
      is.finite(estimate_pct),
      is.finite(pct_lb),
      is.finite(pct_ub)
    )
  
  if (nrow(sub)) {
    sub <- sub |>
      mutate(
        label = paste0(level, " (k=", k, ", m=", m, ")")
      )
    
    make_subgroup_plot <- function(dat, filename, title) {
      if (!nrow(dat)) return(invisible(NULL))
      
      p <- ggplot(
        dat,
        aes(x = estimate_pct, y = reorder(label, estimate_pct))
      ) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_errorbar(
          aes(xmin = pct_lb, xmax = pct_ub),
          orientation = "y", width = 0.15
        ) +
        geom_point(size = 2.8) +
        facet_wrap(~ moderator, scales = "free_y", ncol = 2) +
        labs(
          title = title,
          subtitle = paste(
            "Bars are article-clustered CR2-Satterthwaite 95% CIs."
          ),
          x = paste0("Change in ", label, " (%)"),
          y = NULL
        ) +
        theme_pub()
      
      save_plot_all(
        p, filename,
        width = 11, height = max(6, 0.35 * nrow(dat) + 3)
      )
      invisible(p)
    }
    
    main_sub <- sub |>
      filter(moderator %in% main_text_categorical)
    
    supp_sub <- sub |>
      filter(moderator %in% supplement_categorical)
    
    make_subgroup_plot(
      main_sub,
      file.path(main_fig, "categorical_subgroups_CR2.png"),
      paste("Main-text categorical subgroup estimates:", label)
    )
    
    make_subgroup_plot(
      supp_sub,
      file.path(supp_fig, "categorical_subgroups_CR2.png"),
      paste("Supplementary categorical subgroup estimates:", label)
    )
  }
  
  # Orchard plots (one per eligible categorical moderator).
  orchard_panels <- list()
  
  if (requireNamespace("orchaRd", quietly = TRUE)) {
    for (v in names(b$categorical_models)) {
      p <- tryCatch(
        orchaRd::orchard_plot(
          b$categorical_models[[v]],
          mod = v,
          group = "article_id",
          xlab = "Log response ratio (lnRR)"
        ) +
          labs(title = v) +
          theme_pub(11),
        error = function(e) {
          message("Orchard plot failed for ", v, ": ", e$message)
          NULL
        }
      )
      if (!is.null(p)) orchard_panels[[v]] <- p
    }
  }
  
  if (length(orchard_panels)) {
    main_panels <- orchard_panels[
      intersect(names(orchard_panels), main_text_categorical)
    ]
    supp_panels <- orchard_panels[
      intersect(names(orchard_panels), supplement_categorical)
    ]
    
    if (length(main_panels)) {
      p <- wrap_plots(main_panels, ncol = 2) +
        plot_annotation(
          title = paste("Main-text orchard plots:", label),
          subtitle = paste(
            "Model estimates shown; formal moderator inference",
            "is reported in the CR2 tables."
          )
        )
      save_plot_all(
        p, file.path(main_fig, "categorical_moderators_orchard.png"),
        width = 12, height = max(6, 4.5 * ceiling(length(main_panels) / 2))
      )
    }
    
    if (length(supp_panels)) {
      p <- wrap_plots(supp_panels, ncol = 2) +
        plot_annotation(title = paste("Supplementary orchard plots:", label))
      save_plot_all(
        p, file.path(supp_fig, "categorical_moderators_orchard.png"),
        width = 12, height = max(6, 4.5 * ceiling(length(supp_panels) / 2))
      )
    }
  }
  
  # Continuous meta-regression panels.
  plot_continuous <- function(v, fit) {
    spec <- continuous_spec |> filter(var == v)
    d <- es |> filter(is.finite(.data[[v]]))
    mu <- mean(d[[v]])
    sig <- sd(d[[v]])
    
    if (!is.finite(sig) || sig <= 0) return(NULL)
    
    xgrid <- seq(min(d[[v]]), max(d[[v]]), length.out = 200)
    pr <- predict(
      fit,
      newmods = matrix((xgrid - mu) / sig, ncol = 1)
    )
    
    is_log <- spec$transform == "log10"
    x_plot <- if (is_log) 10^xgrid else xgrid
    x_obs <- if (is_log) 10^d[[v]] else d[[v]]
    
    grid <- tibble(
      x = x_plot,
      estimate = lnrr_to_pct(pr$pred),
      lower = lnrr_to_pct(pr$ci.lb),
      upper = lnrr_to_pct(pr$ci.ub)
    )
    
    obs <- tibble(
      x = x_obs,
      y = lnrr_to_pct(d$yi),
      weight = 1 / d$vi
    )
    
    xlab <- paste0(
      spec$label,
      if (nzchar(spec$unit)) paste0(" (", spec$unit, ")") else ""
    )
    
    p <- ggplot(obs, aes(x = x, y = y)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(size = weight), alpha = 0.55) +
      geom_ribbon(
        data = grid,
        aes(x = x, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        alpha = 0.25
      ) +
      geom_line(
        data = grid,
        aes(x = x, y = estimate),
        inherit.aes = FALSE,
        linewidth = 0.8
      ) +
      scale_size_continuous(guide = "none") +
      labs(
        title = spec$label,
        x = xlab,
        y = paste0("Change in ", label, " (%)")
      ) +
      theme_pub(11)
    
    if (is_log) p <- p + scale_x_log10()
    p
  }
  
  make_continuous_set <- function(vars, destination, filename, title) {
    vars <- intersect(vars, names(b$continuous_models))
    panels <- purrr::compact(lapply(vars, function(v) {
      tryCatch(
        plot_continuous(v, b$continuous_models[[v]]),
        error = function(e) {
          message("Continuous plot failed for ", v, ": ", e$message)
          NULL
        }
      )
    }))
    
    if (!length(panels)) return(invisible(NULL))
    
    p <- wrap_plots(panels, ncol = 2) +
      plot_annotation(
        title = title,
        subtitle = paste(
          "Lines and ribbons are REML model predictions and model-based",
          "95% confidence intervals; they are not CR2 intervals.",
          "Use the exported CR2 tables for formal moderator inference."
        )
      )
    
    save_plot_all(
      p, file.path(destination, filename),
      width = 12, height = max(6, 5.5 * ceiling(length(panels) / 2))
    )
    invisible(p)
  }
  
  make_continuous_set(
    main_text_continuous, main_fig,
    "continuous_moderators_meta_regression.png",
    paste("Main-text continuous meta-regressions:", label)
  )
  
  make_continuous_set(
    supplement_continuous, supp_fig,
    "continuous_moderators_meta_regression.png",
    paste("Supplementary continuous meta-regressions:", label)
  )
  
  # Leave-one-article-out.
  loo <- b$loo
  if (nrow(loo)) {
    p_loo <- ggplot(
      loo,
      aes(x = pct, y = reorder(omitted_article, pct))
    ) +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_errorbar(
        aes(xmin = pct_lb, xmax = pct_ub),
        orientation = "y", width = 0.15
      ) +
      geom_point() +
      labs(
        title = "Leave-one-article-out sensitivity",
        x = paste0("Change in ", label, " (%)"),
        y = "Omitted Study_ID",
        subtitle = "Intervals are article-clustered CR2-Satterthwaite 95% CIs."
      ) +
      theme_pub()
    
    save_plot_all(
      p_loo, file.path(fig, "leave_one_article_out.png"),
      width = 8, height = max(6, 0.22 * nrow(loo) + 2)
    )
  }
  
  # Assumed within-article correlation sensitivity.
  rho <- b$rho_results
  if (nrow(rho)) {
    p_rho <- ggplot(
      rho,
      aes(x = factor(assumed_rho), y = pct)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_errorbar(
        aes(ymin = pct_lb, ymax = pct_ub), width = 0.12
      ) +
      geom_point(size = 2.7) +
      labs(
        title = "Within-article dependence sensitivity",
        subtitle = paste(
          "Rho values are assumed working correlations,",
          "not estimates of the true correlation."
        ),
        x = "Assumed within-article sampling correlation (rho)",
        y = paste0("Change in ", label, " (%)")
      ) +
      theme_pub()
    
    save_plot_all(
      p_rho,
      file.path(fig, "within_article_correlation_sensitivity.png"),
      width = 8, height = 5.5
    )
  }
  
  invisible(TRUE)
}

invisible(lapply(outcomes_requested, make_figures))

# -------------------- 11. Main-text manuscript figures -------

# One consolidated, publication-ready set of Figures 1-4.
# All figures are saved as PNG, PDF and 600-dpi LZW TIFF.
# CR2 intervals are used for pooled and categorical estimates.
# Continuous meta-regression ribbons are model-based REML CIs.

main_fig_dir <- file.path(output_root, "Main_Text_Figures")
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Consistent figure theme.
theme_main <- function(base_size = 12) {
  theme_pub(base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "grey30"),
      axis.text = element_text(color = "grey15"),
      strip.background = element_rect(fill = "grey95", colour = "grey65"),
      strip.text = element_text(face = "bold"),
      panel.spacing = grid::unit(1, "lines"),
      legend.position = "bottom"
    )
}

# Helper for consistent captions saved separately for manuscript use.
figure_captions <- tibble::tribble(
  ~figure, ~title, ~caption,
  "Figure 1",
  "Overall effects of biochar on N2O emissions and global warming potential",
  paste(
    "Pooled log response ratios were estimated using multilevel random-effects",
    "models fitted by REML, with article-level clustering and CR2-Satterthwaite",
    "inference. Points are back-transformed percentage changes relative to controls;",
    "horizontal intervals are article-clustered CR2 95% confidence intervals.",
    "Negative values indicate reductions relative to controls. k denotes usable",
    "effect-size comparisons and m denotes independent articles."
  ),
  "Figure 2",
  "Categorical subgroup estimates for feedstock and biochar application strategy",
  paste(
    "Subgroup estimates are shown for biochar feedstock and application strategy",
    "(biochar alone or combined with other amendments), separately for N2O and GWP.",
    "Points and horizontal intervals show percentage changes and article-clustered",
    "CR2-Satterthwaite 95% confidence intervals. k denotes comparisons and m denotes",
    "independent articles. Overall moderator tests, rather than visual CI overlap,",
    "should be used to assess differences among subgroups."
  ),
  "Figure 3",
  "Categorical subgroup estimates by climate zone and soil texture",
  paste(
    "Subgroup estimates are shown for climate zone and soil texture, separately",
    "for N2O and GWP. Points and horizontal intervals show percentage changes and",
    "article-clustered CR2-Satterthwaite 95% confidence intervals. k denotes",
    "comparisons and m denotes independent articles. Sparse subgroups and low",
    "denominator degrees of freedom should be interpreted cautiously."
  ),
  "Figure 4",
  "Continuous moderator meta-regressions for N2O emissions and global warming potential",
  paste(
    "Panels show pyrolysis temperature, biochar application rate, trial duration,",
    "and initial soil pH. Points represent observed comparison-level percentage",
    "changes. Lines and shaded ribbons show REML meta-regression predictions and",
    "model-based 95% confidence intervals; these are not CR2 confidence intervals.",
    "Each moderator was fitted separately. Robust omnibus inference is reported",
    "in Supplementary Table S5."
  )
)

# Load the reproducibility bundles created by the outcome analyses.
bundles <- setNames(
  lapply(outcomes_requested, function(o) {
    readRDS(file.path(output_root, o, "biochar_meta_models.rds"))
  }),
  outcomes_requested
)

# Figure 1. Overall effects with k and m.
overall_df <- purrr::map_dfr(names(bundles), function(o) {
  bb <- bundles[[o]]
  z <- bb$main_cr2 |>
    dplyr::slice(1)

  tibble::tibble(
    outcome = o,
    estimate = z$estimate[1],
    ci_lb = z$ci_lb[1],
    ci_ub = z$ci_ub[1],
    k = nrow(bb$data),
    m = dplyr::n_distinct(bb$data$article_id)
  )
}) |>
  mutate(
    outcome = factor(outcome, levels = c("N2O", "GWP")),
    estimate_pct = lnrr_to_pct(estimate),
    lower_pct = lnrr_to_pct(ci_lb),
    upper_pct = lnrr_to_pct(ci_ub),
    estimate_label = sprintf("%.1f%%  (k=%d; m=%d)", estimate_pct, k, m)
  )

p1 <- ggplot(overall_df, aes(x = estimate_pct, y = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5,
             colour = "grey45") +
  geom_errorbar(
    aes(xmin = lower_pct, xmax = upper_pct),
    orientation = "y", width = 0.16, linewidth = 0.7
  ) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = estimate_label),
    nudge_y = 0.22, hjust = 0, size = 3.4
  ) +
  scale_y_discrete(labels = c(
    N2O = "N2O emissions",
    GWP = "Global warming potential"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.32))) +
  labs(
    title = "Figure 1. Overall biochar effects",
    subtitle = "Article-clustered CR2-Satterthwaite 95% confidence intervals",
    x = "Change relative to control (%)",
    y = NULL
  ) +
  theme_main(12)

save_plot_all(
  p1, file.path(main_fig_dir, "Figure_1_overall_effects"),
  width = 9, height = 5, dpi = 600
)

# Figures 2 and 3. Split categorical moderators to match Results 3.3.
# Figure 2: feedstock and treatment type. Figure 3: climate and soil texture.
subgroup_both <- purrr::map_dfr(names(bundles), function(o) {
  bundles[[o]]$categorical_terms |>
    mutate(outcome = o)
}) |>
  filter(
    moderator %in% main_text_categorical,
    is.finite(estimate_pct),
    is.finite(pct_lb),
    is.finite(pct_ub)
  ) |>
  mutate(
    outcome = factor(outcome, levels = c("N2O", "GWP")),
    moderator_label = dplyr::recode(
      moderator,
      feedstock = "Feedstock",
      treatment_type = "Treatment type",
      climate_zone = "Climate zone",
      texture = "Soil texture"
    ),
    display_level = paste0(level, " (k=", k, "; m=", m, ")")
  )

make_categorical_figure <- function(data, moderators, figure_title,
                                    output_stub, width = 12) {
  moderator_labels <- c(
    feedstock = "Feedstock",
    treatment_type = "Treatment type",
    climate_zone = "Climate zone",
    texture = "Soil texture"
  )
  plot_data <- data |>
    filter(moderator %in% moderators) |>
    mutate(
      moderator_label = factor(
        moderator_label,
        levels = unname(moderator_labels[moderators])
      )
    )

  if (!nrow(plot_data)) {
    message(figure_title, " skipped: no eligible subgroup estimates.")
    return(invisible(NULL))
  }

  p <- ggplot(
    plot_data,
    aes(x = estimate_pct, y = reorder(display_level, estimate_pct))
  ) +
    geom_vline(
      xintercept = 0, linetype = "dashed",
      linewidth = 0.45, colour = "grey45"
    ) +
    geom_errorbar(
      aes(xmin = pct_lb, xmax = pct_ub),
      orientation = "y", width = 0.14, linewidth = 0.55
    ) +
    geom_point(size = 2.2) +
    facet_grid(
      outcome ~ moderator_label,
      scales = "free_y", space = "free_y"
    ) +
    labs(
      title = figure_title,
      subtitle = paste(
        "Article-clustered CR2-Satterthwaite 95% CIs;",
        "k = comparisons, m = independent articles"
      ),
      x = "Change relative to control (%)",
      y = NULL
    ) +
    theme_main(10) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9),
      panel.spacing = grid::unit(0.8, "lines")
    )

  save_plot_all(
    p, file.path(main_fig_dir, output_stub),
    width = width,
    height = max(6.5, 1.05 * dplyr::n_distinct(plot_data$display_level)),
    dpi = 600
  )
  invisible(p)
}

p2 <- make_categorical_figure(
  subgroup_both,
  moderators = c("feedstock", "treatment_type"),
  figure_title = "Figure 2. Feedstock and biochar application strategy",
  output_stub = "Figure_2_feedstock_and_treatment"
)

p3 <- make_categorical_figure(
  subgroup_both,
  moderators = c("climate_zone", "texture"),
  figure_title = "Figure 3. Climate zone and soil texture",
  output_stub = "Figure_3_climate_and_soil_texture"
)

# Figure 4. Continuous moderator meta-regressions.
continuous_plot_list <- list()

for (oo in names(bundles)) {
  bb <- bundles[[oo]]

  for (v in intersect(main_text_continuous, names(bb$continuous_models))) {
    fit <- bb$continuous_models[[v]]
    dd <- bb$data |>
      filter(is.finite(.data[[v]]))

    if (nrow(dd) < 3L) next

    mu <- mean(dd[[v]], na.rm = TRUE)
    sig <- sd(dd[[v]], na.rm = TRUE)
    if (!is.finite(sig) || sig <= 0) next

    gx <- seq(min(dd[[v]], na.rm = TRUE),
              max(dd[[v]], na.rm = TRUE), length.out = 180)

    pp <- tryCatch(
      predict(fit, newmods = matrix((gx - mu) / sig, ncol = 1)),
      error = function(e) {
        message("Prediction failed for ", oo, " / ", v, ": ", e$message)
        NULL
      }
    )
    if (is.null(pp)) next

    sp <- continuous_spec |> filter(var == v)
    is_log <- identical(sp$transform[1], "log10")
    x_label <- paste0(
      sp$label[1],
      ifelse(is.na(sp$unit[1]) | sp$unit[1] == "",
             "", paste0(" (", sp$unit[1], ")"))
    )

    gd <- tibble(
      x = if (is_log) 10^gx else gx,
      pred = lnrr_to_pct(as.numeric(pp$pred)),
      lo = lnrr_to_pct(as.numeric(pp$ci.lb)),
      hi = lnrr_to_pct(as.numeric(pp$ci.ub))
    )

    od <- tibble(
      x = if (is_log) 10^dd[[v]] else dd[[v]],
      y = lnrr_to_pct(dd$yi),
      weight = 1 / dd$vi
    )

    panel <- ggplot(od, aes(x = x, y = y)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 linewidth = 0.45, colour = "grey45") +
      geom_ribbon(
        data = gd,
        aes(x = x, ymin = lo, ymax = hi),
        inherit.aes = FALSE, alpha = 0.20
      ) +
      geom_line(
        data = gd, aes(x = x, y = pred),
        inherit.aes = FALSE, linewidth = 0.8
      ) +
      geom_point(aes(size = weight), alpha = 0.45) +
      scale_size_continuous(guide = "none") +
      labs(
        title = paste0(oo, ": ", sp$label[1]),
        x = x_label,
        y = "Change relative to control (%)"
      ) +
      theme_main(10)

    if (is_log) panel <- panel + scale_x_log10()
    continuous_plot_list[[paste(oo, v, sep = "_")]] <- panel
  }
}

if (length(continuous_plot_list) > 0) {
  p4 <- patchwork::wrap_plots(continuous_plot_list, ncol = 2) +
    patchwork::plot_annotation(
      title = "Figure 4. Continuous moderator meta-regressions",
      subtitle = paste(
        "Lines and ribbons are REML model predictions and model-based 95% CIs;",
        "they are not CR2 confidence intervals."
      ),
      caption = paste(
        "Each moderator was fitted separately. Points represent observed",
        "comparison-level effects; point size reflects inverse sampling variance."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        plot.caption = element_text(size = 8, hjust = 0)
      )
    )

  save_plot_all(
    p4, file.path(main_fig_dir, "Figure_4_continuous_moderators"),
    width = 14,
    height = max(8, 4.8 * ceiling(length(continuous_plot_list) / 2)),
    dpi = 600
  )
} else {
  message("Figure 4 skipped: no eligible continuous main-text moderator models.")
}

# Verify that every requested main-text figure was exported in all formats.
expected_main_figures <- c(
  "Figure_1_overall_effects",
  "Figure_2_feedstock_and_treatment",
  "Figure_3_climate_and_soil_texture",
  "Figure_4_continuous_moderators"
)

expected_main_files <- unlist(lapply(
  expected_main_figures,
  function(stem) file.path(
    main_fig_dir, paste0(stem, c(".png", ".pdf", ".tiff"))
  )
), use.names = FALSE)

missing_main_files <- expected_main_files[!file.exists(expected_main_files)]

if (length(missing_main_files)) {
  stop(
    "One or more required main-text figures were not exported. ",
    "Check moderator eligibility and the messages above. Missing files:\n",
    paste(missing_main_files, collapse = "\n")
  )
}

# Export captions and a manifest of the files generated.
write_csv_safe(
  figure_captions,
  file.path(main_fig_dir, "Main_text_figure_captions.csv")
)

figure_files <- list.files(
  main_fig_dir,
  pattern = "^Figure_[1-4]_.*\\.(png|pdf|tiff)$",
  full.names = FALSE
)

write_csv_safe(
  tibble(
    file = figure_files,
    format = tools::file_ext(figure_files),
    dpi = ifelse(tools::file_ext(figure_files) == "pdf",
                 NA_integer_, 600L),
    note = ifelse(tools::file_ext(figure_files) == "pdf",
                  "Vector PDF", "600-dpi raster; TIFF uses LZW compression")
  ),
  file.path(main_fig_dir, "Main_text_figure_file_index.csv")
)

# -------------------- 12. Supplementary figures S2-S6 --------

# Supplementary Figure S1 is the PRISMA flow diagram and must be supplied
# from the original screening records.

# S4: article-level funnel plots for small-study-effect diagnostics.
funnel_data <- bind_rows(lapply(names(bundles), function(o) {
  bundles[[o]]$article_level_effects |>
    mutate(outcome = o)
}))
funnel_centres <- bind_rows(lapply(names(bundles), function(o) {
  b <- bundles[[o]]
  tibble::tibble(
    outcome = o,
    centre = if (!is.null(b$article_aggregated_model)) {
      as.numeric(coef(b$article_aggregated_model)[1])
    } else NA_real_
  )
}))

if (nrow(funnel_data) > 0) {
  funnel_data <- funnel_data |>
    left_join(funnel_centres, by = "outcome") |>
    filter(is.finite(yi), is.finite(sei), is.finite(centre))

  funnel_limits <- funnel_data |>
    distinct(outcome, centre, sei) |>
    tidyr::crossing(boundary = c("lower", "upper")) |>
    mutate(
      boundary_x = if_else(
        boundary == "lower",
        centre - 1.96 * sei,
        centre + 1.96 * sei
      )
    )

  pS4 <- ggplot(funnel_data, aes(x = yi, y = sei)) +
    geom_line(
      data = funnel_limits,
      aes(x = boundary_x, y = sei, group = interaction(outcome, boundary)),
      inherit.aes = FALSE, colour = "grey65", linewidth = 0.5
    ) +
    geom_vline(
      data = funnel_centres,
      aes(xintercept = centre),
      linetype = "dashed", inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.75, size = 2) +
    scale_y_reverse() +
    facet_wrap(
      ~ outcome, scales = "free",
      labeller = as_labeller(c(
        N2O = "N2O emissions",
        GWP = "Global warming potential"
      ))
    ) +
    labs(
      title = "Supplementary Figure S4. Article-level funnel plots",
      subtitle = paste(
        "Each point represents an article-level inverse-variance aggregate;",
        "lines show approximate 95% pseudo-confidence limits."
      ),
      x = "Article-level log response ratio",
      y = "Standard error"
    ) +
    theme_pub()

  save_plot_all(
    pS4, file.path(output_root, "Supplementary_Figure_S4"),
    width = 9, height = 5.5, dpi = 600
  )
}

# S5: study-level evidence base by outcome.
evidence_both <- bind_rows(lapply(outcomes_requested, function(o) {
  readr::read_csv(
    file.path(output_root, o, "evidence_base_counts.csv"),
    show_col_types = FALSE
  )
}))

pS1 <- ggplot(
  evidence_both,
  aes(x = outcome, y = k_comparisons)
) +
  geom_col() +
  labs(
    title = "Supplementary Figure S5. Evidence base",
    x = "Outcome", y = "Usable comparisons"
  ) +
  theme_pub()

save_plot_all(
  pS1, file.path(output_root, "Supplementary_Figure_S5"),
  width = 7, height = 5
)

# S2: leave-one-article-out sensitivity, both outcomes.
loo_both <- bind_rows(lapply(outcomes_requested, function(o) {
  readr::read_csv(
    file.path(output_root, o, "leave_one_article_out.csv"),
    show_col_types = FALSE
  ) |>
    mutate(outcome = o)
}))

if (nrow(loo_both)) {
  pS2 <- ggplot(
    loo_both,
    aes(x = pct, y = reorder(omitted_article, pct))
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbar(
      aes(xmin = pct_lb, xmax = pct_ub),
      orientation = "y", width = 0.12
    ) +
    geom_point() +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(
      title = "Supplementary Figure S2. Leave-one-article-out analysis",
      x = "Change (%)", y = "Omitted article"
    ) +
    theme_pub()
  
  save_plot_all(
    pS2, file.path(output_root, "Supplementary_Figure_S2.png"),
    width = 11, height = 8
  )
}

# S3: assumed within-article correlation sensitivity.
rho_both <- bind_rows(lapply(outcomes_requested, function(o) {
  readr::read_csv(
    file.path(output_root, o, "within_article_correlation_sensitivity.csv"),
    show_col_types = FALSE
  ) |>
    mutate(outcome = o)
}))

if (nrow(rho_both)) {
  pS3 <- ggplot(
    rho_both,
    aes(x = assumed_rho, y = pct)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(
      aes(ymin = pct_lb, ymax = pct_ub),
      width = 0.04
    ) +
    geom_point() +
    facet_wrap(~ outcome) +
    labs(
      title = "Supplementary Figure S3. Dependence sensitivity",
      subtitle = "Rho values are assumed, not estimated correlations.",
      x = "Assumed within-article rho",
      y = "Change (%)"
    ) +
    theme_pub()
  
  save_plot_all(
    pS3, file.path(output_root, "Supplementary_Figure_S3.png"),
    width = 9, height = 5.5
  )
}

# S6: moderator support diagnostics.
support_both <- bind_rows(lapply(outcomes_requested, function(o) {
  bb <- readRDS(file.path(
    output_root, o, "biochar_meta_models.rds"
  ))
  bind_rows(
    bb$continuous_support |> mutate(type = "Continuous"),
    bb$categorical_support |> mutate(type = "Categorical")
  ) |>
    mutate(outcome = o)
}))

if (nrow(support_both)) {
  pS4 <- ggplot(
    support_both,
    aes(x = moderator, y = k, shape = support_status)
  ) +
    geom_point(size = 2.8) +
    facet_wrap(~ outcome, scales = "free_x") +
    coord_flip() +
    labs(
      title = "Supplementary Figure S6. Moderator support",
      x = "Moderator", y = "Effect sizes (k)"
    ) +
    theme_pub()
  
  save_plot_all(
    pS4, file.path(output_root, "Supplementary_Figure_S6"),
    width = 10, height = 8
  )
}

# -------------------- 13. Final reproducibility record -------

capture.output(
  sessionInfo(),
  file = file.path(output_root, "sessionInfo.txt")
)

writeLines(
  c(
    "Analysis completed.",
    paste("Outcomes analysed:", paste(outcomes_requested, collapse = ", ")),
    "SOC was not analysed as an outcome.",
    "soc_initial was retained as a continuous moderator.",
    "Primary inference: article-clustered CR2-Satterthwaite.",
    "Moderator model fitting: one moderator at a time.",
    "Nested fixed-effect comparisons: ML.",
    "Meta-regression ribbons: model-based REML confidence intervals.",
    "Within-article rho sensitivity: assumed working covariance structures.",
    "Table S1 search strategy is a template and must be completed from search records.",
    "Supplementary Figure S1 is the PRISMA flow diagram and must be supplied from screening records.",
    "Supplementary Figures S2-S3 show leave-one-article-out and assumed-rho sensitivity.",
    "Supplementary Figure S4 shows article-level funnel plots; asymmetry is not proof of publication bias.",
    "Supplementary Figures S5-S6 show evidence-base and moderator-support summaries.",
    "Table S8 is divided into S8a sensitivity analyses and S8b exploratory small-study-effect diagnostics.",
    "PRISMA counts and search records must be verified against original screening records."
  ),
  file.path(output_root, "RUN_NOTES.txt")
)

message(
  "\nAnalysis and figure workflows completed.\nOutputs: ",
  normalizePath(output_root)
)