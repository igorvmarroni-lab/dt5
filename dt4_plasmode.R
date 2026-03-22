# =============================================================================
# Digital Twin 4.0 (DT4) — Plasmode Simulation · Causal Estimator Validation
# Marroni 2026 | BD-II Primary Prevention
# =============================================================================
# Design (Franklin et al. 2014 Plasmode framework):
#   Real data substrate: MEPS longitudinal cohort
#     - N=28,052 non-cases aged 15-30 (PHQ2<3 at baseline)
#     - Eligible analytic cohort: K6SUM>=7 (N=3,398)
#     - Outcome: PHQ2 transition PHQ2<3 -> PHQ2>=3 within 1 year
#     - Empirical transition rate: 14.9% in K6SUM>=7 stratum
#   Injected causal effect: RR_true = 0.741 = exp(-0.30)
#   Per seed (N=200):
#     1. Resample eligible cohort with replacement (preserve joint distribution)
#     2. Randomise 1:1 treatment/control
#     3. Inject causal effect: P(Y=1|Tx) = P(Y=1|Ctrl) * exp(-0.30)
#     4. Generate potential outcomes from Binomial
#     5. Estimate RR from observed outcomes
#   Primary metric: |RR_estimated - RR_true| (absolute bias)
#
# IMPORTANT: MEPS data must be downloaded from https://meps.ahrq.gov
#   and preprocessed into the format described in DATA_SOURCES.md.
#   A synthetic substitute is provided below for validation purposes.
# =============================================================================

library(tidyverse)

# ── Parameters ───────────────────────────────────────────────────────────────

RR_TRUE        <- 0.741           # = exp(-0.30)
LOG_RR_TRUE    <- log(RR_TRUE)    # -0.2999...
N_SEEDS        <- 200
BASE_RATE      <- 0.149           # empirical PHQ2 transition rate (K6SUM>=7)
N_ELIGIBLE     <- 3398            # MEPS K6SUM>=7 analytic cohort

# ── Load or simulate MEPS substrate ─────────────────────────────────────────
# In production: replace load_meps_data() with your real preprocessed MEPS file
# Required columns:
#   - phq2_baseline  : PHQ2 score at baseline (0-6; <3 required for eligibility)
#   - phq2_followup  : PHQ2 score at 1-year follow-up (0-6)
#   - k6sum          : Kessler-6 sum score
#   - age            : age at baseline (15-30)
#   - sex            : factor
#   - panel_year     : MEPS panel year (1996-2023)
#   - outcome        : binary, PHQ2<3 at baseline -> PHQ2>=3 at follow-up

load_meps_data <- function(path = "data/meps/meps_analytic_cohort.csv") {
  if (file.exists(path)) {
    df <- read.csv(path)
    message("Loaded real MEPS data: N=", nrow(df))
    return(df)
  }

  # ── Synthetic substitute (for code validation only) ──────────────────────
  # Preserves: N, marginal outcome rate, K6SUM>=7 eligibility criterion
  # Does NOT reproduce real MEPS joint covariate distribution
  message("MEPS data file not found at '", path, "'.")
  message("Using SYNTHETIC SUBSTITUTE for code validation only.")
  message("Results with synthetic data are NOT the paper's reported estimates.")
  message("Download real MEPS data from https://meps.ahrq.gov")

  set.seed(2026)
  n <- N_ELIGIBLE

  age        <- round(runif(n, 15, 30))
  sex        <- sample(c("M", "F"), n, replace = TRUE)
  panel_year <- sample(1996:2023, n, replace = TRUE)
  k6sum      <- round(rnorm(n, mean = 10, sd = 3))
  k6sum      <- pmax(k6sum, 7)   # K6SUM>=7 eligibility

  phq2_baseline <- sample(0:2, n, replace = TRUE)

  # Outcome: ~14.9% base rate
  phq2_followup_binary <- rbinom(n, 1, BASE_RATE)
  phq2_followup <- ifelse(phq2_followup_binary == 1,
                          sample(3:6, n, replace = TRUE),
                          sample(0:2, n, replace = TRUE))

  data.frame(
    id             = seq_len(n),
    age            = age,
    sex            = sex,
    panel_year     = panel_year,
    k6sum          = k6sum,
    phq2_baseline  = phq2_baseline,
    phq2_followup  = phq2_followup,
    outcome        = as.integer(phq2_followup >= 3)
  )
}

# ── Single Plasmode seed ──────────────────────────────────────────────────────

run_dt4_seed <- function(seed, meps_df) {
  set.seed(seed)

  # 1. Resample with replacement (preserves joint covariate distribution)
  boot_idx <- sample(nrow(meps_df), nrow(meps_df), replace = TRUE)
  df       <- meps_df[boot_idx, ]

  # 2. Randomise 1:1
  n        <- nrow(df)
  df$arm   <- sample(c("Treatment", "Control"), n, replace = TRUE)

  # 3. Inject causal effect
  #    P(Y=1 | Treatment) = P(Y=1 | Control) * exp(log_rr_true)
  #    Use each individual's empirical outcome as P(Y=1 | Control)
  p_control   <- df$outcome
  p_treatment <- pmin(p_control * exp(LOG_RR_TRUE), 1)

  # 4. Generate potential outcomes
  y_control   <- rbinom(n, 1, p_control)
  y_treatment <- rbinom(n, 1, p_treatment)

  df$y_observed <- ifelse(df$arm == "Treatment", y_treatment, y_control)

  # 5. Estimate RR
  r_tx   <- mean(df$y_observed[df$arm == "Treatment"])
  r_ctrl <- mean(df$y_observed[df$arm == "Control"])

  rr_est <- r_tx / r_ctrl
  rd_est <- r_ctrl - r_tx

  list(
    seed       = seed,
    rr_true    = RR_TRUE,
    rr_est     = rr_est,
    rd_est     = rd_est,
    abs_bias   = abs(rr_est - RR_TRUE),
    rr_lt1     = as.integer(rr_est < 1.0),
    r_ctrl     = r_ctrl,
    r_tx       = r_tx
  )
}

# ── Bootstrap over 200 seeds ─────────────────────────────────────────────────

run_dt4_bootstrap <- function(meps_df, n_seeds = N_SEEDS, verbose = TRUE) {
  message("Running DT4 | Plasmode | N=", nrow(meps_df), " | ", n_seeds, " seeds")

  out <- lapply(seq_len(n_seeds), function(s) {
    if (verbose && s %% 50 == 0) message("  seed ", s, "/", n_seeds)
    run_dt4_seed(seed = s, meps_df = meps_df)
  })

  bind_rows(lapply(out, as.data.frame))
}

# ── Summary statistics ────────────────────────────────────────────────────────

summarise_dt4 <- function(results_df) {
  results_df %>%
    summarise(
      n_seeds          = n(),
      rr_true          = unique(rr_true),
      mean_rr_est      = mean(rr_est, na.rm = TRUE),
      median_rr_est    = median(rr_est, na.rm = TRUE),
      sd_rr_est        = sd(rr_est, na.rm = TRUE),
      mean_abs_bias    = mean(abs_bias, na.rm = TRUE),
      max_abs_bias     = max(abs_bias, na.rm = TRUE),
      pct_rr_lt1       = mean(rr_lt1, na.rm = TRUE) * 100,
      mean_rd          = mean(rd_est, na.rm = TRUE),
      pct_rd_positive  = mean(rd_est > 0, na.rm = TRUE) * 100
    )
}

# ── Main execution ────────────────────────────────────────────────────────────

main_dt4 <- function() {
  message("=== DT4: Plasmode Causal Validation ===")
  meps_df    <- load_meps_data()
  results_df <- run_dt4_bootstrap(meps_df, n_seeds = N_SEEDS)
  summary_df <- summarise_dt4(results_df)

  message("\n>>> RR_true = ", RR_TRUE)
  message(">>> Mean RR estimated = ", round(summary_df$mean_rr_est, 4))
  message(">>> Mean absolute bias = ", round(summary_df$mean_abs_bias, 4),
          " (", round(summary_df$mean_abs_bias / RR_TRUE * 100, 2), "%)")
  message(">>> RD negative in ", round(summary_df$pct_rd_positive, 1),
          "% of seeds (i.e., treatment protective)")
  message(">>> RR < 1.0 in ", round(summary_df$pct_rr_lt1, 1), "% of seeds")

  dir.create("DT4/results", recursive = TRUE, showWarnings = FALSE)
  write.csv(results_df, "DT4/results/dt4_seeds_all.csv",   row.names = FALSE)
  write.csv(summary_df, "DT4/results/dt4_summary_bias.csv", row.names = FALSE)
  message("\nResults written to DT4/results/")

  list(results = results_df, summary = summary_df)
}

# Run
dt4_output <- main_dt4()
print(dt4_output$summary)
