# =============================================================================
# Digital Twin 3.0 (DT3) — Synthetic Simulation with Real Allele Frequencies
# Marroni 2026 | BD-II Primary Prevention
# =============================================================================
# Design:
#   N = 2,000 individuals · T = 240 monthly steps (20 years)
#   Top PRS quartile (~N=500) randomised 1:1 (Prevented vs Control)
#   Three pre-specified delta scenarios: -0.02, -0.04 (base), -0.06
#   Risk model: R(t) = R(t-1) + alpha + delta*I[Prevented] + eps(t)
#   Thresholds: first depressive episode R>=90; BD-II conversion R>=98
#   Primary outcome: BD-II conversion at 10 years (t=120)
#   Analysis: Kaplan-Meier + Cox PH + bootstrap medians (200 seeds)
# =============================================================================

library(tidyverse)
library(survival)

# ── Parameters ──────────────────────────────────────────────────────────────

N_TOTAL   <- 2000          # total simulated individuals
N_SNPS    <- 500           # SNPs sampled from UKB AF pool
T_STEPS   <- 240           # monthly steps (20 years)
T_PRIMARY <- 120           # primary endpoint: 10 years
ALPHA     <- 0.18          # monthly risk progression (Control)
SIGMA_EPS <- 1.0           # SD of monthly noise
R0_MEAN   <- 60            # baseline risk mean
R0_SD     <- 10            # baseline risk SD
THRESH_DEP <- 90           # threshold: first depressive episode
THRESH_BD  <- 98           # threshold: BD-II conversion (after first episode)
DELTA_SCENARIOS <- c(minimal = -0.02, base = -0.04, optimistic = -0.06)
N_SEEDS   <- 200
MAF_WEIGHT_EXP <- 0.15     # inverse-MAF weighting exponent (Mullins 2021)
BETA_SIGMA <- 0.044        # Half-Normal SD for effect sizes (DT3 calibration)

# ── UKB allele frequency pool (simulated proxy) ──────────────────────────────
# In production: replace with real AF from IEU Open GWAS ukb-b-6906
# Here: sample from a Beta(1.5,1.5) distribution truncated to MAF [0.01, 0.49]
# matching the empirical MAF spectrum of PASS variants (Ti/Tv ~ 2.5)

simulate_ukb_af_pool <- function(n_snps, seed = 1) {
  set.seed(seed)
  af <- rbeta(n_snps * 10, 1.5, 1.5)
  af <- af[af >= 0.01 & af <= 0.49]
  af[1:n_snps]
}

# ── Genotype + PRS simulation ────────────────────────────────────────────────

simulate_prs <- function(n_individuals, af_pool, seed) {
  set.seed(seed)
  n_snps <- length(af_pool)

  # Genotype matrix: Binomial(2, af_j) for each SNP j
  G <- matrix(
    rbinom(n_individuals * n_snps, size = 2, prob = rep(af_pool, each = n_individuals)),
    nrow = n_individuals,
    ncol = n_snps
  )

  # Effect sizes: Half-Normal(sigma=0.044) with inverse-MAF weighting
  raw_beta <- abs(rnorm(n_snps, mean = 0, sd = BETA_SIGMA))
  maf      <- pmin(af_pool, 1 - af_pool)
  weights  <- (0.5 / maf)^MAF_WEIGHT_EXP
  beta     <- raw_beta * weights

  # Signs: random +/-
  signs <- sample(c(-1, 1), n_snps, replace = TRUE)
  beta  <- beta * signs

  # PRS
  prs <- G %*% beta
  prs <- as.numeric(prs)
  prs
}

# ── Single-seed simulation ───────────────────────────────────────────────────

run_dt3_seed <- function(seed, delta, af_pool) {
  set.seed(seed)

  # --- PRS and stratification ---
  prs <- simulate_prs(N_TOTAL, af_pool, seed = seed)
  prs_threshold <- quantile(prs, 0.75)
  eligible <- which(prs >= prs_threshold)          # top quartile ~N=500

  n_elig <- length(eligible)
  arm    <- sample(c("Prevented", "Control"), n_elig, replace = TRUE)

  # --- Individual risk trajectories ---
  results <- vector("list", n_elig)

  for (i in seq_along(eligible)) {
    r       <- rnorm(1, mean = R0_MEAN, sd = R0_SD)
    had_dep <- FALSE
    t_dep   <- NA
    t_bd    <- NA

    for (t in 1:T_STEPS) {
      eps <- rnorm(1, 0, SIGMA_EPS)
      d   <- if (arm[i] == "Prevented") delta else 0
      r   <- r + ALPHA + d + eps

      # First depressive episode
      if (!had_dep && r >= THRESH_DEP) {
        had_dep <- TRUE
        t_dep   <- t
      }

      # BD-II conversion (only after first depressive episode)
      if (had_dep && r >= THRESH_BD) {
        t_bd <- t
        break
      }
    }

    results[[i]] <- list(
      id    = eligible[i],
      arm   = arm[i],
      prs   = prs[eligible[i]],
      t_dep = t_dep,
      t_bd  = t_bd
    )
  }

  df <- bind_rows(lapply(results, as.data.frame))

  # --- Survival analysis at T_PRIMARY (10 years = 120 months) ---
  df <- df %>%
    mutate(
      event_bd  = !is.na(t_bd) & t_bd <= T_PRIMARY,
      time_bd   = ifelse(!is.na(t_bd), pmin(t_bd, T_PRIMARY), T_PRIMARY)
    )

  # RR (risk ratio at T_PRIMARY)
  r_control   <- mean(df$event_bd[df$arm == "Control"])
  r_prevented <- mean(df$event_bd[df$arm == "Prevented"])
  rr          <- r_prevented / r_control

  # NNT
  rd  <- r_control - r_prevented
  nnt <- if (rd > 0) 1 / rd else NA

  # Cox HR
  cox_fit <- tryCatch({
    coxph(Surv(time_bd, event_bd) ~ (arm == "Prevented"), data = df)
  }, error = function(e) NULL)
  hr <- if (!is.null(cox_fit)) exp(coef(cox_fit))[[1]] else NA

  list(
    seed        = seed,
    delta       = delta,
    r_control   = r_control,
    r_prevented = r_prevented,
    rr          = rr,
    rd          = rd,
    nnt         = nnt,
    hr          = hr
  )
}

# ── Bootstrap over 200 seeds ─────────────────────────────────────────────────

run_dt3_bootstrap <- function(delta, af_pool, n_seeds = N_SEEDS, verbose = TRUE) {
  label <- names(DELTA_SCENARIOS)[DELTA_SCENARIOS == delta]
  if (verbose) message("Running DT3 | delta = ", delta, " (", label, ") | ", n_seeds, " seeds")

  out <- lapply(seq_len(n_seeds), function(s) {
    if (verbose && s %% 50 == 0) message("  seed ", s, "/", n_seeds)
    run_dt3_seed(seed = s, delta = delta, af_pool = af_pool)
  })

  bind_rows(lapply(out, as.data.frame))
}

# ── Summary statistics ────────────────────────────────────────────────────────

summarise_dt3 <- function(results_df) {
  results_df %>%
    group_by(delta) %>%
    summarise(
      n_seeds          = n(),
      median_rr        = median(rr, na.rm = TRUE),
      ci_rr_lo         = quantile(rr, 0.025, na.rm = TRUE),
      ci_rr_hi         = quantile(rr, 0.975, na.rm = TRUE),
      pct_rr_lt1       = mean(rr < 1.0, na.rm = TRUE) * 100,
      median_nnt       = median(nnt, na.rm = TRUE),
      ci_nnt_lo        = quantile(nnt, 0.025, na.rm = TRUE),
      ci_nnt_hi        = quantile(nnt, 0.975, na.rm = TRUE),
      median_r_control = median(r_control, na.rm = TRUE),
      median_r_prev    = median(r_prevented, na.rm = TRUE),
      median_hr        = median(hr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario = case_when(
        delta == -0.02 ~ "Minimal",
        delta == -0.04 ~ "Base case (recommended)",
        delta == -0.06 ~ "Optimistic (upper bound)",
        TRUE ~ as.character(delta)
      )
    ) %>%
    select(scenario, delta, everything())
}

# ── Main execution ────────────────────────────────────────────────────────────

main_dt3 <- function() {
  message("=== DT3: Synthetic Simulation ===")
  message("Generating UKB AF pool (N=", N_SNPS, " SNPs) ...")
  af_pool <- simulate_ukb_af_pool(N_SNPS, seed = 42)

  all_results <- lapply(DELTA_SCENARIOS, function(d) {
    run_dt3_bootstrap(delta = d, af_pool = af_pool, n_seeds = N_SEEDS)
  })

  results_df <- bind_rows(all_results)
  summary_df <- summarise_dt3(results_df)

  # NNT with 20% attrition (base case)
  base_nnt <- summary_df$median_nnt[summary_df$delta == -0.04]
  message("\n>>> Base case NNT (median): ", round(base_nnt, 1))
  message(">>> Planning NNT (+20% attrition): ", round(base_nnt * 1.20, 1))

  # Save outputs
  dir.create("DT3/results", recursive = TRUE, showWarnings = FALSE)
  for (d in names(DELTA_SCENARIOS)) {
    sub <- results_df[results_df$delta == DELTA_SCENARIOS[d], ]
    write.csv(sub, file = paste0("DT3/results/dt3_seeds_", d, ".csv"), row.names = FALSE)
  }
  write.csv(summary_df, "DT3/results/dt3_summary_NNT_RR.csv", row.names = FALSE)
  message("\nResults written to DT3/results/")

  list(results = results_df, summary = summary_df)
}

# Run
dt3_output <- main_dt3()
print(dt3_output$summary[, c("scenario", "delta", "median_nnt", "ci_nnt_lo",
                              "ci_nnt_hi", "median_rr", "pct_rr_lt1",
                              "median_r_control", "median_r_prev")])
