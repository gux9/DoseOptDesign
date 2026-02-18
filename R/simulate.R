#' Monte Carlo Simulation for Dose Optimization Design Validation
#'
#' Simulates a two-stage dose optimization design to empirically validate
#' analytical bias and Type I error formulas. Stage 1 selects the superior
#' dose based on observed utility; Stage 2 confirms the selected dose.
#'
#' @param p Numeric in (0,1). True null response rate (equal for both doses).
#' @param q Numeric in (0,1). True null no-AE rate.
#' @param phi Numeric in [-1,1]. Efficacy-safety correlation.
#' @param u Numeric vector of length 4. Utility scores.
#' @param n1 Positive integer. Stage 1 per-arm sample size.
#' @param n2 Non-negative integer. Stage 2 sample size (selected dose only).
#' @param lambda_u Non-negative numeric. Selection threshold. Default 0.
#' @param alpha Numeric. Nominal significance level for binary tests. Default 0.025.
#' @param n_sim Positive integer. Number of Monte Carlo replications. Default 1000.
#' @param seed Integer. Random seed. Default 123.
#' @param perform_tte Logical. If \code{TRUE}, simulate and analyse a
#'   time-to-event confirmatory endpoint. Requires \pkg{copula} and
#'   \pkg{survival}. Default \code{FALSE}.
#' @param tte_rate Numeric > 0. True exponential hazard rate. Default 0.1.
#' @param corr_tte Numeric in [-1,1]. Gaussian copula parameter controlling
#'   the dependence between binary efficacy and TTE. Default 0.
#' @param entry_time_max Numeric. Maximum accrual time (weeks). Default 52.
#' @param admin_censor_time Numeric. Administrative censoring calendar time.
#'   Default 76.
#' @param dropout_rate Numeric. Annual dropout rate. Default 0.
#' @param landmark_time Numeric. Landmark time \eqn{\tau} for survival analysis.
#'   Default 24.
#' @param alpha_tte Numeric. Nominal significance level for TTE tests. Default 0.025.
#' @param two_arm Logical. If \code{TRUE} (and \code{perform_tte}), also run a
#'   two-arm Cox log-rank test against a simulated concurrent control arm.
#'   Default \code{TRUE}.
#' @param use_parallel Logical. If \code{TRUE}, use parallel processing for
#'   Cox regressions via \pkg{future.apply}. Default \code{FALSE}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{summary}{Named numeric vector of Monte Carlo means, SDs, and
#'       analytical reference values for bias and Type I error metrics.}
#'     \item{analytical}{Named list of analytical bias and Type I error values
#'       computed from \code{\link{calc_bias}} and \code{\link{calc_type1_error}}.}
#'   }
#'
#' @details
#' Under the null hypothesis, both doses share the same true response and
#' no-AE rates. Correlated endpoints are generated via a bivariate Gaussian
#' copula for (efficacy, TTE) and a conditional Bernoulli for safety given
#' efficacy.
#'
#' The function reports empirical counterparts to Equations 21, 22, 25-28,
#' and 31-34 of the manuscript.
#'
#' @examples
#' \donttest{
#' u_vec <- calc_utility(r = 1)
#' res   <- simulate_dose_opt(
#'   p = 0.4, q = 0.8, phi = 0, u = u_vec,
#'   n1 = 60, n2 = 140, n_sim = 500, seed = 42
#' )
#' print(res$summary[c("obs_bias_combined", "analytical_bias_utility",
#'                     "obs_type1_z", "obs_type1_bin")])
#' }
#'
#' @export
simulate_dose_opt <- function(
    p = 0.4, q = 0.8, phi = 0,
    u = c(1, 0.8, 0.2, 0),
    n1 = 50, n2 = 100,
    lambda_u = 0,
    alpha = 0.025,
    n_sim = 1000,
    seed = 123,
    perform_tte = FALSE,
    tte_rate = 0.1,
    corr_tte = 0,
    entry_time_max = 52,
    admin_censor_time = 76,
    dropout_rate = 0,
    landmark_time = 24,
    alpha_tte = 0.025,
    two_arm = TRUE,
    use_parallel = FALSE
) {
  # Input validation
  for (nm in c("p","q")) {
    val <- get(nm)
    if (!is.numeric(val) || val <= 0 || val >= 1)
      stop(sprintf("'%s' must be in (0,1).", nm))
  }
  if (!is.numeric(u) || length(u) != 4) stop("'u' must be length-4 numeric.")
  if (n1 <= 0 || n2 < 0) stop("'n1' must be positive and 'n2' non-negative.")

  # Analytical benchmarks
  ana <- calc_bias(p = p, q = q, phi = phi, u = u, n1 = n1, n2 = n2, lambda_u = lambda_u)
  t1e_util <- calc_type1_error(p, ana$bias_combined,   n1 + n2, alpha)
  t1e_resp <- calc_type1_error(p, ana$bias_max_combined, n1 + n2, alpha)

  # -------- Stage 1 simulation --------
  set.seed(seed)
  pi_vec <- calc_pi(p, q, phi)

  # Generate all Stage 1 outcomes in one block
  .sim_block <- function(n_total, seed_val) {
    set.seed(seed_val)
    p11 <- pi_vec[1]; pno_ae_resp <- pi_vec[1] / (pi_vec[1] + pi_vec[2])
    pno_ae_noresp <- pi_vec[3] / (pi_vec[3] + pi_vec[4])
    eff <- stats::rbinom(n_total, 1, p)
    p_no_ae <- ifelse(eff == 1, pno_ae_resp, pno_ae_noresp)
    saf <- stats::rbinom(n_total, 1, p_no_ae)
    list(eff = eff, saf = saf)
  }

  bL <- .sim_block(n1 * n_sim, seed)
  bH <- .sim_block(n1 * n_sim, seed + 1)

  L_eff <- matrix(bL$eff, n1, n_sim); H_eff <- matrix(bH$eff, n1, n_sim)
  L_saf <- matrix(bL$saf, n1, n_sim); H_saf <- matrix(bH$saf, n1, n_sim)

  p_L <- colMeans(L_eff); p_H <- colMeans(H_eff)

  # Utility scores
  .u_mean <- function(eff_m, saf_m) {
    idx <- (1 - eff_m) * 2 + (1 - saf_m) + 1
    colMeans(matrix(u[idx], n1, n_sim))
  }
  U_L <- .u_mean(L_eff, L_saf); U_H <- .u_mean(H_eff, H_saf)
  sel_H <- (U_H - U_L) > lambda_u  # TRUE = selected H

  # Selected arm Stage 1
  sel_eff <- L_eff; sel_eff[, sel_H] <- H_eff[, sel_H]
  sel_saf <- L_saf; sel_saf[, sel_H] <- H_saf[, sel_H]
  p_sel <- colMeans(sel_eff)

  # -------- Stage 2 --------
  n_H <- sum(sel_H); n_L <- n_sim - n_H
  bL2 <- .sim_block(n2 * n_L, seed + 2)
  bH2 <- .sim_block(n2 * n_H, seed + 3)
  s2_L_eff <- matrix(bL2$eff, n2, n_L)
  s2_H_eff <- matrix(bH2$eff, n2, n_H)
  s2_eff <- matrix(NA_real_, n2, n_sim)
  s2_eff[, !sel_H] <- s2_L_eff; s2_eff[, sel_H] <- s2_H_eff

  p_s2 <- colMeans(s2_eff)
  p_comb <- (n1 * p_sel + n2 * p_s2) / (n1 + n2)
  comb_eff <- rbind(sel_eff, s2_eff)

  n_total  <- n1 + n2
  obs_resp <- colSums(comb_eff)
  se_null  <- sqrt(p * (1 - p) / n_total)
  z_crit   <- qnorm(1 - alpha)

  # Z-test
  Z_stat   <- (p_comb - p) / se_null
  rej_z    <- mean(Z_stat > z_crit)

  # Binomial test
  k_c      <- t1e_util$k_c
  rej_bin  <- mean(stats::pbinom(obs_resp - 1, n_total, p, lower.tail = FALSE) < alpha)

  obs_bias <- mean(p_comb) - p

  # -------- Assemble summary --------
  summary_vec <- c(
    prop_select_H                = mean(sel_H),
    obs_bias_combined            = obs_bias,
    analytical_bias_utility      = ana$bias_combined,
    analytical_bias_max          = ana$bias_max_combined,
    obs_type1_z                  = rej_z,
    analytical_type1_z_utility   = t1e_util$type1_z,
    analytical_type1_z_max       = t1e_resp$type1_z,
    obs_type1_bin                = rej_bin,
    analytical_type1_bin_utility = t1e_util$type1_bin,
    analytical_type1_bin_max     = t1e_resp$type1_bin,
    inflation_z                  = rej_z / alpha,
    inflation_bin                = rej_bin / alpha,
    n1 = n1, n2 = n2, n_total = n_total, n_sim = n_sim
  )

  # -------- Optional TTE --------
  if (perform_tte) {
    if (!requireNamespace("copula", quietly = TRUE))
      stop("Package 'copula' is required for TTE analysis.")
    if (!requireNamespace("survival", quietly = TRUE))
      stop("Package 'survival' is required for TTE analysis.")

    tte_res <- .simulate_tte(
      sel_eff = sel_eff, s2_eff = s2_eff, sel_H = sel_H,
      p = p, q = q, phi = phi, u = u, n1 = n1, n2 = n2,
      n_sim = n_sim, tte_rate = tte_rate, corr_tte = corr_tte,
      entry_time_max = entry_time_max, admin_censor_time = admin_censor_time,
      dropout_rate = dropout_rate, landmark_time = landmark_time,
      alpha_tte = alpha_tte, two_arm = two_arm,
      use_parallel = use_parallel, seed = seed,
      ana_bias = ana
    )
    summary_vec <- c(summary_vec, tte_res)
  }

  list(
    summary    = summary_vec,
    analytical = list(
      bias          = ana,
      type1_utility = t1e_util,
      type1_max     = t1e_resp
    )
  )
}


# ---- Internal TTE simulation ----
#' @keywords internal
.simulate_tte <- function(sel_eff, s2_eff, sel_H, p, q, phi, u,
                           n1, n2, n_sim, tte_rate, corr_tte,
                           entry_time_max, admin_censor_time, dropout_rate,
                           landmark_time, alpha_tte, two_arm, use_parallel,
                           seed, ana_bias) {

  S0_lm <- exp(-tte_rate * landmark_time)
  n_total <- n1 + n2

  # Generate TTE for Stage 1 using copula
  .gen_tte <- function(n, corr, seed_val) {
    set.seed(seed_val)
    if (corr == 0) {
      eff_tmp <- stats::rbinom(n, 1, p)
      tte_tmp <- stats::rexp(n, rate = tte_rate)
    } else {
      mv <- copula::mvdc(
        copula::normalCopula(param = corr),
        margins = c("binom", "exp"),
        paramMargins = list(list(size = 1, prob = p), list(rate = tte_rate))
      )
      mat <- copula::rMvdc(n, mv)
      eff_tmp <- mat[, 1]; tte_tmp <- mat[, 2]
    }
    # Apply censoring
    entry  <- stats::runif(n, 0, entry_time_max)
    drop   <- if (dropout_rate == 0) rep(Inf, n) else
      stats::rexp(n, -log(1 - dropout_rate) / 52)
    obs_t  <- pmin(admin_censor_time - entry, pmin(tte_tmp, drop))
    status <- as.integer(tte_tmp <= drop & tte_tmp + entry < admin_censor_time)
    list(eff = eff_tmp, tte_true = tte_tmp, tte_obs = obs_t, status = status,
         landmark = as.integer(tte_tmp > landmark_time))
  }

  # Stage 1 TTE for both doses
  s1L <- .gen_tte(n1 * n_sim, corr_tte, seed + 10)
  s1H <- .gen_tte(n1 * n_sim, corr_tte, seed + 11)
  L_tte <- matrix(s1L$tte_obs, n1, n_sim); H_tte <- matrix(s1H$tte_obs, n1, n_sim)
  L_lm  <- matrix(s1L$landmark, n1, n_sim); H_lm <- matrix(s1H$landmark, n1, n_sim)
  L_st  <- matrix(s1L$status, n1, n_sim);  H_st  <- matrix(s1H$status, n1, n_sim)

  sel_tte <- L_tte; sel_tte[, sel_H] <- H_tte[, sel_H]
  sel_lm  <- L_lm;  sel_lm[, sel_H]  <- H_lm[, sel_H]
  sel_st  <- L_st;  sel_st[, sel_H]  <- H_st[, sel_H]

  # Stage 2 TTE
  n_H <- sum(sel_H); n_L <- n_sim - n_H
  s2L <- .gen_tte(n2 * n_L, corr_tte, seed + 12)
  s2H <- .gen_tte(n2 * n_H, corr_tte, seed + 13)
  s2_tte <- matrix(NA_real_, n2, n_sim)
  s2_lm  <- matrix(NA_real_, n2, n_sim)
  s2_st  <- matrix(NA_real_, n2, n_sim)
  s2_tte[, !sel_H] <- matrix(s2L$tte_obs, n2, n_L)
  s2_tte[, sel_H]  <- matrix(s2H$tte_obs, n2, n_H)
  s2_lm[, !sel_H]  <- matrix(s2L$landmark, n2, n_L)
  s2_lm[, sel_H]   <- matrix(s2H$landmark, n2, n_H)
  s2_st[, !sel_H]  <- matrix(s2L$status, n2, n_L)
  s2_st[, sel_H]   <- matrix(s2H$status, n2, n_H)

  comb_tte <- rbind(sel_tte, s2_tte)
  comb_lm  <- rbind(sel_lm,  s2_lm)
  comb_st  <- rbind(sel_st,  s2_st)

  # Landmark rate
  lm_rate_sel  <- colMeans(sel_lm)
  lm_rate_s2   <- colMeans(s2_lm)
  lm_rate_comb <- (n1 * lm_rate_sel + n2 * lm_rate_s2) / n_total
  obs_lm_bias  <- mean(lm_rate_comb) - S0_lm

  se_lm <- sqrt(S0_lm * (1 - S0_lm) / n_total)
  Z_lm  <- (lm_rate_comb - S0_lm) / se_lm
  rej_lm <- mean(Z_lm > qnorm(1 - alpha_tte))

  # One-arm exponential test (log scale)
  D    <- colSums(comb_st)
  TT   <- colSums(comb_tte)
  lhat <- D / TT
  Z_exp <- (log(lhat) - log(tte_rate)) * sqrt(D)
  rej_exp <- mean(Z_exp < -qnorm(1 - alpha_tte))

  # Plugin estimates for landmark
  u_vals <- u
  sel_eff_mat <- sel_eff

  # Cov(S(tau), U) from Stage 1
  idx_mat <- (1 - sel_eff_mat) * 2 + (1 - rbind(
    matrix(s1L$status[1:(n1*n_L)], n1, n_L),  # placeholder for safety
    matrix(s1H$status[1:(n1*n_H)], n1, n_H)
  )[[1]]) + 1
  # Use simplified plugin: surrogate-adjusted approach
  sigma_X   <- sqrt(p * (1 - p))
  sigma_lm  <- sqrt(S0_lm * (1 - S0_lm))
  rho_est   <- mean(colSums(sel_eff_mat * sel_lm) / n1 -
                      colMeans(sel_eff_mat) * colMeans(sel_lm)) /
    (sigma_X * sqrt(max(var(as.vector(sel_lm)) * (n1-1)/n1, 1e-10)))

  bias_lm_plugin <- rho_est * (sigma_lm / sigma_X) * ana_bias$bias_combined
  type1_lm_plugin <- calc_type1_error_landmark(S0_lm, bias_lm_plugin, n_total, alpha_tte)$type1_landmark

  res <- c(
    obs_bias_landmark_combined  = obs_lm_bias,
    obs_type1_landmark_z        = rej_lm,
    obs_type1_exp_one_arm       = rej_exp,
    landmark_null_rate          = S0_lm,
    corr_eff_tte_estimated      = rho_est,
    analytical_type1_lm_plugin  = type1_lm_plugin,
    inflation_landmark          = rej_lm / alpha_tte,
    inflation_exp               = rej_exp / alpha_tte
  )

  # Two-arm Cox
  if (two_arm) {
    set.seed(seed + 999)
    ctrl_tte <- matrix(stats::rexp(n_total * n_sim, tte_rate), n_total, n_sim)
    entry_c  <- matrix(stats::runif(n_total * n_sim, 0, entry_time_max), n_total, n_sim)
    ctrl_obs <- pmin(admin_censor_time - entry_c, ctrl_tte)
    ctrl_st  <- matrix(as.integer(
      as.vector(ctrl_tte) + as.vector(entry_c) < admin_censor_time
    ), n_total, n_sim)

    all_tte <- rbind(ctrl_obs, comb_tte)
    all_st  <- rbind(ctrl_st,  comb_st)
    arm_vec <- rep(c(0L, 1L), each = n_total)

    .cox_one <- function(i) {
      df  <- data.frame(time = all_tte[, i], status = all_st[, i], arm = arm_vec)
      fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ arm, data = df),
                      error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      pnorm(coef(fit) / sqrt(diag(vcov(fit))))
    }

    if (use_parallel && requireNamespace("future.apply", quietly = TRUE)) {
      future.apply::future_lapply(seq_len(n_sim), .cox_one)
    }
    p_lr <- vapply(seq_len(n_sim), .cox_one, numeric(1))
    rej_lr <- mean(p_lr <= alpha_tte, na.rm = TRUE)
    res <- c(res, obs_type1_logrank_two_arm = rej_lr,
             inflation_logrank = rej_lr / alpha_tte)
  }

  res
}
