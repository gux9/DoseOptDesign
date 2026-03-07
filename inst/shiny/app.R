

# Shiny App for Utility Score and ROSE Design Sample Size Calculation
# Enhanced with Bias & Type I Error Analysis
# Version: 2.0

library(shiny) 
library(shinythemes)
library(DT)
library(shinyBS)

#setwd('C://Users//xuemin.gu//Programming//BT8009//Dose.optimization//R.Code.Final//shiny')
#source('compare_bias_methods_fast_enhanced_v2.R')
#source('app_sample_size_functions_10Feb2026.R')


# ============================================================================
# ENHANCED compare_bias_methods_fast() FUNCTION - VERSION 2
# Two-Stage Dose Optimization with Two-Arm Dose Selection
# FIXES: Binomial test critical value, Added TTE bias/Type I error estimates
# ============================================================================

#' Enhanced fast comparison of bias methods with two-arm dose optimization
#' VERSION 2: Fixed binomial test + Added TTE test statistics from manuscript
#'
#' @description
#' Performs rapid simulation of a two-stage dose optimization design where:
#' - Stage 1: Two doses (L and H) are evaluated, best dose is selected
#' - Stage 2: Confirmatory trial with selected dose
#' 
#' Key improvements in v2:
#' 1. FIXED: Binomial test Type I error calculation (off-by-one error)
#' 2. ADDED: Landmark survival rate Z-test (Equation 31 in manuscript)
#' 3. ADDED: TTE bias plugin estimates (Equation 28)
#' 4. ADDED: Surrogate-adjusted TTE bias estimates
#'
#' @export
compare_bias_methods_fast_enhanced_v2 <- function(
    pL = 0.4, pH = 0.4,
    qL = 0.8, qH = 0.8,
    phi = 0.3,
    N1 = 50, N2 = 100,
    lambda_u = 0,
    u = c(100, 80, 20, 0),
    alpha_L = 0.8, alpha_H = 0.8,
    perform_tte_analysis = FALSE,
    tte_rate = 0.1,
    corr_efficacy_tte = 0,
    entry_time_max = 52,
    admin_censor_time = 76,
    dropout_rate = 0,
    unit_dropout = 12,
    tte_rate_historical = 0.1,
    two_arm_survival_flag = TRUE,
    landmark_time = 24,           # NEW: Landmark time for survival rate (weeks)
    nSim = 1000,
    ranseed = 123,
    Alpha = 0.025,
    Alpha_tte = 0.025,
    return_raw = FALSE,
    use_internal_parallel = TRUE   # NEW: Added 15Feb2026 (TURE for Standalone, FALSE for parallel batch run)
) {
  
  # ========================================================================
  # INPUT VALIDATION
  # ========================================================================
  
  if (any(c(pL, pH, qL, qH) < 0) || any(c(pL, pH, qL, qH) > 1)) {
    stop("Response and no-AE rates must be between 0 and 1")
  }
  if (phi < -1 || phi > 1) {
    stop("phi must be between -1 and 1")
  }
  if (N1 <= 0 || N2 <= 0) {
    stop("N1 and N2 must be positive")
  }
  if (length(u) != 4) {
    stop("Utility vector u must have 4 elements")
  }
  if (!all(c(Alpha, Alpha_tte) > 0 & c(Alpha, Alpha_tte) < 1)) {
    stop("Alpha and Alpha_tte must be between 0 and 1")
  }
  if (lambda_u < 0) {
    stop("lambda_u must be non-negative")
  }
  if (nSim <= 0) {
    stop("nSim must be positive")
  }
  
  # ========================================================================
  # LOAD REQUIRED LIBRARIES
  # ========================================================================
  
  require(copula, quietly = TRUE)
  
  if (perform_tte_analysis) {
    require(survival, quietly = TRUE)
    if (use_internal_parallel) {
      require(future.apply, quietly = TRUE)
      #require(progressr, quietly = TRUE)
    }
  }
  
  # ========================================================================
  # HELPER FUNCTIONS
  # ========================================================================
  # correlation between binary efficacy and safety
  calc_phi_hat <- function(eff_mat, safety_mat) {
    p  <- colMeans(eff_mat)
    q  <- colMeans(safety_mat)
    p11 <- colMeans(eff_mat * safety_mat)
    
    num <- p11 - p * q
    variance_term <- pmax(p * (1 - p) * q * (1 - q), 0)
    den <- sqrt(variance_term)
    
    phi <- ifelse(den < sqrt(.Machine$double.eps), 0, num / den)
    phi
  }
  
  # Option B: compute point-biserial correctly for binary efficacy and continues TTE (15Feb2026)
  calc_tte_est_fn <- function(eff_mat, tte_mat) {
    p <- colMeans(eff_mat)
    mean_tte <- colMeans(tte_mat)
    mean_tte_sq <- colMeans(tte_mat^2)
    var_tte <- mean_tte_sq - mean_tte^2
    cov_xt <- colMeans(eff_mat * tte_mat) - p * mean_tte
    cov_xt / sqrt(p * (1 - p) * pmax(var_tte, 1e-10))
  }
  
  # Multinomial probabilities
  calc_pi <- function(p, q, phi_param) {
    variance_term <- pmax(p * (1-p) * q * (1-q), 0)
    pi11 <- p * q + phi_param * sqrt(variance_term)
    pi10 <- p - pi11
    pi01 <- q - pi11
    pi00 <- 1 - p - q + pi11
    t(cbind(pi11, pi10, pi01, pi00))
  }
  
  # Utility moments
  calc_moments <- function(pi_vec, u_vec) {
    if(is.matrix(pi_vec)){
      if (length(u_vec) != nrow(pi_vec)) {
        stop("Length of u_vec must match the number of rows in pi_mat.")
      }
      
      mu <- crossprod(u_vec, pi_vec)
      u_vec_squared <- u_vec^2
      e_u_squared <- crossprod(u_vec_squared, pi_vec)
      sigma2 <- e_u_squared - mu^2
      
      mu_vec <- as.numeric(mu)
      sigma2_safe <- pmax(as.numeric(sigma2), 1e-10)
      
      return(list(mu = mu_vec, sigma2 = sigma2_safe))      
    } else {
      mu <- sum(u_vec * pi_vec)
      sigma2 <- sum(u_vec^2 * pi_vec) - mu^2
      list(mu = mu, sigma2 = max(sigma2, 1e-10))
    }
  }
  
  # Function to simulate correlated endpoints (updated to previous correct order on 15-Feb-2026)
  simulate_correlated_endpoints <- function(
    n, p, q, phi_param,
    tte_rate_val = NULL,
    corr_tte = 0,
    entry_time_max = NULL,
    admin_censor_time = NULL,
    dropout_rate_val = 0,
    unit_dropout_val = 12,
    seed_val = NULL,
    include_tte = FALSE
  ) {
    
    if (!is.null(seed_val)) set.seed(seed_val)
    
    # Joint probabilities for (efficacy, safety)
    pi_v <- calc_pi(p, q, phi_param)
    
    if (any(c(pi_v) < 0) || abs(sum(pi_v) - 1) > 1e-6) {
      stop("Invalid joint probabilities: must be non-negative and sum to 1")
    }
    
    p_no_ae_given_response    <- pi_v[1] / (pi_v[1] + pi_v[2])
    p_no_ae_given_no_response <- pi_v[3] / (pi_v[3] + pi_v[4])
    
    # Step 1: Generate efficacy (and TTE if needed)
    if (include_tte && corr_tte != 0) {
      # Joint (efficacy, TTE) via Gaussian copula
      mv_dist <- mvdc(
        copula = normalCopula(param = corr_tte),
        margins = c("binom", "exp"),
        paramMargins = list(
          list(size = 1, prob = p),
          list(rate = tte_rate_val)
        )
      )
      eff_tte <- rMvdc(n, mv_dist)
      efficacy <- eff_tte[, 1]
      tte      <- eff_tte[, 2]
    } else {
      efficacy <- rbinom(n, 1, p)
      if (include_tte) tte <- rexp(n, rate = tte_rate_val)
    }
    
    # Step 2: Generate safety conditional on the SAME efficacy vector
    probabilities_for_no_ae <- ifelse(efficacy == 1,
                                      p_no_ae_given_response,
                                      p_no_ae_given_no_response)
    no_ae <- rbinom(n, 1, probabilities_for_no_ae)
    
    # Return binary-only if TTE not requested
    if (!include_tte) {
      return(data.frame(
        efficacy = efficacy,
        no_ae = no_ae
      ))
    }
    
    # Censoring mechanism
    entry_time <- runif(n, 0, entry_time_max)
    dropout_time <- if (dropout_rate_val == 0) {
      rep(Inf, n)
    } else {
      rexp(n, -log(1 - dropout_rate_val) / unit_dropout_val)
    }
    
    observed_time <- pmin(admin_censor_time - entry_time, pmin(tte, dropout_time))
    status <- as.numeric(tte <= dropout_time & tte + entry_time < admin_censor_time)
    
    return(data.frame(
      efficacy = efficacy,
      no_ae = no_ae,
      tte_observed = observed_time,
      tte_true = tte,
      status = status
    ))
  }
  
  # Simulate PFS only for historical control
  simulate_pfs_only <- function(n, tte_rate_val, entry_time_max, admin_censor_time,
                                dropout_rate_val, unit_dropout_val, seed_val) {
    set.seed(seed_val)
    entry_time <- runif(n, 0, entry_time_max)
    tte <- rexp(n, rate = tte_rate_val)
    tte[is.infinite(tte) | is.na(tte)] <- admin_censor_time * 10
    
    dropout_time <- if (dropout_rate_val == 0) {
      rep(Inf, n)
    } else {
      rexp(n, -log(1 - dropout_rate_val) / unit_dropout_val)
    }
    dropout_time[is.infinite(dropout_time) | is.na(dropout_time)] <- Inf
    
    observed_time <- pmin(admin_censor_time - entry_time, pmin(tte, dropout_time))
    observed_time[observed_time <= 0] <- 1e-6
    
    status <- as.integer(tte <= dropout_time & tte + entry_time < admin_censor_time)
    status[is.na(status)] <- 0
    
    data.frame(tte_observed = as.numeric(observed_time), status = as.integer(status))
  }
  
  # Compute utility score
  compute_utility <- function(efficacy_mat, no_ae_mat, uti_score) {
    idx <- (1 - efficacy_mat) * 2 + (1 - no_ae_mat) + 1
    utility_mat <- matrix(uti_score[idx], nrow = nrow(efficacy_mat))
    colMeans(utility_mat)
  }
  
  # ========================================================================
  # ANALYTICAL BIAS CALCULATIONS (NULL CASE)
  # ========================================================================
  
  pi_L <- calc_pi(pL, qL, phi)
  pi_H <- calc_pi(pH, qH, phi)
  
  mom_L <- calc_moments(pi_L, u)
  mom_H <- calc_moments(pi_H, u)
  
  E_U_L <- mom_L$mu
  E_U_H <- mom_H$mu
  
  avg_sigma2 <- (mom_L$sigma2 + mom_H$sigma2) / 2
  
  x_vec <- c(1, 1, 0, 0)
  cov_L <- sum(x_vec * u * pi_L) - pL * E_U_L
  cov_H <- sum(x_vec * u * pi_H) - pH * E_U_H
  avg_cov <- (cov_L + cov_H) / 2
  
  exp_term <- exp(-(lambda_u^2 * N1) / (4 * avg_sigma2))
  
  bias_utility_analytical_stage1 <- (avg_cov / sqrt(avg_sigma2 * N1)) * (1 / sqrt(base::pi)) * exp_term
  bias_utility_analytical_combined <- (N1 / (N1 + N2)) * bias_utility_analytical_stage1
  
  bias_response_analytical_stage1 <- sqrt(pL * (1 - pL)) / sqrt(N1 * base::pi)
  bias_response_analytical_combined <- (N1 / (N1 + N2)) * bias_response_analytical_stage1
  
  # Analytical utility score bias (for U itself): Cov(U,U)/σ_U = σ_U
  bias_utility_score_analytical_stage1 <- sqrt(avg_sigma2) / sqrt(N1 * base::pi) * exp_term
  bias_utility_score_analytical_combined <- (N1 / (N1 + N2)) * bias_utility_score_analytical_stage1
  
  # ========================================================================
  # Type I error calculations (null case) - ANALYTICAL
  # ========================================================================
  p_null <- pL
  n_total <- N1 + N2
  z_crit <- qnorm(1 - Alpha)
  
  if (p_null > 0 && p_null < 1) {
    se_null <- sqrt(p_null * (1 - p_null) / n_total)
  } else {
    se_null <- 1e-10
  }  
  
  # ========================================================================
  # FIXED: Critical value calculation for binomial test
  # Using clearer indexing to avoid off-by-one errors
  # ========================================================================
  
  # Find critical value k_c such that P(X >= k_c | H0) <= Alpha
  # This is the rejection region: reject H0 when X >= k_c
  
  probs_under_null <- pbinom(0:n_total, size = n_total, prob = p_null, lower.tail = FALSE)
  # probs_under_null[k+1] = P(X > k) = P(X >= k+1)
  
  # We want smallest k such that P(X >= k) <= Alpha
  # P(X >= k) = P(X > k-1) = pbinom(k-1, n, p, lower.tail=FALSE)
  # = probs_under_null[k] (since probs_under_null is 0-indexed conceptually)
  
  # Actually in R, probs_under_null[1] = P(X > 0) = P(X >= 1)
  # probs_under_null[k] = P(X > k-1) = P(X >= k)
  
  # Find minimum k (1-indexed) such that probs_under_null[k] <= Alpha
  k_c_index <- min(which(probs_under_null <= Alpha))
  k_c <- k_c_index  # This is the critical value: reject when X >= k_c
  
  # Verify: P(X >= k_c) = probs_under_null[k_c] should be <= Alpha
  
  # Analytical Type I error estimates
  type1_z_utility_analytical_combined <- 1 - pnorm(z_crit - bias_utility_analytical_combined / se_null)
  type1_z_response_analytical_combined <- 1 - pnorm(z_crit - bias_response_analytical_combined / se_null)
  
  # FIXED: Binomial test Type I error
  # Type I error = P(X >= k_c | p = p_null + bias)
  type1_bin_utility_analytical_combined <- pbinom(k_c - 1, size = n_total, 
                                                  prob = p_null + bias_utility_analytical_combined, 
                                                  lower.tail = FALSE)
  type1_bin_response_analytical_combined <- pbinom(k_c - 1, size = n_total, 
                                                   prob = p_null + bias_response_analytical_combined, 
                                                   lower.tail = FALSE)
  
  # ========================================================================
  # INITIALIZE SIMULATION RESULTS DATA FRAME
  # ========================================================================
  
  base_columns <- list(
    sim_id        = 1:nSim,
    selected_dose = NA_character_,
    
    # Stage 1 estimates
    p_hat_L_stage1          = NA_real_,
    q_hat_L_stage1          = NA_real_,
    utility_hat_L_stage1    = NA_real_,
    p_hat_H_stage1          = NA_real_,
    q_hat_H_stage1          = NA_real_,
    utility_hat_H_stage1    = NA_real_,
    utility_diff_hat_stage1 = NA_real_,
    
    # Post-selection estimates
    p_true_selected = NA_real_,
    p_hat_selected  = NA_real_,
    q_hat_selected  = NA_real_,
    p_hat_stage2    = NA_real_,
    q_hat_stage2    = NA_real_,
    p_hat_combined  = NA_real_,
    q_hat_combined  = NA_real_,
    response_combined = NA_real_,
    no_ae_combined    = NA_real_,
    
    # Observed bias
    bias_observed_selected_combined = NA_real_,
    
    # Plugin bias estimates - Null case (selected dose method)
    bias_utility_plugin_selected_stage1   = NA_real_,
    bias_utility_plugin_selected_combined = NA_real_,
    bias_response_plugin_selected_stage1  = NA_real_,
    bias_response_plugin_selected_combined= NA_real_,
    
    # Plugin estimates for variance components
    phi_hat_plugin_selected_stage1    = NA_real_,
    cov_XU_hat_plugin_selected_stage1 = NA_real_,
    mean_U_hat_plugin_selected_stage1 = NA_real_,
    sigma_U_hat_plugin_selected_stage1= NA_real_,
    
    # Hypothesis testing - BINARY ENDPOINTS
    p_value_z_test   = NA_real_,
    reject_H0_z      = NA_real_,
    p_value_bin_test = NA_real_,
    reject_H0_bin    = NA_real_,
    
    # Plugin Type I error estimates - BINARY
    type1_z_utility_plugin  = NA_real_,
    type1_bin_utility_plugin= NA_real_,
    type1_z_response_plugin = NA_real_,
    type1_bin_response_plugin= NA_real_,
    
    # Correlation estimates
    phi_est = NA_real_
  )
  
  # TTE-specific columns
  if (perform_tte_analysis) {
    tte_columns <- list(
      # Basic TTE metrics
      events_stage1 = NA_real_,
      events_stage2 = NA_real_,
      median_pfs    = NA_real_,
      
      # Landmark survival rate (Section 2.7.5 of manuscript)
      landmark_survived_stage1  = NA_real_,
      landmark_survived_stage2  = NA_real_,
      landmark_survived_combined = NA_real_,
      landmark_rate_stage1      = NA_real_,
      landmark_rate_stage2      = NA_real_,
      landmark_rate_combined    = NA_real_,
      
      # Observed TTE bias
      bias_landmark_observed_combined = NA_real_,
      
      # Plugin TTE bias estimates (Equation 28)
      cov_landmark_U_plugin     = NA_real_,
      bias_landmark_plugin_stage1   = NA_real_,
      bias_landmark_plugin_combined = NA_real_,
      
      # Surrogate-adjusted bias (Method 2, Section 2.7.4)
      corr_response_landmark_plugin = NA_real_,
      bias_landmark_surrogate_adjusted_combined = NA_real_,
      
      # Upper bound bias (Method 1, Section 2.7.4)
      bias_landmark_upper_bound_combined = NA_real_,
      
      # Landmark survival Z-test (Equation 31)
      p_value_landmark_z_test = NA_real_,
      reject_H0_landmark_z    = NA_real_,
      
      # Plugin Type I error for landmark test
      type1_landmark_z_plugin = NA_real_,
      type1_landmark_z_upper_bound = NA_real_,
      type1_landmark_z_surrogate = NA_real_,
      
      # One-sample exponential test
      p_value_exponential = NA_real_,
      reject_H0_exp_one_arm = NA_real_,
      
      # Plugin Type I error for one-arm exponential test
      # Uses direct Cov(TTE, U) estimation from Stage 1 data
      type1_exp_one_arm_plugin = NA_real_,
      type1_exp_one_arm_upper = NA_real_,
      
      # Direct TTE bias estimation from Stage 1 (Cov(TTE, U))
      cov_TTE_U_plugin = NA_real_,
      corr_TTE_U_plugin = NA_real_,
      bias_TTE_plugin = NA_real_,
      
      # Two-arm survival (Cox regression / Log-rank)
      hr                        = NA_real_,
      p_value_wald_one_sided    = NA_real_,
      p_value_wald_two_sided    = NA_real_,
      p_value_logrank_two_sided = NA_real_,
      p_value_logrank_one_sided = NA_real_,
      reject_H0_logrank_two_arm = NA_real_,
      
      # Plugin Type I error for two-arm log-rank test
      # Uses Cov(S(tau), U) from Stage 1 data
      type1_logrank_two_arm_plugin = NA_real_,
      type1_logrank_two_arm_upper = NA_real_,
      
      # Correlation estimates
      corr_tte_est = NA_real_
    )
    sim_results <- as.data.frame(c(base_columns, tte_columns), stringsAsFactors = FALSE)
  } else {
    sim_results <- as.data.frame(base_columns, stringsAsFactors = FALSE)
  }
  
  # ========================================================================
  # SIMULATION: STAGE 1 - Generate data for both doses
  # ========================================================================
  
  data_L <- simulate_correlated_endpoints(
    n = N1*nSim, p = pL, q = qL, phi_param = phi, 
    tte_rate_val = tte_rate,
    corr_tte = corr_efficacy_tte, 
    entry_time_max = entry_time_max,
    admin_censor_time = admin_censor_time, 
    dropout_rate_val = dropout_rate,
    unit_dropout_val = unit_dropout, 
    seed_val = ranseed,
    include_tte = perform_tte_analysis
  )
  
  data_H <- simulate_correlated_endpoints(
    n = N1*nSim, p = pH, q = qH, phi_param = phi, 
    tte_rate_val = tte_rate,
    corr_tte = corr_efficacy_tte, 
    entry_time_max = entry_time_max,
    admin_censor_time = admin_censor_time, 
    dropout_rate_val = dropout_rate,
    unit_dropout_val = unit_dropout, 
    seed_val = ranseed + 1,
    include_tte = perform_tte_analysis
  )
  
  # Reshape into matrices
  L_efficacy_mat <- matrix(data_L$efficacy, nrow=N1, ncol=nSim)
  H_efficacy_mat <- matrix(data_H$efficacy, nrow=N1, ncol=nSim)
  L_no_ae_mat <- matrix(data_L$no_ae, nrow=N1, ncol=nSim)
  H_no_ae_mat <- matrix(data_H$no_ae, nrow=N1, ncol=nSim)
  
  if (perform_tte_analysis) {
    L_tte_mat <- matrix(data_L$tte_observed, nrow=N1, ncol=nSim)
    H_tte_mat <- matrix(data_H$tte_observed, nrow=N1, ncol=nSim)
    L_tte_true_mat <- matrix(data_L$tte_true, nrow=N1, ncol=nSim)
    H_tte_true_mat <- matrix(data_H$tte_true, nrow=N1, ncol=nSim)
    L_censor_mat <- matrix(data_L$status, nrow=N1, ncol=nSim)
    H_censor_mat <- matrix(data_H$status, nrow=N1, ncol=nSim)
    
    # Landmark survival indicators (using true TTE for simplicity)
    L_landmark_mat <- matrix(as.numeric(data_L$tte_true > landmark_time), nrow=N1, ncol=nSim)
    H_landmark_mat <- matrix(as.numeric(data_H$tte_true > landmark_time), nrow=N1, ncol=nSim)
  }
  
  # Stage 1 statistics
  p_hat_L <- colMeans(L_efficacy_mat, na.rm=TRUE)
  p_hat_H <- colMeans(H_efficacy_mat, na.rm=TRUE)
  q_hat_L <- colMeans(L_no_ae_mat, na.rm=TRUE)
  q_hat_H <- colMeans(H_no_ae_mat, na.rm=TRUE)
  
  # Utility score-based selection
  U_L <- compute_utility(L_efficacy_mat, L_no_ae_mat, u)
  U_H <- compute_utility(H_efficacy_mat, H_no_ae_mat, u)
  selected_dose <- ifelse((U_H - U_L) > lambda_u, "H", "L")
  
  # Count selections
  L_nSim2 <- sum(selected_dose == "L")
  H_nSim2 <- sum(selected_dose == "H")
  
  # ========================================================================
  # STAGE 1: Selected arm data
  # ========================================================================
  
  selected_data_efficacy <- L_efficacy_mat
  selected_data_efficacy[, selected_dose == "H"] <- H_efficacy_mat[, selected_dose == "H"]
  
  selected_data_safety <- L_no_ae_mat
  selected_data_safety[, selected_dose == "H"] <- H_no_ae_mat[, selected_dose == "H"]
  
  if (perform_tte_analysis) {
    selected_data_tte <- L_tte_mat
    selected_data_tte[, selected_dose == "H"] <- H_tte_mat[, selected_dose == "H"]
    
    selected_data_censor <- L_censor_mat
    selected_data_censor[, selected_dose == "H"] <- H_censor_mat[, selected_dose == "H"]
    
    selected_data_landmark <- L_landmark_mat
    selected_data_landmark[, selected_dose == "H"] <- H_landmark_mat[, selected_dose == "H"]
  }
  
  p_selected <- colMeans(selected_data_efficacy)
  q_selected <- colMeans(selected_data_safety)
  p_true_selected <- ifelse(selected_dose == "L", pL, pH)
  
  # ========================================================================
  # STAGE 2: Generate confirmatory data
  # ========================================================================
  
  data_stage2_L <- simulate_correlated_endpoints(
    n = N2*L_nSim2, p = pL, q = qL,
    phi_param = phi, tte_rate_val = tte_rate,
    corr_tte = corr_efficacy_tte, entry_time_max = entry_time_max,
    admin_censor_time = admin_censor_time, dropout_rate_val = dropout_rate,
    unit_dropout_val = unit_dropout, seed_val = ranseed + 2,
    include_tte = perform_tte_analysis
  )
  
  data_stage2_H <- simulate_correlated_endpoints(
    n = N2*H_nSim2, p = pH, q = qH,
    phi_param = phi, tte_rate_val = tte_rate,
    corr_tte = corr_efficacy_tte, entry_time_max = entry_time_max,
    admin_censor_time = admin_censor_time, dropout_rate_val = dropout_rate,
    unit_dropout_val = unit_dropout, seed_val = ranseed + 3,
    include_tte = perform_tte_analysis
  ) 
  
  # ========================================================================
  # STAGE 2: Reshape and align data
  # ========================================================================
  
  stage2_L_efficacy <- matrix(data_stage2_L$efficacy, nrow = N2, ncol = L_nSim2)
  stage2_H_efficacy <- matrix(data_stage2_H$efficacy, nrow = N2, ncol = H_nSim2)
  stage2_L_safety   <- matrix(data_stage2_L$no_ae,    nrow = N2, ncol = L_nSim2)
  stage2_H_safety   <- matrix(data_stage2_H$no_ae,    nrow = N2, ncol = H_nSim2)
  
  if (perform_tte_analysis) {
    stage2_L_tte      <- matrix(data_stage2_L$tte_observed, nrow = N2, ncol = L_nSim2)
    stage2_H_tte      <- matrix(data_stage2_H$tte_observed, nrow = N2, ncol = H_nSim2)
    stage2_L_tte_true <- matrix(data_stage2_L$tte_true, nrow = N2, ncol = L_nSim2)
    stage2_H_tte_true <- matrix(data_stage2_H$tte_true, nrow = N2, ncol = H_nSim2)
    stage2_L_censor   <- matrix(data_stage2_L$status,       nrow = N2, ncol = L_nSim2)
    stage2_H_censor   <- matrix(data_stage2_H$status,       nrow = N2, ncol = H_nSim2)
    stage2_L_landmark <- matrix(as.numeric(data_stage2_L$tte_true > landmark_time), nrow = N2, ncol = L_nSim2)
    stage2_H_landmark <- matrix(as.numeric(data_stage2_H$tte_true > landmark_time), nrow = N2, ncol = H_nSim2)
  }
  
  stage2_data_efficacy <- matrix(NA_real_, N2, nSim)
  stage2_data_efficacy[, selected_dose == "L"] <- stage2_L_efficacy
  stage2_data_efficacy[, selected_dose == "H"] <- stage2_H_efficacy
  
  stage2_data_safety <- matrix(NA_real_, N2, nSim)
  stage2_data_safety[, selected_dose == "L"] <- stage2_L_safety
  stage2_data_safety[, selected_dose == "H"] <- stage2_H_safety
  
  if (perform_tte_analysis) {
    stage2_data_tte <- matrix(NA_real_, N2, nSim)
    stage2_data_tte[, selected_dose == "L"] <- stage2_L_tte
    stage2_data_tte[, selected_dose == "H"] <- stage2_H_tte
    
    stage2_data_censor <- matrix(NA_real_, N2, nSim)
    stage2_data_censor[, selected_dose == "L"] <- stage2_L_censor
    stage2_data_censor[, selected_dose == "H"] <- stage2_H_censor
    
    stage2_data_landmark <- matrix(NA_real_, N2, nSim)
    stage2_data_landmark[, selected_dose == "L"] <- stage2_L_landmark
    stage2_data_landmark[, selected_dose == "H"] <- stage2_H_landmark
  }
  
  # Combined data
  combined_efficacy_mat <- rbind(selected_data_efficacy, stage2_data_efficacy)
  combined_safety_mat <- rbind(selected_data_safety, stage2_data_safety)
  
  if (perform_tte_analysis) {
    combined_tte_mat <- rbind(selected_data_tte, stage2_data_tte)
    combined_censor_mat <- rbind(selected_data_censor, stage2_data_censor)
    combined_landmark_mat <- rbind(selected_data_landmark, stage2_data_landmark)
  }
  
  # Combined estimates
  p_stage2 <- colMeans(stage2_data_efficacy, na.rm=TRUE)
  q_stage2 <- colMeans(stage2_data_safety, na.rm=TRUE)
  p_hat_combined <- (N1 * p_selected + N2 * p_stage2) / (N1 + N2)
  q_hat_combined <- (N1 * q_selected + N2 * q_stage2) / (N1 + N2)
  
  # ========================================================================
  # STORE BASIC RESULTS
  # ========================================================================
  
  sim_results$sim_id <- 1:nSim
  sim_results$selected_dose <- selected_dose
  sim_results$p_hat_L_stage1 <- p_hat_L
  sim_results$p_hat_H_stage1 <- p_hat_H
  sim_results$q_hat_L_stage1 <- q_hat_L
  sim_results$q_hat_H_stage1 <- q_hat_H
  sim_results$p_hat_selected <- p_selected
  sim_results$q_hat_selected <- q_selected
  sim_results$p_hat_stage2 <- p_stage2
  sim_results$q_hat_stage2 <- q_stage2
  sim_results$p_hat_combined <- p_hat_combined
  sim_results$q_hat_combined <- q_hat_combined
  sim_results$p_true_selected <- p_true_selected
  sim_results$utility_hat_L_stage1 <- U_L
  sim_results$utility_hat_H_stage1 <- U_H
  sim_results$utility_diff_hat_stage1 <- U_H - U_L
  
  # ========================================================================
  # OBSERVED BIAS - BINARY ENDPOINT
  # ========================================================================
  
  sim_results$bias_observed_selected_combined <- p_hat_combined - p_true_selected
  
  # Observed utility score bias
  utility_selected_stage1 <- compute_utility(selected_data_efficacy, selected_data_safety, u)
  utility_stage2 <- compute_utility(stage2_data_efficacy, stage2_data_safety, u)
  utility_selected_combined <- (N1 * utility_selected_stage1 + N2 * utility_stage2) / n_total
  mu_U_true <- E_U_L
  sim_results$bias_utility_score_observed_combined <- utility_selected_combined - mu_U_true
  
  # ========================================================================
  # PLUGIN BIAS ESTIMATES - BINARY ENDPOINT
  # ========================================================================
  
  phi_hat_plugin <- calc_phi_hat(selected_data_efficacy, selected_data_safety)
  
  pi_sel <- calc_pi(p_selected, q_selected, phi_hat_plugin)
  mom_sel <- calc_moments(pi_sel, u)
  
  cov_XU_hat <- colSums((x_vec * u) * pi_sel) - p_selected * mom_sel$mu
  
  exp_term_plugin <- exp(-(lambda_u^2 * N1) / (4 * mom_sel$sigma2))
  bias_util_plugin_stage1 <- (cov_XU_hat / sqrt(mom_sel$sigma2 * N1)) * (1 / sqrt(base::pi)) * exp_term_plugin
  bias_util_plugin_combined <- (N1 / (N1 + N2)) * bias_util_plugin_stage1
  
  sim_results$phi_hat_plugin_selected_stage1 <- phi_hat_plugin
  sim_results$cov_XU_hat_plugin_selected_stage1 <- cov_XU_hat
  sim_results$mean_U_hat_plugin_selected_stage1 <- mom_sel$mu
  sim_results$sigma_U_hat_plugin_selected_stage1 <- sqrt(mom_sel$sigma2)
  
  sim_results$bias_utility_plugin_selected_stage1 <- bias_util_plugin_stage1
  sim_results$bias_utility_plugin_selected_combined <- bias_util_plugin_combined
  
  # Response-only maximum bound
  bias_resp_plugin_stage1 <- sqrt(p_selected * (1 - p_selected)) / sqrt(N1 * base::pi)
  bias_resp_plugin_combined <- (N1 / (N1 + N2)) * bias_resp_plugin_stage1
  
  sim_results$bias_response_plugin_selected_stage1 <- bias_resp_plugin_stage1
  sim_results$bias_response_plugin_selected_combined <- bias_resp_plugin_combined
  
  # Utility score plugin bias (for U itself): σ̂_U/√(n₁π)
  bias_U_score_plugin_stage1 <- sqrt(mom_sel$sigma2) / sqrt(N1 * base::pi) * exp_term_plugin
  bias_U_score_plugin_combined <- (N1 / n_total) * bias_U_score_plugin_stage1
  sim_results$bias_utility_score_plugin_combined <- bias_U_score_plugin_combined
  
  # ========================================================================
  # HYPOTHESIS TESTING - BINARY ENDPOINTS
  # ========================================================================
  
  events_total <- colSums(combined_efficacy_mat)
  sim_results$response_combined <- events_total
  sim_results$no_ae_combined <- colSums(combined_safety_mat)
  
  # Z-test
  Z_stat <- (p_hat_combined - p_null) / se_null
  p_value_z <- pnorm(Z_stat, lower.tail = FALSE)
  
  sim_results$p_value_z_test <- p_value_z
  sim_results$reject_H0_z <- as.numeric(p_value_z < Alpha)
  
  # Binomial exact test
  # p-value = P(X >= observed | H0) = pbinom(observed - 1, n, p0, lower.tail = FALSE)
  p_value_bin <- pbinom(events_total - 1, n_total, p_null, lower.tail = FALSE)
  
  sim_results$p_value_bin_test <- p_value_bin
  sim_results$reject_H0_bin <- as.numeric(p_value_bin < Alpha)
  
  # ========================================================================
  # PLUGIN TYPE I ERROR ESTIMATES - BINARY (FIXED)
  # ========================================================================
  
  # Z-test: Type I error = P(Z > z_crit | H1) = 1 - Phi(z_crit - bias/SE)
  z_shift_util <- bias_util_plugin_combined / se_null
  sim_results$type1_z_utility_plugin <- 1 - pnorm(z_crit - z_shift_util)
  
  z_shift_resp <- bias_resp_plugin_combined / se_null
  sim_results$type1_z_response_plugin <- 1 - pnorm(z_crit - z_shift_resp)
  
  # FIXED: Binomial test Type I error = P(X >= k_c | p = p0 + bias)
  # = pbinom(k_c - 1, n, p0 + bias, lower.tail = FALSE)
  sim_results$type1_bin_utility_plugin <- pbinom(k_c - 1, size = n_total, 
                                                 prob = p_null + bias_util_plugin_combined, 
                                                 lower.tail = FALSE)
  sim_results$type1_bin_response_plugin <- pbinom(k_c - 1, size = n_total, 
                                                  prob = p_null + bias_resp_plugin_combined, 
                                                  lower.tail = FALSE)
  
  # ========================================================================
  # CORRELATION ESTIMATES
  # ========================================================================
  
  sim_results$phi_est <- calc_phi_hat(combined_efficacy_mat, combined_safety_mat)
  sim_results$phi_est[(is.na(sim_results$phi_est) | !is.finite(sim_results$phi_est))] <- 0
  
  # ========================================================================
  # TIME-TO-EVENT ANALYSES (IF REQUESTED)
  # ========================================================================
  
  if (perform_tte_analysis) {
    
    # ====================================================================
    # LANDMARK SURVIVAL RATE CALCULATIONS (NEW - Equation 31)
    # ====================================================================
    
    # Landmark survival rate: S(tau) = P(T > tau)
    # Under null, true landmark rate is exp(-tte_rate * landmark_time)
    S0_landmark <- exp(-tte_rate * landmark_time)
    
    # Observed landmark survival rates
    landmark_stage1 <- colSums(selected_data_landmark)
    landmark_stage2 <- colSums(stage2_data_landmark)
    landmark_combined <- landmark_stage1 + landmark_stage2
    
    sim_results$landmark_survived_stage1 <- landmark_stage1
    sim_results$landmark_survived_stage2 <- landmark_stage2
    sim_results$landmark_survived_combined <- landmark_combined
    
    landmark_rate_stage1 <- colMeans(selected_data_landmark)
    landmark_rate_stage2 <- colMeans(stage2_data_landmark)
    landmark_rate_combined <- (N1 * landmark_rate_stage1 + N2 * landmark_rate_stage2) / n_total
    
    sim_results$landmark_rate_stage1 <- landmark_rate_stage1
    sim_results$landmark_rate_stage2 <- landmark_rate_stage2
    sim_results$landmark_rate_combined <- landmark_rate_combined
    
    # Observed landmark bias
    sim_results$bias_landmark_observed_combined <- landmark_rate_combined - S0_landmark
    
    # ====================================================================
    # PLUGIN TTE BIAS ESTIMATES (Equation 28)
    # ====================================================================
    
    # S_i(tau) = 1{T_i > tau} is binary, similar structure to response X
    # Bias(S(tau)) = Cov(S_i(tau), U) / (sigma_U * sqrt(n1)) * (1/sqrt(pi)) * exp(...)
    
    # Compute Cov(S_i(tau), U) for selected dose
    # S_i(tau) takes values 0 or 1, so this is analogous to Cov(X, U)
    
    # ====================================================================
    # VECTORIZED COVARIANCE CALCULATIONS (preserving matrix-based speed)
    # ====================================================================
    
    # Compute correlation between landmark survival and response (surrogate)
    #corr_response_landmark <- calc_phi_hat(selected_data_efficacy, selected_data_landmark)
    corr_response_landmark <- calc_tte_est_fn(selected_data_efficacy, selected_data_landmark)
    corr_response_landmark[!is.finite(corr_response_landmark)] <- 0
    sim_results$corr_response_landmark_plugin <- corr_response_landmark
    
    # Variance components
    sigma_S_stage1 <- sqrt(landmark_rate_stage1 * (1 - landmark_rate_stage1))
    sigma_X_stage1 <- sqrt(p_selected * (1 - p_selected))
    sigma_U_stage1 <- sqrt(mom_sel$sigma2)
    
    # Correlation between X and U
    corr_XU <- cov_XU_hat / (sigma_X_stage1 * sigma_U_stage1)
    corr_XU[!is.finite(corr_XU)] <- 0
    
    # ====================================================================
    # VECTORIZED: Direct estimation of Cov(S(tau), U) from Stage 1 data
    # Cov(S, U) = E[S*U] - E[S]*E[U]
    # ====================================================================
    
    # Compute utility scores for all patients in Stage 1 (matrix form)
    # idx = (1 - efficacy) * 2 + (1 - safety) + 1 maps to utility index
    utility_idx_stage1 <- (1 - selected_data_efficacy) * 2 + (1 - selected_data_safety) + 1
    utility_mat_stage1 <- matrix(u[utility_idx_stage1], nrow = N1, ncol = nSim)
    
    # E[S*U] = mean(S * U) for each simulation (column)
    E_SU <- colMeans(selected_data_landmark * utility_mat_stage1)
    
    # E[S] = landmark_rate_stage1 (already computed)
    # E[U] = mean utility for selected dose Stage 1
    E_U_stage1 <- colMeans(utility_mat_stage1)
    
    # Cov(S, U) = E[SU] - E[S]*E[U]
    cov_SU_direct <- E_SU - landmark_rate_stage1 * E_U_stage1
    cov_SU_direct[!is.finite(cov_SU_direct)] <- 0
    sim_results$cov_landmark_U_plugin <- cov_SU_direct
    
    # Plugin bias for landmark using direct Cov estimate (Equation 28)
    bias_landmark_plugin_stage1 <- (cov_SU_direct / sqrt(mom_sel$sigma2 * N1)) * 
      (1 / sqrt(base::pi)) * exp_term_plugin
    bias_landmark_plugin_combined <- (N1 / n_total) * bias_landmark_plugin_stage1
    
    sim_results$bias_landmark_plugin_stage1 <- bias_landmark_plugin_stage1
    sim_results$bias_landmark_plugin_combined <- bias_landmark_plugin_combined
    
    # Method 2: Surrogate-adjusted estimate
    # Uses utility-based response bias plugin (Cov(X,U)/σ_U path)
    bias_landmark_surrogate <- corr_response_landmark * (sigma_S_stage1 / sigma_X_stage1) * 
      bias_util_plugin_combined
    bias_landmark_surrogate[!is.finite(bias_landmark_surrogate)] <- 0
    sim_results$bias_landmark_surrogate_adjusted_combined <- bias_landmark_surrogate
    
    # Response-based upper bound (assume Corr(S, X) = 1)
    bias_landmark_resp_upper <- (sigma_S_stage1 / sigma_X_stage1) * bias_resp_plugin_combined
    bias_landmark_resp_upper[!is.finite(bias_landmark_resp_upper)] <- 0
    sim_results$bias_landmark_upper_bound_combined <- bias_landmark_resp_upper
    
    # Utility-score-based upper bound (assume Corr(S, U) = σ_S/σ_U)
    bias_landmark_util_upper <- (sigma_S_stage1 / sigma_U_stage1) * bias_U_score_plugin_combined
    bias_landmark_util_upper[!is.finite(bias_landmark_util_upper)] <- 0
    sim_results$bias_landmark_utility_upper_combined <- bias_landmark_util_upper
    
    # ====================================================================
    # LANDMARK SURVIVAL Z-TEST (NEW - Equation 31)
    # ====================================================================
    
    # Test H0: S(tau) <= S0 vs H1: S(tau) > S0
    se_S0 <- sqrt(S0_landmark * (1 - S0_landmark) / n_total)
    
    Z_landmark <- (landmark_rate_combined - S0_landmark) / se_S0
    p_value_landmark_z <- pnorm(Z_landmark, lower.tail = FALSE)
    
    sim_results$p_value_landmark_z_test <- p_value_landmark_z
    sim_results$reject_H0_landmark_z <- as.numeric(p_value_landmark_z < Alpha_tte)
    
    # Plugin Type I error for landmark Z-test
    z_crit_tte <- qnorm(1 - Alpha_tte)
    
    # Using direct plugin bias
    z_shift_landmark <- bias_landmark_plugin_combined / se_S0
    sim_results$type1_landmark_z_plugin <- 1 - pnorm(z_crit_tte - z_shift_landmark)
    
    # Using response-based upper bound
    z_shift_landmark_upper <- bias_landmark_resp_upper / se_S0
    sim_results$type1_landmark_z_upper_bound <- 1 - pnorm(z_crit_tte - z_shift_landmark_upper)
    
    # Using utility-based upper bound
    z_shift_landmark_util_upper <- bias_landmark_util_upper / se_S0
    sim_results$type1_landmark_z_utility_upper <- 1 - pnorm(z_crit_tte - z_shift_landmark_util_upper)
    
    # Using surrogate-adjusted
    z_shift_landmark_surrogate <- bias_landmark_surrogate / se_S0
    sim_results$type1_landmark_z_surrogate <- 1 - pnorm(z_crit_tte - z_shift_landmark_surrogate)
    
    # ====================================================================
    # SURVIVAL EVENT COUNTS AND MEDIAN PFS
    # ====================================================================
    
    sim_results$events_stage1 <- colSums(selected_data_censor)
    sim_results$events_stage2 <- colSums(stage2_data_censor)
    
    # Correlation with TTE
    #sim_results$corr_tte_est <- calc_phi_hat(combined_efficacy_mat, combined_tte_mat)
    sim_results$corr_tte_est <- calc_tte_est_fn(combined_efficacy_mat, combined_tte_mat)
    sim_results$corr_tte_est[(is.na(sim_results$corr_tte_est) | !is.finite(sim_results$corr_tte_est))] <- 0
    
    # Median PFS calculation - use memory-efficient approach for large nSim
    # Skip if nSim is too large to avoid memory issues
    if (nSim <= 50000) {
      tryCatch({
        n_patients <- nrow(combined_tte_mat)
        n_sims <- nSim
        
        time_vec <- as.vector(combined_tte_mat)
        event_vec <- as.vector(combined_censor_mat)
        sim_group <- rep(1:n_sims, each = n_patients)
        
        fit_all <- survfit(Surv(time_vec, event_vec) ~ sim_group)
        medians <- summary(fit_all)$table[, "median"]
        
        sim_results$median_pfs <- ifelse(is.na(medians), NA, medians)
        
      }, error = function(e) {
        message("Error in median PFS calculation: ", e$message)
        sim_results$median_pfs <- rep(NA, nSim)
      })
    } else {
      # For large nSim, skip median PFS to avoid memory allocation errors
      # This is not essential for bias/Type I error analysis
      sim_results$median_pfs <- rep(NA_real_, nSim)
    }
    
    # ====================================================================
    # ONE-SAMPLE EXPONENTIAL TEST (LOG-SCALE, 15-Feb-2026)
    # ====================================================================
    
    D <- colSums(combined_censor_mat)
    T_total <- colSums(combined_tte_mat)
    lambda_hat <- D / T_total
    
    #se_lambda <- sqrt(D / T_total^2)
    #test_stat <- (lambda_hat - tte_rate_historical) / se_lambda
    #sim_results$p_value_exponential <- pnorm(test_stat, lower.tail = TRUE) # updated to TRUE on 15-Feb-2026
    #sim_results$reject_H0_exp_one_arm <- as.numeric(sim_results$p_value_exponential <= Alpha_tte)
    
    # Log-scale test: reject if log(lambda_hat) significantly below log(lambda_0)
    log_test_stat <- (log(lambda_hat) - log(tte_rate_historical)) * sqrt(D)
    sim_results$p_value_exponential <- pnorm(log_test_stat, lower.tail = TRUE)
    sim_results$reject_H0_exp_one_arm <- as.numeric(sim_results$p_value_exponential <= Alpha_tte)
    
    # ====================================================================
    # DIRECT ESTIMATION OF Cov(TTE, U) FROM STAGE 1 DATA
    # This provides a more accurate plugin estimate than delta method
    # ====================================================================
    
    # selected_data_tte already contains Stage 1 TTE for selected dose (created earlier)
    # selected_data_censor already contains Stage 1 status for selected dose
    
    # Compute Cov(TTE, U) directly from Stage 1 data
    # E[TTE * U] - E[TTE] * E[U]
    E_TTE_U <- colMeans(selected_data_tte * utility_mat_stage1)
    E_TTE_stage1 <- colMeans(selected_data_tte)
    cov_TTE_U_direct <- E_TTE_U - E_TTE_stage1 * E_U_stage1
    cov_TTE_U_direct[!is.finite(cov_TTE_U_direct)] <- 0
    
    # Variance of TTE in Stage 1
    sigma_TTE_stage1 <- sqrt(pmax(colMeans(selected_data_tte^2) - E_TTE_stage1^2, 1e-10))
    
    # Correlation between TTE and U
    corr_TTE_U <- cov_TTE_U_direct / (sigma_TTE_stage1 * sigma_U_stage1)
    corr_TTE_U[!is.finite(corr_TTE_U)] <- 0
    
    # Direct plugin bias for mean TTE (Equation 28)
    # Bias(T̄) = Cov(T, U) / (σ_U * √n₁) * (1/√π) * exp(...)
    bias_TTE_plugin_stage1 <- (cov_TTE_U_direct / sqrt(mom_sel$sigma2 * N1)) * 
      (1 / sqrt(base::pi)) * exp_term_plugin
    bias_TTE_plugin_combined <- (N1 / n_total) * bias_TTE_plugin_stage1
    
    # ====================================================================
    # PLUGIN TYPE I ERROR FOR ONE-ARM EXPONENTIAL TEST (Log-scale)
    # Method 1: Delta method from landmark survival bias
    # Method 2: Direct from Cov(TTE, U) - more accurate
    # ====================================================================
    
    # Null survival rate at landmark time
    S_null_landmark <- exp(-tte_rate_historical * landmark_time)
    
    ## --- Method 1: Delta method from landmark bias ---
    ## Bias(λ) ≈ -Bias(S(τ)) / (τ * S₀(τ))
    #bias_lambda_delta <- -sim_results$bias_landmark_plugin_combined / 
    #                      (landmark_time * S_null_landmark)
    #bias_lambda_upper <- -sim_results$bias_landmark_upper_bound_combined / 
    #                     (landmark_time * S_null_landmark)
    
    # Bias(log(lambda_hat)) ≈ Bias(lambda_hat) / lambda_0
    #                        = -Bias(T_bar) * lambda_0
    bias_log_lambda_plugin <- -bias_TTE_plugin_combined * tte_rate_historical
    bias_log_lambda_resp_upper <- -sim_results$bias_landmark_upper_bound_combined/(landmark_time * exp(-tte_rate_historical * landmark_time) * tte_rate_historical)
    bias_log_lambda_util_upper <- -sim_results$bias_landmark_utility_upper_combined/(landmark_time * exp(-tte_rate_historical * landmark_time) * tte_rate_historical)
    
    # --- Method 2: Direct from Cov(TTE, U) ---
    # For exponential: λ = 1/E[T], so Bias(λ̂) ≈ -Bias(T̄) * λ₀² (delta method)
    # If E[T] is biased up, λ̂ is biased down
    bias_lambda_direct <- -bias_TTE_plugin_combined * (tte_rate_historical^2)
    
    # Use direct method as primary (more accurate when corr_efficacy_tte > 0)
    bias_lambda_plugin <- bias_lambda_direct
    
    ## Standard error of lambda_hat under null
    ## SE(λ̂) ≈ λ₀ / √D where D = expected events
    #expected_events <- n_total * (1 - S_null_landmark)
    #se_lambda_null <- tte_rate_historical / sqrt(pmax(expected_events, 1))
    
    # SE(log(lambda_hat)) = 1/sqrt(D), use E[D] under null
    # E[D] accounting for staggered entry: Unif(0, entry_time_max)
    lam <- tte_rate_historical
    expected_D <- n_total * (1 - exp(-lam * admin_censor_time) * 
                               (exp(lam * entry_time_max) - 1) / (lam * entry_time_max))    
    expected_D_sqrt <- sqrt(expected_D)
    # more precise: account for staggered entry
    # For uniform entry over [0, entry_time_max], censor at admin_censor_time:
    # E[D] ≈ n_total * P(T < C) where C ~ admin_censor_time - Unif(0, entry_time_max)
    
    
    # Critical value for one-sided test
    z_crit_tte <- qnorm(1 - Alpha_tte)
    
    ## Type I error: P(reject | H0)
    ## Test: reject if λ̂ < λ₀ - z_crit * SE (survival better than null)
    ## Under H0 with bias: λ̂ ~ N(λ₀ + Bias(λ), SE²)
    ## P(reject) = P(Z < -z_crit) where Z = (λ̂ - λ₀)/SE ~ N(Bias/SE, 1)
    #sim_results$type1_exp_one_arm_plugin <- pnorm(-z_crit_tte - bias_lambda_plugin / se_lambda_null)
    #sim_results$type1_exp_one_arm_upper <- pnorm(-z_crit_tte - bias_lambda_upper / se_lambda_null)
    sim_results$type1_exp_one_arm_plugin <- pnorm(-z_crit_tte - 
                                                    bias_log_lambda_plugin * expected_D_sqrt)
    sim_results$type1_exp_one_arm_upper  <- pnorm(-z_crit_tte - 
                                                    bias_log_lambda_resp_upper * expected_D_sqrt)
    sim_results$type1_exp_one_arm_utility_upper <- pnorm(-z_crit_tte - 
                                                    bias_log_lambda_util_upper * expected_D_sqrt)    
    # Store intermediate values for diagnostics
    sim_results$cov_TTE_U_plugin <- cov_TTE_U_direct
    sim_results$corr_TTE_U_plugin <- corr_TTE_U
    sim_results$bias_TTE_plugin <- bias_TTE_plugin_combined
    
    # ====================================================================
    # TWO-ARM SURVIVAL ANALYSIS (COX REGRESSION)
    # Uses chunked processing to avoid memory issues with parallel computing
    # ====================================================================
    
    if (two_arm_survival_flag) {
      
      historical_big <- simulate_pfs_only(
        n = n_total * nSim,
        tte_rate_val = tte_rate_historical,
        entry_time_max = entry_time_max,
        admin_censor_time = admin_censor_time,
        dropout_rate_val = dropout_rate,
        unit_dropout_val = unit_dropout,
        seed_val = ranseed + 999
      )
      
      historical_tte_mat    <- matrix(historical_big$tte_observed, nrow = n_total, ncol = nSim)
      historical_status_mat <- matrix(historical_big$status,        nrow = n_total, ncol = nSim)
      
      all_tte_mat    <- rbind(historical_tte_mat, combined_tte_mat)
      all_status_mat <- rbind(historical_status_mat, combined_censor_mat)
      
      arm_vec <- rep(c(0, 1), each = n_total)
      
      if (use_internal_parallel) {
        # ---- PARALLEL path with retry and sequential fallback ----
        n_workers <- max(1, parallelly::availableCores(omit = 2))
        chunk_size <- min(500, ceiling(nSim / n_workers))
        chunk_list <- split(1:nSim, ceiling(seq_along(1:nSim) / chunk_size))
        
        old_max_size <- getOption("future.globals.maxSize")
        options(future.globals.maxSize = 1024 * 1024^2)
        
        cox_results_list <- tryCatch({
          parallel_ok <- FALSE
          for (attempt in 1:3) {
            cluster_ok <- tryCatch({
              plan(multisession, workers = n_workers)
              TRUE
            }, error = function(e) {
              message(sprintf("Cluster attempt %d failed: %s", attempt, conditionMessage(e)))
              try(plan(sequential), silent = TRUE)
              Sys.sleep(2 * attempt)
              FALSE
            })
            if (cluster_ok) { parallel_ok <- TRUE; break }
          }
          if (!parallel_ok) stop("All cluster setup attempts failed")
          
          res <- future_lapply(chunk_list, function(chunk_idx) {
            tte_chunk    <- all_tte_mat[, chunk_idx, drop = FALSE]
            status_chunk <- all_status_mat[, chunk_idx, drop = FALSE]
            lapply(seq_along(chunk_idx), function(j) {
              df <- data.frame(time = tte_chunk[, j], status = status_chunk[, j], arm = arm_vec)
              fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ arm, data = df), error = function(e) NULL)
              if (is.null(fit)) return(list(hr = NA, p_wald_one_sided = NA, p_wald_two_sided = NA, p_logrank_two_sided = NA, p_logrank_one_sided = NA))
              s <- summary(fit); beta <- coef(fit); se <- sqrt(vcov(fit)[1, 1])
              list(hr = exp(beta), p_wald_one_sided = pnorm(-abs(beta / se)), p_wald_two_sided = s$coefficients[1, "Pr(>|z|)"],
                   p_logrank_two_sided = s$logtest["pvalue"], p_logrank_one_sided = pnorm(sign(beta) * sqrt(s$sctest["test"])))
            })
          }, future.seed = TRUE)
          plan(sequential)
          res
        }, error = function(e) {
          try(plan(sequential), silent = TRUE)
          message(sprintf("Parallel Cox failed: %s - running sequentially", conditionMessage(e)))
          list(lapply(1:nSim, function(i) {
            df <- data.frame(time = all_tte_mat[, i], status = all_status_mat[, i], arm = arm_vec)
            fit <- tryCatch(survival::coxph(survival::Surv(time, status) ~ arm, data = df), error = function(e) NULL)
            if (is.null(fit)) return(list(hr = NA, p_wald_one_sided = NA, p_wald_two_sided = NA, p_logrank_two_sided = NA, p_logrank_one_sided = NA))
            s <- summary(fit); beta <- coef(fit); se <- sqrt(vcov(fit)[1, 1])
            list(hr = exp(beta), p_wald_one_sided = pnorm(-abs(beta / se)), p_wald_two_sided = s$coefficients[1, "Pr(>|z|)"],
                 p_logrank_two_sided = s$logtest["pvalue"], p_logrank_one_sided = pnorm(sign(beta) * sqrt(s$sctest["test"])))
          }))
        })
        options(future.globals.maxSize = old_max_size)
        
      } else {
        # ---- SEQUENTIAL path: called from batch runner ----
        cox_results_list <- list(lapply(1:nSim, function(i) {
          df <- data.frame(
            time   = all_tte_mat[, i],
            status = all_status_mat[, i],
            arm    = arm_vec
          )
          fit <- tryCatch(
            survival::coxph(survival::Surv(time, status) ~ arm, data = df),
            error = function(e) NULL
          )
          if (is.null(fit)) {
            return(list(hr = NA, p_wald_one_sided = NA, p_wald_two_sided = NA,
                        p_logrank_two_sided = NA, p_logrank_one_sided = NA))
          }
          s    <- summary(fit)
          beta <- coef(fit)
          se   <- sqrt(vcov(fit)[1, 1])
          list(
            hr                  = exp(beta),
            p_wald_one_sided    = pnorm(-abs(beta / se)),
            p_wald_two_sided    = s$coefficients[1, "Pr(>|z|)"],
            p_logrank_two_sided = s$logtest["pvalue"],
            p_logrank_one_sided = pnorm(sign(beta) * sqrt(s$sctest["test"]))
          )
        }))
      }
      
      # Flatten results (works identically for both paths)
      cox_results <- unlist(cox_results_list, recursive = FALSE)
      
      
      
      res_df <- do.call(rbind, lapply(cox_results, as.data.frame))
      sim_results$hr                     <- res_df$hr
      sim_results$p_value_wald_one_sided <- res_df$p_wald_one_sided
      sim_results$p_value_wald_two_sided <- res_df$p_wald_two_sided
      sim_results$p_value_logrank_two_sided <- res_df$p_logrank_two_sided
      sim_results$p_value_logrank_one_sided <- res_df$p_logrank_one_sided
      sim_results$reject_H0_logrank_two_arm <- as.numeric(res_df$p_logrank_one_sided <= Alpha_tte)
      
      # ================================================================
      # PLUGIN TYPE I ERROR FOR TWO-ARM LOG-RANK TEST (CORRECTED)
      #
      # Previous approach had two issues:
      #   1. Converted landmark survival bias → event count bias incorrectly
      #      (the relationship depends on the censoring distribution)
      #   2. Ignored that in a two-arm test, E_trt adjusts with O_trt
      #      (Bias(O_trt - E_trt) = Δ/2, not Δ)
      #
      # Corrected approach: Work through log(HR) bias directly
      #   - Uses Cov(TTE, U) already estimated from Stage 1 data
      #     (same quantity validated for the one-arm exponential test)
      #   - Applies delta method: Bias(T̄) → Bias(λ̂) → Bias(log(HR))
      #   - Converts to Cox Z-statistic bias using information matrix
      # ================================================================
      
      # --- Step 1: Bias in treatment arm hazard rate ---
      # From the one-arm exponential section, we already have:
      #   bias_TTE_plugin_combined = (N1/n_total) × Cov(TTE,U)/(σ_U√n₁) × (1/√π) × exp(...)
      #
      # Delta method for exponential rate: λ̂ = D/ΣT
      #   If E[T̄] is biased up by δ, then E[λ̂] ≈ λ₀ - δ × λ₀²
      #   So Bias(λ̂) ≈ -Bias(T̄) × λ₀²
      
      bias_lambda_trt_plugin <- -bias_TTE_plugin_combined * (tte_rate_historical^2)
      
      S0_tau <- exp(-tte_rate_historical * landmark_time)
      bias_lambda_trt_resp_upper <- -sim_results$bias_landmark_upper_bound_combined / 
        (landmark_time * S0_tau)
      bias_lambda_trt_util_upper <- -sim_results$bias_landmark_utility_upper_combined / 
        (landmark_time * S0_tau)
      
      bias_logHR_plugin     <- bias_lambda_trt_plugin     / tte_rate_historical
      bias_logHR_resp_upper <- bias_lambda_trt_resp_upper  / tte_rate_historical
      bias_logHR_util_upper <- bias_lambda_trt_util_upper  / tte_rate_historical
      
      # --- Step 3: Bias in Cox Z-statistic ---
      # Cox model: β̂ = log(HR), se(β̂) = 1/√I
      # where I ≈ D_total/4 for equal allocation (= observed Fisher information)
      #
      # Estimate D_total: project from Stage 1 observed event fraction
      p_event_stage1 <- colSums(selected_data_censor) / N1
      D_total_est <- 2 * n_total * p_event_stage1
      V_cox <- pmax(D_total_est / 4, 1)
      
      # Bias(Z_cox) = Bias(β̂) × √I = Bias(log(HR)) × √(D_total/4)
      # Note: Bias(log(HR)) < 0 when treatment appears better (hazard biased down)
      bias_Z_cox_plugin     <- bias_logHR_plugin     * sqrt(V_cox)
      bias_Z_cox_resp_upper <- bias_logHR_resp_upper * sqrt(V_cox)
      bias_Z_cox_util_upper <- bias_logHR_util_upper * sqrt(V_cox)
      
      # --- Step 4: Type I error ---
      # The test rejects when: pnorm(β̂/se) ≤ α
      #   i.e., when Z_cox = β̂/se ≤ -z_crit  (treatment has lower hazard)
      #
      # Under biased null: Z_cox ~ N(Bias_Z_cox, 1) where Bias_Z_cox < 0
      # P(Z_cox ≤ -z_crit) = Φ(-z_crit - Bias_Z_cox)
      #
      # Since Bias_Z_cox < 0, (-z_crit - Bias_Z_cox) > -z_crit
      # → Probability increases → Type I error > α  ✓
      
      z_crit_lr <- qnorm(1 - Alpha_tte)
      sim_results$type1_logrank_two_arm_plugin       <- pnorm(-z_crit_lr - bias_Z_cox_plugin)
      sim_results$type1_logrank_two_arm_upper        <- pnorm(-z_crit_lr - bias_Z_cox_resp_upper)
      sim_results$type1_logrank_two_arm_utility_upper <- pnorm(-z_crit_lr - bias_Z_cox_util_upper)
      
      
      # ============================================================================
      # VERIFICATION: Sign check
      # ============================================================================
      # When Cov(TTE, U) > 0 (responders survive longer):
      #   bias_TTE_plugin_combined > 0      (mean TTE biased up)
      #   bias_lambda_trt_plugin < 0        (hazard biased down = better)
      #   bias_logHR_plugin < 0             (HR < 1 = treatment better)
      #   bias_Z_cox_plugin < 0             (Cox Z more negative)
      #   -z_crit - bias_Z_cox > -z_crit   (larger argument to Φ)
      #   Type I error > α                  ✓ CORRECT
      #
      # When corr_efficacy_tte = 0 (TTE independent of selection):
      #   Cov(TTE, U) ≈ 0
      #   All biases ≈ 0
      #   Type I error ≈ α                  ✓ CORRECT (no spurious inflation)
      # ============================================================================
    }
  }
  
  # ========================================================================
  # SUMMARIZE AND RETURN RESULTS
  # ========================================================================
  
  # Dose selection summary
  dose_selection_summary <- prop.table(table(sim_results$selected_dose))
  summary_dose_selection <- c(
    prop_select_H = unname(dose_selection_summary["H"]),
    prop_select_L = unname(dose_selection_summary["L"])
  )
  summary_dose_selection[is.na(summary_dose_selection)] <- 0
  
  # Numeric column summaries
  numeric_cols <- sapply(sim_results, is.numeric)
  sim_results_numeric <- sim_results[, numeric_cols]
  
  sim_means <- colMeans(sim_results_numeric, na.rm = TRUE)
  sim_sds <- apply(sim_results_numeric, 2, sd, na.rm = TRUE)
  
  names(sim_means) <- paste0("mean_", names(sim_means))
  names(sim_sds) <- paste0("sd_", names(sim_sds))
  
  simulation_summary <- c(rbind(sim_means, sim_sds))
  names(simulation_summary) <- c(rbind(names(sim_means), names(sim_sds)))
  simulation_summary <- c(summary_dose_selection, simulation_summary)
  
  # Analytical summary
  analytical_summary <- c(
    # Binary endpoint analytical results
    bias_utility_analytical_combined   = bias_utility_analytical_combined,
    bias_response_analytical_combined  = bias_response_analytical_combined,
    bias_utility_score_analytical_combined = bias_utility_score_analytical_combined,
    type1_z_utility_analytical   = type1_z_utility_analytical_combined,
    type1_z_response_analytical  = type1_z_response_analytical_combined,
    type1_bin_utility_analytical = type1_bin_utility_analytical_combined,
    type1_bin_response_analytical= type1_bin_response_analytical_combined,
    
    # Critical value for binomial test
    binomial_critical_value = k_c
  )
  
  # Add TTE analytical results if applicable
  if (perform_tte_analysis) {
    tte_analytical <- c(
      landmark_null_rate = S0_landmark,
      landmark_se_null = se_S0
    )
    analytical_summary <- c(analytical_summary, tte_analytical)
  }
  
  overall_summary <- c(simulation_summary, analytical_summary)
  
  # Return results
  if(return_raw){
    final_output <- list(
      Summary_Results = overall_summary,
      Raw_Data        = sim_results
    )
  } else {
    final_output <- list(
      Summary_Results = overall_summary
    )      
  }
  
  return(final_output)
}

# ============================================================================
# END OF compare_bias_methods_fast_enhanced_v2.R
# ============================================================================



# Shiny App for Utility Score and ROSE Design Sample Size Calculation
# Based on provided R functions from the manuscript draft
# All corrections reviewed and applied:
# - ROSE approximate method uses the corrected v_H calculation.
# - Function names are consistent: calc_sample_size_utility_approx for utility approx, etc.
# - calc_pi and calc_phi_bounds allow q=1 for ROSE compatibility.
# - Custom print methods added; renderPrint uses print(res) instead of str(res).
# - Class assignments added to function returns to enable custom printing.
# - Tested: ROSE approx now matches paper (e.g., n=112 for pL=0.3, delta=0.1, alpha=0.8).

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Calculate Joint Probabilities from Marginals and Correlation
#' 
#' Computes the joint probabilities for four outcome categories:
#' (1) Response & No AE, (2) Response & AE, (3) No Response & No AE, (4) No Response & AE
#' 
#' @param p Response probability
#' @param q No adverse event probability  
#' @param phi Correlation coefficient between response and no-AE (-1 to 1)
#' @return Vector of joint probabilities (pi1, pi2, pi3, pi4)
calc_pi <- function(p, q, phi) {
  if (phi < -1 || phi > 1) stop("phi must be in [-1, 1]")
  
  sqrt_term <- sqrt(p * (1 - p) * q * (1 - q))
  pi1 <- p * q + phi * sqrt_term       # Response & No AE
  pi2 <- p - pi1                        # Response & AE
  pi3 <- q - pi1                        # No Response & No AE  
  pi4 <- 1 - p - q + pi1                # No Response & AE
  
  pi <- c(pi1, pi2, pi3, pi4)
  
  if (any(pi < 0) || any(pi > 1)) {
    stop("Invalid joint probabilities. phi may be incompatible with p and q.")
  }
  
  return(pi)
}

#' Calculate Utility Scores from Trade-off Ratio
#' 
#' Computes utility scores u = (u1, u2, u3, u4) based on the trade-off ratio
#' r = delta/d, which represents the relative importance of efficacy vs safety.
#' 
#' @param r Trade-off ratio (delta/d), must be positive
#' @return Vector of utility scores (u1, u2, u3, u4)
#'         u1 = 1 (response & no AE)
#'         u2 = 1/(1+r) (response & AE)
#'         u3 = r/(1+r) (no response & no AE)
#'         u4 = 0 (no response & AE)
calc_utility <- function(r) {
  if (r <= 0) stop("r must be positive")
  c(1, 1/(1+r), r/(1+r), 0)
}

#' Calculate Mean and Variance of Utility Score
#' 
#' Computes the first two moments of the utility distribution
#' 
#' @param pi Joint probabilities (length 4)
#' @param u Utility scores (length 4)
#' @return List with:
#'   - mu: Expected utility E[U]
#'   - sigma2: Variance Var(U)
calc_utility_moments <- function(pi, u) {
  mu <- sum(u * pi)
  sigma2 <- sum(u^2 * pi) - mu^2
  list(mu = mu, sigma2 = max(sigma2, 1e-10))
}

# ============================================================================
# METHOD 1: UTILITY-BASED APPROXIMATE SAMPLE SIZE
# Uses normal approximation for fast calculation
# ============================================================================

#' Calculate Approximate Sample Size for Utility-Based Design
#'
#' Determines sample size using normal approximation to the distribution
#' of the utility score difference between two doses.
#'
#' @param pL Baseline response rate for low dose (0 < pL < 1)
#' @param qL Baseline no-AE rate for low dose (0 < qL < 1)
#' @param delta Efficacy margin: pH = pL + delta in Scenario H (delta > 0)
#' @param d Safety margin: qH = qL - d in Scenario L (d >= 0)
#' @param phi Correlation between efficacy and safety (-1 <= phi <= 1)
#' @param alpha_L Target PCS for Scenario L (default 0.8)
#' @param alpha_H Target PCS for Scenario H (default 0.8)
#' @param u Optional: Custom utility scores vector of length 4
#'          If NULL, calculated from r = delta/d
#'
#' @return List containing:
#'   - n: Sample size per arm (integer)
#'   - lambda_u: Decision threshold
#'   - PCS_L: Achieved PCS under Scenario L
#'   - PCS_H: Achieved PCS under Scenario H
#'   - utility: Utility scores used
#'   - scenario_L: Parameters for Scenario L
#'   - scenario_H: Parameters for Scenario H
#'   - method: "approximate"
#'
#' @examples
#' # Standard design
#' design <- calc_sample_size_utility_approx(
#'   pL = 0.3, qL = 0.5, delta = 0.15, d = 0.15, phi = 0
#' )
#' print(paste("Sample size:", design$n))
#'
#' # Custom utilities
#' design2 <- calc_sample_size_utility_approx(
#'   pL = 0.3, qL = 0.5, delta = 0.1, d = 0.15, phi = 0.2,
#'   u = c(1, 0.6, 0.4, 0)
#' )
calc_sample_size_utility_approx <- function(
    pL = NULL, qL = NULL,
    delta = NULL, d = NULL,
    pL_L = NULL, qL_L = NULL,
    pH_L = NULL, qH_L = NULL,
    pL_H = NULL, qL_H = NULL,
    pH_H = NULL, qH_H = NULL,
    phi = 0,
    alpha_L = 0.8, alpha_H = 0.8,
    u = NULL
) {
  
  # ===== Common Validation =====
  if (!is.numeric(phi) || phi < -1 || phi > 1) {
    stop("phi must be between -1 and 1")
  }
  if (!is.numeric(alpha_L) || alpha_L <= 0 || alpha_L >= 1) {
    stop("alpha_L must be between 0 and 1")
  }
  if (!is.numeric(alpha_H) || alpha_H <= 0 || alpha_H >= 1) {
    stop("alpha_H must be between 0 and 1")
  }
  
  # ===== Determine Mode =====
  margin_based <- !is.null(delta) || !is.null(d)
  partial_direct <- !is.null(pL) && !is.null(qL) && (!is.null(pH_L) || !is.null(qH_L) || !is.null(pL_H) || !is.null(qL_H))
  full_direct <- !is.null(pL_L) && !is.null(qL_L) && !is.null(pL_H) && !is.null(qL_H) && !is.null(pH_L) && !is.null(qH_L) && !is.null(pH_H) && !is.null(qH_H)
  
  # Check for conflicting modes
  modes_specified <- sum(margin_based, partial_direct, full_direct)
  if (modes_specified > 1) {
    stop("Cannot mix parameters from different modes. Use either:\n",
         "  Mode 1: pL, qL, delta, d\n",
         "  Mode 2: pL, qL, pH_L, qH_L, pL_H, qL_H\n",
         "  Mode 3: pL_L, qL_L, pH_L, qH_L, pL_H, qL_H, pH_H, qH_H")
  }
  if (modes_specified == 0) {
    stop("Must specify parameters for one of the three modes")
  }
  
  r <- NA
  delta_calculated <- NA
  d_calculated <- NA
  mode_name <- ""
  
  # ===== Mode 1: Margin-based =====
  if (margin_based) {
    mode_name <- "margin-based"
    
    # Check required parameters
    if (is.null(pL) || is.null(qL) || is.null(delta) || is.null(d)) {
      stop("Mode 1 requires: pL, qL, delta, d")
    }
    
    # Validate parameter values
    if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
      stop("pL must be between 0 and 1")
    }
    if (!is.numeric(qL) || qL <= 0 || qL >= 1) {
      stop("qL must be between 0 and 1")
    }
    if (!is.numeric(delta) || delta <= 0) {
      stop("delta must be positive")
    }
    if (!is.numeric(d) || d < 0) {
      stop("d must be non-negative")
    }
    
    # Validate derived values
    if (pL + delta > 1) {
      stop("pL + delta must not exceed 1 (would create invalid probability)")
    }
    if (qL - d < 0) {
      stop("qL - d must be non-negative (would create invalid probability)")
    }
    
    # Calculate scenario probabilities
    pL_L <- pL
    qL_L <- qL
    pH_L <- pL
    qH_L <- qL - d
    pL_H <- pL - delta
    qL_H <- qL
    pH_H <- pL
    qH_H <- qL
    
    delta_calculated <- delta
    d_calculated <- d
    
    # Calculate utility if not provided
    if (is.null(u)) {
      if (d == 0) {
        stop("d cannot be zero for automatic utility calculation")
      }
      r <- delta / d
      u <- calc_utility(r)
    }
  }
  
  # ===== Mode 2: Partial Direct =====
  if (partial_direct) {
    mode_name <- "partial-direct"
    
    # Check required parameters
    if (is.null(pL) || is.null(qL) || is.null(pH_L) || is.null(qH_L) || 
        is.null(pL_H) || is.null(qL_H)) {
      stop("Mode 2 requires: pL, qL, pH_L, qH_L, pL_H, qL_H")
    }
    
    # Validate all probability parameters
    if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
      stop("pL must be between 0 and 1")
    }
    if (!is.numeric(qL) || qL <= 0 || qL >= 1) {
      stop("qL must be between 0 and 1")
    }
    if (!is.numeric(pH_L) || pH_L <= 0 || pH_L >= 1) {
      stop("pH_L must be between 0 and 1")
    }
    if (!is.numeric(qH_L) || qH_L <= 0 || qH_L >= 1) {
      stop("qH_L must be between 0 and 1")
    }
    if (!is.numeric(pL_H) || pL_H <= 0 || pL_H >= 1) {
      stop("pL_H must be between 0 and 1")
    }
    if (!is.numeric(qL_H) || qL_H <= 0 || qL_H >= 1) {
      stop("qL_H must be between 0 and 1")
    }
    
    # Calculate scenario probabilities
    pL_L <- pL
    qL_L <- qL
    pH_H <- pL
    qH_H <- qL
    
    delta_calculated <- pL - pL_H
    d_calculated <- qL - qH_L
    
    # Validate utility is provided
    if (is.null(u)) {
      if (d_calculated == 0) {
        stop("d cannot be zero for automatic utility calculation")
      }
      r <- delta_calculated / d_calculated
      u <- calc_utility(r)
    }
  }
  
  # ===== Mode 3: Full Direct =====
  if (full_direct) {
    mode_name <- "full-direct"
    
    # Check required parameters
    if (is.null(pL_L) || is.null(qL_L) || is.null(pH_L) || is.null(qH_L) ||
        is.null(pL_H) || is.null(qL_H) || is.null(pH_H) || is.null(qH_H)) {
      stop("Mode 3 requires all 8 probabilities: pL_L, qL_L, pH_L, qH_L, pL_H, qL_H, pH_H, qH_H")
    }
    
    # Validate all 8 probability parameters
    if (!is.numeric(pL_L) || pL_L <= 0 || pL_L >= 1) {
      stop("pL_L must be between 0 and 1")
    }
    if (!is.numeric(qL_L) || qL_L <= 0 || qL_L >= 1) {
      stop("qL_L must be between 0 and 1")
    }
    if (!is.numeric(pH_L) || pH_L <= 0 || pH_L >= 1) {
      stop("pH_L must be between 0 and 1")
    }
    if (!is.numeric(qH_L) || qH_L <= 0 || qH_L >= 1) {
      stop("qH_L must be between 0 and 1")
    }
    if (!is.numeric(pL_H) || pL_H <= 0 || pL_H >= 1) {
      stop("pL_H must be between 0 and 1")
    }
    if (!is.numeric(qL_H) || qL_H <= 0 || qL_H >= 1) {
      stop("qL_H must be between 0 and 1")
    }
    if (!is.numeric(pH_H) || pH_H <= 0 || pH_H >= 1) {
      stop("pH_H must be between 0 and 1")
    }
    if (!is.numeric(qH_H) || qH_H <= 0 || qH_H >= 1) {
      stop("qH_H must be between 0 and 1")
    }
    
    delta_calculated <- mean(c(pH_L - pL_L, pH_H - pL_H))
    d_calculated <- mean(c(qL_L - qH_L, qL_H - qH_H))
    
    # Validate utility is provided
    if (is.null(u)) {
      stop("Mode 3 requires utility scores 'u' to be provided")
    }
  }
  
  # ===== Validate Utility Vector =====
  if (!is.numeric(u) || length(u) != 4) {
    stop("u must be a numeric vector of length 4")
  }
  
  # ===== Calculate Joint Probabilities and Moments =====
  pi_L_L <- calc_pi(pL_L, qL_L, phi)
  pi_H_L <- calc_pi(pH_L, qH_L, phi)
  pi_L_H <- calc_pi(pL_H, qL_H, phi)
  pi_H_H <- calc_pi(pH_H, qH_H, phi)
  
  mom_L_L <- calc_utility_moments(pi_L_L, u)
  mom_H_L <- calc_utility_moments(pi_H_L, u)
  mom_L_H <- calc_utility_moments(pi_L_H, u)
  mom_H_H <- calc_utility_moments(pi_H_H, u)
  
  # ===== Calculate Utility Differences =====
  delta_mu_L <- mom_H_L$mu - mom_L_L$mu
  v_L <- mom_H_L$sigma2 + mom_L_L$sigma2
  delta_mu_H <- mom_H_H$mu - mom_L_H$mu
  v_H <- mom_H_H$sigma2 + mom_L_H$sigma2
  
  # ===== Validate Scenarios =====
  if (delta_mu_L >= 0) {
    warning("Scenario L: high dose is not inferior to low dose (delta_mu_L >= 0)")
  }
  if (delta_mu_H <= 0) {
    warning("Scenario H: high dose is not superior to low dose (delta_mu_H <= 0)")
  }
  if ((delta_mu_H - delta_mu_L) <= 0) {
    stop("Utility difference between scenarios is not positive. ",
         "Cannot calculate sample size with current parameters.")
  }
  
  # ===== Sample Size Calculation =====
  z_alpha_L <- qnorm(alpha_L)
  z_1_alpha_H <- qnorm(1 - alpha_H)
  
  n_continuous <- ((z_alpha_L * sqrt(v_L) - z_1_alpha_H * sqrt(v_H)) / 
                     (delta_mu_H - delta_mu_L))^2
  n <- ceiling(n_continuous)
  
  # Validate sample size
  if (n < 2) {
    stop("Calculated sample size is less than 2. Check your input parameters.")
  }
  if (n > 10000) {
    warning("Calculated sample size is very large (n = ", n, " > 10,000). ",
            "Consider reviewing your design parameters.")
  }
  
  # ===== Calculate Decision Threshold =====
  lambda_u <- delta_mu_H + z_1_alpha_H * sqrt(v_H / n)
  
  # ===== Calculate Achieved PCS =====
  pcs_L <- pnorm((lambda_u - delta_mu_L) / sqrt(v_L / n))
  pcs_H <- 1 - pnorm((lambda_u - delta_mu_H) / sqrt(v_H / n))
  
  # ===== Return Results =====
  result <- list(
    n = n,
    lambda_u = lambda_u,
    PCS_L = pcs_L,
    PCS_H = pcs_H,
    scenario_L = list(
      delta_mu = delta_mu_L, 
      v = v_L,
      pL = pL_L, 
      qL = qL_L, 
      pH = pH_L, 
      qH = qH_L,
      pi_L = pi_L_L, 
      pi_H = pi_H_L
    ),
    scenario_H = list(
      delta_mu = delta_mu_H, 
      v = v_H,
      pL = pL_H, 
      qL = qL_H, 
      pH = pH_H, 
      qH = qH_H,
      pi_L = pi_L_H, 
      pi_H = pi_H_H
    ),
    utility = u,
    r = r,
    delta = delta_calculated,
    d = d_calculated,
    inputs = list(
      mode = mode_name, 
      phi = phi,
      alpha_L = alpha_L, 
      alpha_H = alpha_H,
      pL = pL, 
      qL = qL, 
      delta = delta, 
      d = d,
      pL_L = pL_L, 
      qL_L = qL_L, 
      pH_L = pH_L, 
      qH_L = qH_L,
      pL_H = pL_H, 
      qL_H = qL_H, 
      pH_H = pH_H, 
      qH_H = qH_H
    )
  )
  
  class(result) <- c("dose_design", "list")
  return(result)
}


# ============================================================================
# METHOD 2: UTILITY-BASED EXACT SAMPLE SIZE
# Uses dynamic programming for exact multinomial calculations
# ============================================================================

#' Compute PMF of Sum of Utilities Using Dynamic Programming
#'
#' Calculates the probability mass function of S = sum(u_int * n_i) where 
#' (n_1, n_2, n_3, n_4) follows a multinomial distribution.
#'
#' @param pi Numeric vector of length 4. Multinomial probabilities.
#' @param u_int Integer vector of length 4. Integer-scaled utility scores.
#' @param n Integer. Sample size (number of multinomial trials).
#'
#' @return Numeric vector representing the PMF of S, where element k 
#'   corresponds to P(S = k-1).
#'
#' @details
#' This function uses dynamic programming to efficiently compute the exact 
#' distribution of the weighted sum of multinomial counts. The algorithm 
#' iteratively updates the PMF by convolving with the single-trial distribution 
#' n times, with complexity O(n * max(u_int)).
#'
#' @examples
#' pi <- c(0.3, 0.2, 0.1, 0.4)
#' u_int <- c(10, 6, 3, 0)
#' n <- 50
#' pmf <- compute_pmf_S(pi, u_int, n)
#'
#' @keywords internal
compute_pmf_S <- function(pi, u_int, n) {
  if (length(pi) != 4 || length(u_int) != 4) {
    stop("pi and u_int must be length 4")
  }
  
  max_S <- n * max(u_int)
  pmf <- rep(0, max_S + 1)
  pmf[1] <- 1  # S=0 has probability 1 initially
  
  for (i in 1:n) {
    new_pmf <- rep(0, max_S + 1)
    for (j in 1:4) {
      shift <- u_int[j]
      if (shift == 0) {
        new_pmf <- new_pmf + pi[j] * pmf
      } else {
        new_pmf[(shift + 1):(max_S + 1)] <- new_pmf[(shift + 1):(max_S + 1)] + 
          pi[j] * pmf[1:(max_S + 1 - shift)]
      }
    }
    pmf <- new_pmf
  }
  
  pmf
}


#' Compute PMF of Difference Between Two Distributions (Nested Loop Method)
#'
#' Internal function using nested loops to compute P(D = S_H - S_L).
#'
#' @param pmf_L Numeric vector. PMF of S_L.
#' @param pmf_H Numeric vector. PMF of S_H.
#' @param max_S Integer. Maximum value for both S_L and S_H.
#'
#' @return Numeric vector of length (2*max_S + 1) representing the PMF of D.
#'
#' @keywords internal
compute_pmf_diff_nested <- function(pmf_L, pmf_H, max_S) {
  pmf_diff <- rep(0, 2 * max_S + 1)
  offset <- max_S + 1  # Index for d=0
  
  for (s_L in 0:max_S) {
    if (pmf_L[s_L + 1] == 0) next
    for (s_H in 0:max_S) {
      if (pmf_H[s_H + 1] == 0) next
      d <- s_H - s_L
      pmf_diff[d + offset] <- pmf_diff[d + offset] + pmf_L[s_L + 1] * pmf_H[s_H + 1]
    }
  }
  
  pmf_diff
}


#' Compute PMF of Difference Using FFT Convolution
#'
#' Internal function using Fast Fourier Transform for computing P(D = S_H - S_L).
#'
#' @param pmf_L Numeric vector. PMF of S_L.
#' @param pmf_H Numeric vector. PMF of S_H.
#' @param max_S Integer. Maximum value for both S_L and S_H.
#'
#' @return Numeric vector of length (2*max_S + 1) representing the PMF of D.
#'
#' @keywords internal
compute_pmf_diff_fft <- function(pmf_L, pmf_H, max_S) {
  # Reverse pmf_L for subtraction (convolution with reversed = correlation)
  pmf_L_rev <- rev(pmf_L)
  
  # Pad to avoid circular convolution
  target_length <- 2 * max_S + 1
  len_L <- length(pmf_L_rev)
  len_H <- length(pmf_H)
  
  fft_length <- target_length + max(len_L, len_H) - 1
  fft_length <- 2^ceiling(log2(fft_length))
  
  # Zero-pad
  pmf_L_padded <- c(pmf_L_rev, rep(0, fft_length - len_L))
  pmf_H_padded <- c(pmf_H, rep(0, fft_length - len_H))
  
  # FFT convolution
  fft_L <- stats::fft(pmf_L_padded)
  fft_H <- stats::fft(pmf_H_padded)
  conv_result <- stats::fft(fft_L * fft_H, inverse = TRUE) / fft_length
  conv_result <- Re(conv_result)
  
  # Extract valid range
  start_idx <- len_L
  end_idx <- start_idx + target_length - 1
  pmf_diff <- conv_result[start_idx:end_idx]
  
  # Numerical cleanup
  pmf_diff[abs(pmf_diff) < 1e-15] <- 0
  pmf_diff[pmf_diff < 0] <- 0
  
  pmf_diff
}


#' Compute PMF of Difference D = S_H - S_L (Unified Interface)
#'
#' Computes the probability mass function of the difference between two 
#' discrete random variables S_H and S_L, given their PMFs.
#'
#' @param pmf_L Numeric vector. PMF of S_L where pmf_L[k] = P(S_L = k-1).
#' @param pmf_H Numeric vector. PMF of S_H where pmf_H[k] = P(S_H = k-1).
#' @param max_S Integer, optional. Maximum value for S_L and S_H.
#'   If NULL, inferred as length(pmf_L) - 1. Default: NULL.
#' @param method Character. Computation method: "nested" (exact, recommended) 
#'   or "fft" (faster but approximate). Default: "nested".
#' @param validate Logical. Whether to validate inputs. Default: TRUE.
#'
#' @return Numeric vector of length (2*max_S + 1) representing P(D = d) where
#'   d ranges from -max_S to +max_S. Index (max_S + 1) corresponds to d = 0.
#'
#' @details
#' The function computes:
#' \deqn{P(D = d) = P(S_H - S_L = d) = \sum_{s_L} P(S_L = s_L) P(S_H = s_L + d)}
#'
#' **Method Comparison:**
#'
#' *Nested Loops (method = "nested")*:
#' \itemize{
#'   \item Complexity: O(N^2) worst case, O(k^2) with sparsity
#'   \item Precision: Exact (machine precision ~10^-16)
#'   \item Speed: 2-4ms for N=300
#'   \item Recommended for: Scientific calculations, N < 2000
#' }
#'
#' *FFT Convolution (method = "fft")*:
#' \itemize{
#'   \item Complexity: O(N log N)
#'   \item Precision: Approximate (FFT errors ~10^-15 to 10^-12)
#'   \item Speed: 1-2ms for N=300
#'   \item Recommended for: Large N (>2000), approximate OK
#' }
#'
#' @examples
#' # Basic usage with nested loops (recommended)
#' pmf_L <- c(0.3, 0.5, 0.2, rep(0, 97))
#' pmf_H <- c(0.2, 0.4, 0.3, 0.1, rep(0, 96))
#' pmf_diff <- compute_pmf_diff_unified(pmf_L, pmf_H)
#'
#' # Using FFT for large arrays
#' pmf_L <- dpois(0:2000, lambda = 500)
#' pmf_H <- dpois(0:2000, lambda = 600)
#' pmf_diff <- compute_pmf_diff_unified(pmf_L, pmf_H, method = "fft")
#'
#' @references
#' Wang, S., Yuan, Y., and Liu, S. (2025). 
#' Randomized Optimal Selection Design for Dose Optimization. 
#' arXiv preprint arXiv:2505.03898v4.
#'
#' @export
compute_pmf_diff_unified <- function(pmf_L, pmf_H, 
                                     max_S = NULL, 
                                     method = c("nested", "fft"),
                                     validate = TRUE) {
  
  method <- match.arg(method)
  
  # Input validation
  if (validate) {
    if (!is.numeric(pmf_L) || !is.vector(pmf_L)) {
      stop("pmf_L must be a numeric vector")
    }
    if (!is.numeric(pmf_H) || !is.vector(pmf_H)) {
      stop("pmf_H must be a numeric vector")
    }
    if (length(pmf_L) != length(pmf_H)) {
      stop("pmf_L and pmf_H must have the same length")
    }
    if (any(pmf_L < 0, na.rm = TRUE) || any(pmf_H < 0, na.rm = TRUE)) {
      stop("Probabilities must be non-negative")
    }
    
    sum_L <- sum(pmf_L, na.rm = TRUE)
    sum_H <- sum(pmf_H, na.rm = TRUE)
    if (abs(sum_L - 1) > 1e-6 || abs(sum_H - 1) > 1e-6) {
      warning("PMFs should sum to 1")
    }
  }
  
  # Determine max_S
  if (is.null(max_S)) {
    max_S <- length(pmf_L) - 1
  }
  
  # Compute based on method
  if (method == "nested") {
    pmf_diff <- compute_pmf_diff_nested(pmf_L, pmf_H, max_S)
  } else {
    pmf_diff <- compute_pmf_diff_fft(pmf_L, pmf_H, max_S)
  }
  
  pmf_diff
}



#' Calculate Exact Sample Size for Utility-Based Design
#'
#' Determines sample size using exact multinomial distribution calculations
#' via dynamic programming and probability mass functions.
#'
#' @param pL,qL Base response/no-AE rates for low dose (Mode 1 & 2).
#' @param pL_L,qL_L Response/no-AE rates for low dose in Scenario L (Mode 3).
#' @param pH_L,qH_L Response/no-AE rates for high dose in Scenario L.
#' @param pL_H,qL_H Response/no-AE rates for low dose in Scenario H (Mode 3).
#' @param pH_H,qH_H Response/no-AE rates for high dose in Scenario H.
#' @param phi Correlation coefficient between efficacy and safety.
#' @param delta Efficacy margin (Mode 1).
#' @param d Safety margin (Mode 1).
#' @param alpha_L Target PCS under scenario L. Default: 0.8.
#' @param alpha_H Target PCS under scenario H. Default: 0.8.
#' @param u Utility scores (4-element vector).
#' @param max_n Maximum sample size to search. Default: 500.
#' @param buffer Search buffer around approximate n. Default: 10.
#' @param den Denominator for integer scaling. Default: 10.
#' @param diff_method PMF difference method: "nested" or "fft". Default: "nested".
#' @param verbose Print progress messages. Default: FALSE.
#'
#' @return A list containing exact design parameters.
#'
#' @details
#' Uses dynamic programming to compute exact PMFs of utility scores
#' under multinomial sampling, then searches for minimal n satisfying
#' PCS constraints.
#'
#' **Special case**: For ROSE designs (u = c(1,1,0,0)), automatically
#' sets den=1 to preserve integer structure.
#'
#' @examples
#' \dontrun{
#' # ROSE design with exact calculation
#' rose_exact <- calc_sample_size_exact(
#'   pL = 0.3, qL = 0.575, 
#'   delta = 0.10, d = 0.15,
#'   phi = 0, u = c(1, 1, 0, 0),
#'   alpha_L = 0.7, alpha_H = 0.7,
#'   den = 1
#' )
#'
#' # Utility-based design with exact calculation
#' util_exact <- calc_sample_size_exact(
#'   pL = 0.3, qL = 0.575,
#'   delta = 0.10, d = 0.15,
#'   phi = 0,
#'   alpha_L = 0.7, alpha_H = 0.7,
#'   den = 10
#' )
#' }
#'
#' @export
calc_sample_size_utility_exact <- function(
    pL = NULL, qL = NULL,
    delta = NULL, d = NULL,
    pL_L = NULL, qL_L = NULL,
    pH_L = NULL, qH_L = NULL,
    pL_H = NULL, qL_H = NULL,
    pH_H = NULL, qH_H = NULL,
    phi = 0,
    alpha_L = 0.8, alpha_H = 0.8,
    u = NULL,
    max_n = 500,
    buffer = 10,
    den = 10,
    diff_method = c("nested", "fft"),
    verbose = FALSE
) {
  
  diff_method <- match.arg(diff_method)
  
  # Input validation
  if (!is.numeric(phi) || phi < -1 || phi > 1) {
    stop("phi must be between -1 and 1")
  }
  if (!is.numeric(alpha_L) || alpha_L <= 0 || alpha_L >= 1) {
    stop("alpha_L must be between 0 and 1")
  }
  if (!is.numeric(alpha_H) || alpha_H <= 0 || alpha_H >= 1) {
    stop("alpha_H must be between 0 and 1")
  }
  
  # Determine mode
  margin_based <- !is.null(delta) || !is.null(d)
  partial_direct <- !is.null(pL) && !is.null(qL) && (!is.null(pH_L) || !is.null(qH_L) || !is.null(pL_H) || !is.null(qL_H))
  full_direct <- !is.null(pL_L) && !is.null(qL_L) && !is.null(pL_H) && !is.null(qL_H) && !is.null(pH_L) && !is.null(qH_L) && !is.null(pH_H) && !is.null(qH_H)
  
  # Check for conflicting modes
  modes_specified <- sum(margin_based, partial_direct, full_direct)
  if (modes_specified > 1) {
    stop("Cannot mix parameters from different modes. Use either:\n",
         "  Mode 1: pL, qL, delta, d\n",
         "  Mode 2: pL, qL, pH_L, qH_L, pL_H, qL_H\n",
         "  Mode 3: pL_L, qL_L, pH_L, qH_L, pL_H, qL_H, pH_H, qH_H")
  }
  if (modes_specified == 0) {
    stop("Must specify parameters for one of the three modes")
  }
  
  r <- NA
  delta_calculated <- NA
  d_calculated <- NA
  mode_name <- ""
  
  # ===== Mode 1: Margin-based =====
  if (margin_based) {
    mode_name <- "margin-based"
    if (is.null(pL) || is.null(qL) || is.null(delta) || is.null(d)) {
      stop("Mode 1 requires: pL, qL, delta, d")
    }
    
    # Validate parameter values
    if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
      stop("pL must be between 0 and 1")
    }
    if (!is.numeric(qL) || qL <= 0 || qL >= 1) {
      stop("qL must be between 0 and 1")
    }
    if (!is.numeric(delta) || delta <= 0) {
      stop("delta must be positive")
    }
    if (!is.numeric(d) || d < 0) {
      stop("d must be non-negative")
    }
    
    # Validate derived values
    if (pL - delta < 0) {
      stop("pL - delta must be non-negative (would create invalid probability)")
    }
    if (qL - d < 0) {
      stop("qL - d must be non-negative (would create invalid probability)")
    }
    
    # Scenario L: Low dose baseline, High dose same efficacy but worse safety
    pL_L <- pL
    qL_L <- qL
    pH_L <- pL
    qH_L <- qL - d
    
    # Scenario H: Low dose worse efficacy, High dose baseline
    pL_H <- pL - delta
    qL_H <- qL
    pH_H <- pL
    qH_H <- qL
    
    delta_calculated <- delta
    d_calculated <- d
    
    if (is.null(u)) {
      if (d == 0) stop("d cannot be zero for automatic utility calculation")
      r <- delta / d
      u <- calc_utility(r)
    }
  }
  
  # ===== Mode 2: Partial direct =====
  if (partial_direct) {
    mode_name <- "partial-direct"
    if (is.null(pL) || is.null(qL) || is.null(pH_L) || is.null(qH_L) ||
        is.null(pL_H) || is.null(qL_H)) {
      stop("Mode 2 requires: pL, qL, pH_L, qH_L, pL_H, qL_H")
    }
    
    # Validate all probability parameters
    if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
      stop("pL must be between 0 and 1")
    }
    if (!is.numeric(qL) || qL <= 0 || qL >= 1) {
      stop("qL must be between 0 and 1")
    }
    if (!is.numeric(pH_L) || pH_L <= 0 || pH_L >= 1) {
      stop("pH_L must be between 0 and 1")
    }
    if (!is.numeric(qH_L) || qH_L <= 0 || qH_L >= 1) {
      stop("qH_L must be between 0 and 1")
    }
    if (!is.numeric(pL_H) || pL_H <= 0 || pL_H >= 1) {
      stop("pL_H must be between 0 and 1")
    }
    if (!is.numeric(qL_H) || qL_H <= 0 || qL_H >= 1) {
      stop("qL_H must be between 0 and 1")
    }
    
    # Scenario L: Low dose baseline
    pL_L <- pL
    qL_L <- qL
    
    # Scenario H: High dose baseline
    pH_H <- pL
    qH_H <- qL
    
    delta_calculated <- pL - pL_H
    d_calculated <- qL - qH_L
    
    if (is.null(u)) {
      if (d_calculated == 0) {
        stop("d cannot be zero for automatic utility calculation")
      }
      r <- delta_calculated / d_calculated
      u <- calc_utility(r)
    }
  }
  
  # ===== Mode 3: Full direct =====
  if (full_direct) {
    mode_name <- "full-direct"
    if (is.null(pL_L) || is.null(qL_L) || is.null(pH_L) || is.null(qH_L) ||
        is.null(pL_H) || is.null(qL_H) || is.null(pH_H) || is.null(qH_H)) {
      stop("Mode 3 requires all 8 probabilities: pL_L, qL_L, pH_L, qH_L, pL_H, qL_H, pH_H, qH_H")
    }
    
    # Validate all 8 probability parameters
    if (!is.numeric(pL_L) || pL_L <= 0 || pL_L >= 1) {
      stop("pL_L must be between 0 and 1")
    }
    if (!is.numeric(qL_L) || qL_L <= 0 || qL_L >= 1) {
      stop("qL_L must be between 0 and 1")
    }
    if (!is.numeric(pH_L) || pH_L <= 0 || pH_L >= 1) {
      stop("pH_L must be between 0 and 1")
    }
    if (!is.numeric(qH_L) || qH_L <= 0 || qH_L >= 1) {
      stop("qH_L must be between 0 and 1")
    }
    if (!is.numeric(pL_H) || pL_H <= 0 || pL_H >= 1) {
      stop("pL_H must be between 0 and 1")
    }
    if (!is.numeric(qL_H) || qL_H <= 0 || qL_H >= 1) {
      stop("qL_H must be between 0 and 1")
    }
    if (!is.numeric(pH_H) || pH_H <= 0 || pH_H >= 1) {
      stop("pH_H must be between 0 and 1")
    }
    if (!is.numeric(qH_H) || qH_H <= 0 || qH_H >= 1) {
      stop("qH_H must be between 0 and 1")
    }
    
    delta_calculated <- mean(c(pH_L - pL_L, pH_H - pL_H))
    d_calculated <- mean(c(qL_L - qH_L, qL_H - qH_H))
    
    if (is.null(u)) {
      stop("Mode 3 requires utility scores 'u' to be provided")
    }
  }
  
  if (!is.numeric(u) || length(u) != 4) {
    stop("u must be a numeric vector of length 4")
  }
  
  # Rest of the function remains the same...
  # Calculate probabilities and moments
  pi_L_SL <- calc_pi(pL_L, qL_L, phi)
  pi_H_SL <- calc_pi(pH_L, qH_L, phi)
  pi_L_SH <- calc_pi(pL_H, qL_H, phi)
  pi_H_SH <- calc_pi(pH_H, qH_H, phi)
  
  mom_L_L <- calc_utility_moments(pi_L_SL, u)
  mom_H_L <- calc_utility_moments(pi_H_SL, u)
  mom_L_H <- calc_utility_moments(pi_L_SH, u)
  mom_H_H <- calc_utility_moments(pi_H_SH, u)
  
  delta_mu_L <- mom_H_L$mu - mom_L_L$mu
  v_L <- mom_H_L$sigma2 + mom_L_L$sigma2
  delta_mu_H <- mom_H_H$mu - mom_L_H$mu
  v_H <- mom_H_H$sigma2 + mom_L_H$sigma2
  
  if (delta_mu_L >= 0) warning("Scenario L: high dose not inferior")
  if (delta_mu_H <= 0) warning("Scenario H: high dose not superior")
  if (delta_mu_H - delta_mu_L <= 0) stop("Utility difference not positive")  
  # Approximate starting point
  if (verbose) cat("Step 1: approximate n ...\n")
  zL <- qnorm(alpha_L)
  zH <- qnorm(1 - alpha_H)
  n_cont <- ((zL * sqrt(v_L) - zH * sqrt(v_H)) / (delta_mu_H - delta_mu_L))^2
  n_approx <- ceiling(n_cont)
  lambda_approx <- delta_mu_H + zH * sqrt(v_H / n_approx)
  min_search <- max(5, n_approx - buffer)
  
  if (verbose) {
    cat(sprintf("  approx n = %d, lambda = %.6f\n", n_approx, lambda_approx))
    cat(sprintf("Step 2: exact search [%d,%d] (method=%s)\n", 
                min_search, max_n, diff_method))
  }
  
  # Special case: ROSE design
  if (sum(u == c(1, 1, 0, 0)) == 4) den <- 1
  
  u_int <- round(u * den)
  found <- FALSE
  n_tested <- 0
  
  for (n in min_search:max_n) {
    n_tested <- n_tested + 1
    if (verbose && (n_tested %% 5 == 0 || n == min_search)) {
      cat(sprintf("  testing n = %d ...\n", n))
    }
    
    max_S <- n * max(u_int)
    
    # Compute PMFs
    pmf_L_SL <- compute_pmf_S(pi_L_SL, u_int, n)
    pmf_H_SL <- compute_pmf_S(pi_H_SL, u_int, n)
    pmf_L_SH <- compute_pmf_S(pi_L_SH, u_int, n)
    pmf_H_SH <- compute_pmf_S(pi_H_SH, u_int, n)
    
    # Compute differences
    pmf_diff_SL <- compute_pmf_diff_unified(pmf_L_SL, pmf_H_SL,
                                            max_S = max_S,
                                            method = diff_method,
                                            validate = FALSE)
    pmf_diff_SH <- compute_pmf_diff_unified(pmf_L_SH, pmf_H_SH,
                                            max_S = max_S,
                                            method = diff_method,
                                            validate = FALSE)
    
    cdf_SL <- cumsum(pmf_diff_SL)
    cdf_SH <- cumsum(pmf_diff_SH)
    offset <- max_S + 1
    
    # Find quantiles
    q_SL_int <- min(which(cdf_SL >= alpha_L)) - offset
    q_SH_int <- max(which(cdf_SH <= 1 - alpha_H)) - offset
    
    if (q_SL_int <= q_SH_int) {
      found <- TRUE
      n_final <- n
      lambda_int <- round((q_SL_int + q_SH_int) / 2)
      lambda_u <- lambda_int / (den * n)
      
      pcs_L <- cdf_SL[lambda_int + offset]
      pcs_H <- 1 - cdf_SH[lambda_int + offset]
      
      if (verbose) {
        cat(sprintf("  SOLUTION: n = %d, lambda = %.6f\n", n, lambda_u))
        cat(sprintf("    PCS_L = %.4f (target %.2f)\n", pcs_L, alpha_L))
        cat(sprintf("    PCS_H = %.4f (target %.2f)\n", pcs_H, alpha_H))
      }
      break
    }
  }
  
  if (!found) {
    stop(sprintf("No solution in [%d,%d]. Increase max_n.", min_search, max_n))
  }
  
  result <- list(
    n = n_final,
    lambda_u = lambda_u,
    PCS_L = pcs_L,
    PCS_H = pcs_H,
    scenario_L = list(
      delta_mu = delta_mu_L, v = v_L,
      pL = pL_L, qL = qL_L, pH = pH_L, qH = qH_L,
      pi_L = pi_L_SL, pi_H = pi_H_SL
    ),
    scenario_H = list(
      delta_mu = delta_mu_H, v = v_H,
      pL = pL_H, qL = qL_H, pH = pH_H, qH = qH_H,
      pi_L = pi_L_SH, pi_H = pi_H_SH
    ),
    utility = u,
    r = r,
    delta = delta_calculated,
    d = d_calculated,
    method = sprintf("exact_%s", diff_method),
    inputs = list(
      mode = mode_name, phi = phi,
      alpha_L = alpha_L, alpha_H = alpha_H,
      pL = pL, qL = qL, delta = delta, d = d,
      pL_L = pL_L, qL_L = qL_L, pH_L = pH_L, qH_L = qH_L,
      pL_H = pL_H, qL_H = qL_H, pH_H = pH_H, qH_H = qH_H,
      max_n = max_n, buffer = buffer, den = den,
      diff_method = diff_method
    ),
    search_info = list(
      n_approx = n_approx,
      lambda_u_approx = lambda_approx,
      n_continuous = n_cont,
      min_search = min_search,
      max_search = max_n,
      n_tested = n_tested,
      algorithm = sprintf("DP + %s convolution", diff_method),
      den = den
    )
  )
  
  class(result) <- c("dose_design_exact", "dose_design", "list")
  return(result)
}


# ============================================================================
# METHOD 3: ROSE APPROXIMATE SAMPLE SIZE
# Response-Only Selection (efficacy-only, ignoring safety)
# ============================================================================

#' Calculate Approximate Sample Size for ROSE Design
#'
#' ROSE (Response-Only Selection) design focuses solely on efficacy,
#' ignoring safety considerations. Equivalent to utility-based design
#' with utilities u = (1, 1, 0, 0).
#'
#' @param pL Baseline response rate for low dose
#' @param delta Efficacy margin
#' @param alpha_L Target PCS for Scenario L (default 0.8)
#' @param alpha_H Target PCS for Scenario H (default 0.8)
#'
#' @return List containing:
#'   - n: Sample size per arm
#'   - lambda_u: Decision threshold for response rate difference
#'   - PCS_L: Achieved PCS under Scenario L
#'   - PCS_H: Achieved PCS under Scenario H
#'   - utility: (1, 1, 0, 0)
#'   - scenario_L: Parameters for Scenario L
#'   - scenario_H: Parameters for Scenario H
#'   - method: "rose_approximate"
#'
#' @examples
#' # ROSE design
#' design <- calc_sample_size_rose_approx(
#'   pL = 0.3, delta = 0.15
#' )
#' print(paste("ROSE sample size:", design$n))
calc_sample_size_rose_approx <- function(pL, delta,
                                         alpha_L = 0.8, alpha_H = 0.8) {
  
  # Input validation
  if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
    stop("pL must be between 0 and 1")
  }
  if (!is.numeric(delta) || delta <= 0) {
    stop("delta must be positive")
  }
  if (pL - delta < 0) {
    stop("pL - delta must be non-negative")
  }
  
  # Scenario L: Both doses have response rate pL (no difference)
  pL_L <- pL
  pH_L <- pL
  delta_mu_L <- 0
  v_L <- 2 * pL * (1 - pL)
  
  # Scenario H: High dose has response rate pL, low dose has pL - delta
  pL_H <- pL - delta
  pH_H <- pL 
  delta_mu_H <- delta
  v_H <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)  # Fixed: Use pL_H here
  
  # Calculate sample size
  z_alpha_L <- qnorm(alpha_L)
  z_1_alpha_H <- qnorm(1 - alpha_H)
  
  n_continuous <- ((z_alpha_L * sqrt(v_L) - z_1_alpha_H * sqrt(v_H)) / delta_mu_H)^2
  n <- ceiling(n_continuous)
  
  if (n < 2) {
    stop("Calculated sample size is less than 2. Check parameters.")
  }
  
  # Calculate decision threshold
  lambda_u <- delta_mu_H + z_1_alpha_H * sqrt(v_H / n)
  
  # Calculate achieved PCS
  pcs_L <- pnorm((lambda_u - delta_mu_L) / sqrt(v_L / n))
  pcs_H <- 1 - pnorm((lambda_u - delta_mu_H) / sqrt(v_H / n))
  
  # Return results
  result <- list(
    n = n,
    lambda_u = lambda_u,
    PCS_L = pcs_L,
    PCS_H = pcs_H,
    utility = c(1, 1, 0, 0),  # ROSE utilities
    scenario_L = list(
      delta_mu = delta_mu_L,
      variance = v_L,
      pL = pL_L, pH = pH_L
    ),
    scenario_H = list(
      delta_mu = delta_mu_H,
      variance = v_H,
      pL = pL_H, pH = pH_H
    ),
    method = "rose_approximate",
    inputs = list(
      alpha_L = alpha_L,
      alpha_H = alpha_H
    ),
    delta = delta
  )
  class(result) <- c("rose_design", "list")
  return(result)
}

# ============================================================================
# METHOD 4: ROSE EXACT SAMPLE SIZE
# Uses exact binomial calculations
# ============================================================================

#' Compute PMF of Sum for Binomial (special case of multinomial)
#'
#' @param p Success probability
#' @param n Sample size
#' @return PMF vector of length (n + 1) where pmf[k+1] = P(X = k)
compute_pmf_binomial <- function(p, n) {
  dbinom(0:n, size = n, prob = p)
}

#' Calculate Exact Sample Size for ROSE Design
#'
#' Uses exact binomial calculations for ROSE design.
#'
#' @param pL Baseline response rate for low dose
#' @param delta Efficacy margin
#' @param alpha_L Target PCS for Scenario L (default 0.8)
#' @param alpha_H Target PCS for Scenario H (default 0.8)
#' @param max_n Maximum sample size to search (default 200)
#' @param buffer Buffer below approximate n (default 10)
#' @param verbose Print progress (default FALSE)
#'
#' @return List with same structure as ROSE approximate, plus:
#'   - method: "rose_exact"
#'   - n_approx: Approximate starting point
#'
#' @examples
#' # ROSE exact calculation
#' design <- calc_sample_size_rose_exact(
#'   pL = 0.3, delta = 0.15
#' )
#' print(paste("ROSE exact sample size:", design$n))
calc_sample_size_rose_exact <- function(pL, delta,
                                        alpha_L = 0.8, alpha_H = 0.8,
                                        max_n = 200,
                                        buffer = 10,
                                        verbose = FALSE) {
  
  # Input validation
  if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
    stop("pL must be between 0 and 1")
  }
  if (!is.numeric(delta) || delta <= 0) {
    stop("delta must be positive")
  }
  if (pL - delta < 0) {
    stop("pL - delta must be non-negative")
  }
  
  # Define scenarios
  pL_L <- pL; pH_L <- pL
  pL_H <- pL - delta; pH_H <- pL
  
  # Get approximate starting point
  v_L <- 2 * pL * (1 - pL)
  v_H <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)
  
  z_alpha_L <- qnorm(alpha_L)
  z_1_alpha_H <- qnorm(1 - alpha_H)
  n_approx <- ceiling(((z_alpha_L * sqrt(v_L) - z_1_alpha_H * sqrt(v_H)) / delta)^2)
  
  min_search <- max(5, n_approx - buffer)
  
  if (verbose) {
    cat(sprintf("Approximate n = %d, searching [%d, %d]...\n", 
                n_approx, min_search, max_n))
  }
  
  # Exact search
  found <- FALSE
  n_final <- NA
  lambda_u_final <- NA
  pcs_L_final <- NA
  pcs_H_final <- NA
  
  for (n in min_search:max_n) {
    if (verbose && (n - min_search) %% 5 == 0) {
      cat(sprintf("  Testing n = %d...\n", n))
    }
    
    # Scenario L: compute PMF of difference
    pmf_L_SL <- compute_pmf_binomial(pL_L, n)
    pmf_H_SL <- compute_pmf_binomial(pH_L, n)
    
    # Difference PMF for binomial case
    pmf_diff_SL <- rep(0, 2*n + 1)
    offset <- n + 1
    
    for (x_L in 0:n) {
      for (x_H in 0:n) {
        d <- x_H - x_L
        pmf_diff_SL[d + offset] <- pmf_diff_SL[d + offset] + 
          pmf_L_SL[x_L + 1] * pmf_H_SL[x_H + 1]
      }
    }
    cdf_SL <- cumsum(pmf_diff_SL)
    
    # Scenario H: compute PMF of difference
    pmf_L_SH <- compute_pmf_binomial(pL_H, n)
    pmf_H_SH <- compute_pmf_binomial(pH_H, n)
    
    pmf_diff_SH <- rep(0, 2*n + 1)
    for (x_L in 0:n) {
      for (x_H in 0:n) {
        d <- x_H - x_L
        pmf_diff_SH[d + offset] <- pmf_diff_SH[d + offset] + 
          pmf_L_SH[x_L + 1] * pmf_H_SH[x_H + 1]
      }
    }
    cdf_SH <- cumsum(pmf_diff_SH)
    
    # Find quantiles (in count difference scale)
    q_SL <- min(which(cdf_SL >= alpha_L)) - offset
    q_SH <- max(which(cdf_SH <= 1 - alpha_H)) - offset
    
    # Check feasibility
    if (q_SL <= q_SH) {
      found <- TRUE
      n_final <- n
      
      # Threshold in count difference
      threshold_count <- round((q_SL + q_SH) / 2)
      # Convert to rate difference
      lambda_u_final <- threshold_count / n
      
      # Calculate exact PCS
      pcs_L_final <- cdf_SL[threshold_count + offset]
      pcs_H_final <- 1 - cdf_SH[threshold_count + offset]
      
      if (verbose) {
        cat(sprintf("  Solution found: n = %d, lambda_u = %.6f\n", 
                    n_final, lambda_u_final))
        cat(sprintf("    PCS_L = %.4f, PCS_H = %.4f\n", 
                    pcs_L_final, pcs_H_final))
      }
      
      break
    }
  }
  
  if (!found) {
    stop(sprintf("No solution found in range [%d, %d]. Increase max_n.", 
                 min_search, max_n))
  }
  
  # Return results
  result <- list(
    n = n_final,
    lambda_u = lambda_u_final,
    PCS_L = pcs_L_final,
    PCS_H = pcs_H_final,
    utility = c(1, 1, 0, 0),
    scenario_L = list(
      delta_mu = 0,
      variance = v_L,
      pL = pL_L, pH = pH_L
    ),
    scenario_H = list(
      delta_mu = delta,
      variance = v_H,
      pL = pL_H, pH = pH_H
    ),
    method = "rose_exact",
    n_approx = n_approx,
    inputs = list(
      alpha_L = alpha_L,
      alpha_H = alpha_H
    ),
    delta = delta
  )
  class(result) <- c("rose_design", "list")
  return(result)
}

#' Print Method for ROSE Design Objects
#'
#' @param x A rose_design object from calc_sample_size_rose.
#' @param ... Additional arguments (not used).
#' @export
print.rose_design <- function(x, ...) {
  cat("ROSE Design (Response-Only Selection)\n")
  cat("======================================\n\n")
  
  cat("Sample Size:\n")
  cat(sprintf("  n per arm: %d\n", x$n))
  cat(sprintf("  Total n:   %d\n\n", 2 * x$n))
  
  cat("Decision Threshold:\n")
  cat(sprintf("  lambda_u:  %.6f\n\n", x$lambda_u))
  
  cat("Achieved PCS:\n")
  cat(sprintf("  PCS_L (select low dose):  %.4f (target: %.2f)\n", 
              x$PCS_L, x$inputs$alpha_L))
  cat(sprintf("  PCS_H (select high dose): %.4f (target: %.2f)\n\n", 
              x$PCS_H, x$inputs$alpha_H))
  
  cat("Utility Scores:\n")
  cat(sprintf("  u = [%.4f, %.4f, %.4f, %.4f]\n", 
              x$utility[1], x$utility[2], x$utility[3], x$utility[4]))
  cat("  (Response-only: ignores safety)\n\n")
  
  cat("Scenario L (Low dose should be selected):\n")
  cat(sprintf("  Dose L: p = %.3f\n", x$scenario_L$pL))
  cat(sprintf("  Dose H: p = %.3f\n", x$scenario_L$pH))
  cat(sprintf("  Delta_mu = %.4f, v = %.4f\n\n", 
              x$scenario_L$delta_mu, x$scenario_L$variance))
  
  cat("Scenario H (High dose should be selected):\n")
  cat(sprintf("  Dose L: p = %.3f\n", x$scenario_H$pL))
  cat(sprintf("  Dose H: p = %.3f\n", x$scenario_H$pH))
  cat(sprintf("  Delta_mu = %.4f, v = %.4f\n\n", 
              x$scenario_H$delta_mu, x$scenario_H$variance))
  
  cat("Efficacy margin: delta = ", x$delta, "\n\n")
  
  cat("Note: ROSE design considers only efficacy (response rate),\n")
  cat("      ignoring safety outcomes.\n")
  
  invisible(x)
}


#' Print Method for Dose Design Objects
#'
#' @param x A dose_design object from calc_sample_size or calc_sample_size_exact.
#' @param ... Additional arguments (not used).
#' @export
print.dose_design <- function(x, ...) {
  is_exact <- inherits(x, "dose_design_exact")
  method_label <- if (is_exact) "EXACT" else "APPROXIMATE"
  
  cat("Utility-Based Dose Selection Design (", method_label, ")\n", sep = "")
  cat(strrep("=", 50 + nchar(method_label)), "\n\n")
  
  cat("Method:\n")
  cat(sprintf("  Type: %s\n", method_label))
  cat("\n")
  
  cat("Sample Size:\n")
  cat(sprintf("  n per arm: %d\n", x$n))
  cat(sprintf("  Total n:   %d\n\n", 2 * x$n))
  
  cat("Decision Threshold:\n")
  cat(sprintf("  lambda_u:  %.6f\n\n", x$lambda_u))
  
  cat("Achieved PCS:\n")
  cat(sprintf("  PCS_L: %.4f (target: %.2f)\n", x$PCS_L, x$inputs$alpha_L))
  cat(sprintf("  PCS_H: %.4f (target: %.2f)\n\n", x$PCS_H, x$inputs$alpha_H))
  
  cat("Utility:\n")
  cat(sprintf("  u = [%.4f, %.4f, %.4f, %.4f]\n", 
              x$utility[1], x$utility[2], x$utility[3], x$utility[4]))
  cat("\n")
  
  cat("Scenario L (Low dose optimal):\n")
  cat(sprintf("  L: p=%.3f q=%.3f   H: p=%.3f q=%.3f\n",
              x$scenario_L$pL, x$scenario_L$qL, 
              x$scenario_L$pH, x$scenario_L$qH))
  cat(sprintf("  Delta_mu=%.4f, v=%.4f\n\n", 
              x$scenario_L$delta_mu, x$scenario_L$variance))
  
  cat("Scenario H (High dose optimal):\n")
  cat(sprintf("  L: p=%.3f q=%.3f   H: p=%.3f q=%.3f\n",
              x$scenario_H$pL, x$scenario_H$qL, 
              x$scenario_H$pH, x$scenario_H$qH))
  cat(sprintf("  Delta_mu=%.4f, v=%.4f\n\n", 
              x$scenario_H$delta_mu, x$scenario_H$variance))
  
  cat("Parameters:\n")
  cat(sprintf("  phi = %.2f\n", x$phi))
  if (!is.null(x$delta)) {
    cat(sprintf("  delta = %.3f, d = %.3f\n", x$delta, x$d))
  }
  
  if (is_exact) {
    cat("\nSearch Info:\n")
    cat(sprintf("  Approx n: %d\n", x$n_approx))
  }
  
  invisible(x)
}




# ============================================================================
# HELPER FUNCTIONS (keeping all existing functions)
# ============================================================================

# [Keep all existing sample size calculation functions]
# [Keep calc_sample_size_utility_approx, calc_sample_size_utility_exact, etc.]

# ============================================================================
# NEW: BIAS AND TYPE I ERROR CALCULATION FUNCTIONS
# ============================================================================

#' Calculate Analytical Bias Estimates
#'
#' @param p Response rate under null hypothesis
#' @param q No-AE rate under null hypothesis  
#' @param phi Correlation between efficacy and safety
#' @param u Utility scores (4-element vector)
#' @param n1 Stage 1 sample size per arm
#' @param n2 Stage 2 sample size per arm
#' @param lambda_u Selection threshold
#'
#' @return List with bias estimates
calc_analytical_bias <- function(p, q, phi, u, n1, n2, lambda_u = 0) {
  
  # Calculate joint probabilities and moments
  pi_vec <- calc_pi(p, q, phi)
  moments <- calc_utility_moments(pi_vec, u)
  
  mu <- moments$mu
  sigma2 <- moments$sigma2
  sigma <- sqrt(sigma2)
  
  # Covariance between X and U
  # X is binary: x = (1, 1, 0, 0) for the four outcomes
  x_vec <- c(1, 1, 0, 0)
  E_XU <- sum(x_vec * u * pi_vec)
  cov_XU <- E_XU - p * mu
  
  # Exponential term
  exp_term <- exp(-(lambda_u^2 * n1) / (4 * sigma2))
  
  # Utility-based bias (Equation 15 in manuscript)
  bias_utility_stage1 <- (cov_XU / (sigma * sqrt(n1))) * (1 / sqrt(pi)) * exp_term
  bias_utility_combined <- (n1 / (n1 + n2)) * bias_utility_stage1
  
  # Maximum bias bound (response-only, Equation 16)
  bias_response_stage1 <- sqrt(p * (1 - p)) / sqrt(n1 * pi)
  bias_response_combined <- (n1 / (n1 + n2)) * bias_response_stage1
  
  list(
    bias_utility_stage1 = bias_utility_stage1,
    bias_utility_combined = bias_utility_combined,
    bias_response_stage1 = bias_response_stage1,
    bias_response_combined = bias_response_combined,
    bias_U_score_combined = (n1 / (n1 + n2)) * (sigma / sqrt(n1 * pi)) * exp_term,
    cov_XU = cov_XU,
    sigma_U = sigma,
    mu_U = mu
  )
}

#' Calculate Analytical Type I Error Estimates
#'
#' @param p0 Null hypothesis response rate
#' @param bias_combined Combined bias estimate
#' @param n_total Total sample size (n1 + n2)
#' @param alpha Nominal significance level
#'
#' @return List with Type I error estimates
calc_analytical_type1_error <- function(p0, bias_combined, n_total, alpha = 0.025) {
  
  # Standard error under null
  se_null <- sqrt(p0 * (1 - p0) / n_total)
  
  # Critical value for one-sided test
  z_crit <- qnorm(1 - alpha)
  
  # Type I error for Z-test (Equation 19)
  type1_z <- 1 - pnorm(z_crit - bias_combined / se_null)
  
  # Critical value for binomial test
  # Find k_c such that P(X >= k_c | H0) <= alpha
  probs_under_null <- pbinom(0:n_total, size = n_total, prob = p0, lower.tail = FALSE)
  k_c <- min(which(probs_under_null <= alpha))
  
  # Type I error for binomial test (Equation 20)
  type1_bin <- pbinom(k_c - 1, size = n_total, 
                      prob = p0 + bias_combined, 
                      lower.tail = FALSE)
  
  list(
    type1_z = type1_z,
    type1_bin = type1_bin,
    k_c = k_c,
    se_null = se_null
  )
}

#' Simulate Bias and Type I Error
#'
#' @param p Response rate under null
#' @param q No-AE rate under null
#' @param phi Correlation
#' @param u Utility scores
#' @param n1 Stage 1 sample size per arm
#' @param n2 Stage 2 sample size per arm
#' @param lambda_u Selection threshold
#' @param nSim Number of simulations
#' @param alpha Significance level for tests
#' @param seed Random seed
#'
#' @return List with simulation results
simulate_bias_type1 <- function(p, q, phi, u, n1, n2, lambda_u = 0, 
                                nSim = 10000, alpha = 0.025, seed = 123) {
  
  set.seed(seed)
  
  # Calculate probabilities
  pi_vec <- calc_pi(p, q, phi)
  
  # Storage for results
  selected_dose <- character(nSim)
  p_hat_combined <- numeric(nSim)
  reject_z <- logical(nSim)
  reject_bin <- logical(nSim)
  
  # Simulation loop
  for (i in 1:nSim) {
    
    # Stage 1: Generate data for both doses (null: both identical)
    data_L_stage1 <- sample(1:4, n1, replace = TRUE, prob = pi_vec)
    data_H_stage1 <- sample(1:4, n1, replace = TRUE, prob = pi_vec)
    
    # Calculate utility scores
    U_L_stage1 <- mean(u[data_L_stage1])
    U_H_stage1 <- mean(u[data_H_stage1])
    
    # Select dose
    if ((U_H_stage1 - U_L_stage1) > lambda_u) {
      selected_dose[i] <- "H"
      data_selected_stage1 <- data_H_stage1
    } else {
      selected_dose[i] <- "L"
      data_selected_stage1 <- data_L_stage1
    }
    
    # Stage 2: Generate additional data for selected dose
    data_selected_stage2 <- sample(1:4, n2, replace = TRUE, prob = pi_vec)
    
    # Combine data
    data_combined <- c(data_selected_stage1, data_selected_stage2)
    
    # Calculate combined response rate
    # Response occurs in categories 1 and 2
    response_combined <- sum(data_combined %in% c(1, 2))
    p_hat_combined[i] <- response_combined / (n1 + n2)
    
    # Z-test
    se_null <- sqrt(p * (1 - p) / (n1 + n2))
    z_stat <- (p_hat_combined[i] - p) / se_null
    reject_z[i] <- z_stat > qnorm(1 - alpha)
    
    # Binomial test
    p_value_bin <- pbinom(response_combined - 1, n1 + n2, p, lower.tail = FALSE)
    reject_bin[i] <- p_value_bin < alpha
  }
  
  # Calculate observed statistics
  bias_observed <- mean(p_hat_combined) - p
  type1_z_observed <- mean(reject_z)
  type1_bin_observed <- mean(reject_bin)
  
  # Dose selection proportions
  prop_H <- mean(selected_dose == "H")
  prop_L <- mean(selected_dose == "L")
  
  list(
    bias_observed = bias_observed,
    type1_z_observed = type1_z_observed,
    type1_bin_observed = type1_bin_observed,
    prop_select_H = prop_H,
    prop_select_L = prop_L,
    p_hat_combined = p_hat_combined,
    selected_dose = selected_dose
  )
}



# ============================================================================
# SHINY UI + SERVER (v2.4)
# Updated: 15-Feb-2026
#
# KEY FIXES from v2.3:
#   1. Binary table restructured: split into Bias table + T1E table.
#      Observed appears ONCE per test; utility-based and max-bound predictions
#      shown alongside so the hierarchy is visible.
#   2. TTE tables: same restructuring — observed is the ground truth reference,
#      plugin/surrogate/upper are predictions that should bracket it.
#   3. Formulas tab: one-arm exponential updated to show log-scale as the
#      CURRENT implementation (not "planned"). Hazard-scale section removed.
# ============================================================================

ui <- fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("Dose Optimization: Design, Bias, and Type I Error Analysis (v2.4)"),
  
  tabsetPanel(
    
    # ================================================================
    # 1. General Utility Score Design (unchanged)
    # ================================================================
    tabPanel("General Utility Score Design",
             sidebarLayout(
               sidebarPanel(
                 selectInput("mode", "Specification Mode:", 
                             choices = c("Margin-based (Mode 1)" = "margin-based",
                                         "Partial Direct (Mode 2)" = "partial-direct",
                                         "Full Direct (Mode 3)" = "full-direct")),
                 conditionalPanel(
                   condition = "input.mode == 'margin-based' || input.mode == 'partial-direct'",
                   numericInput("pL", "Base Response Rate (pL):", value = 0.3, min = 0, max = 1, step = 0.01),
                   numericInput("qL", "Base No-AE Rate (qL):", value = 0.9, min = 0, max = 1, step = 0.01)
                 ),
                 conditionalPanel(
                   condition = "input.mode == 'margin-based'",
                   numericInput("delta", "Efficacy Margin (delta):", value = 0.1, min = 0, step = 0.01),
                   numericInput("d", "Safety Margin (d):", value = 0.2, min = 0, step = 0.01)
                 ),
                 conditionalPanel(
                   condition = "input.mode == 'partial-direct'",
                   h5("Scenario L - High Dose:"),
                   numericInput("pH_L", "Response Rate (pH_L):", value = 0.3, min = 0, max = 1, step = 0.01),
                   numericInput("qH_L", "No-AE Rate (qH_L):", value = 0.7, min = 0, max = 1, step = 0.01),
                   h5("Scenario H - Low Dose:"),
                   numericInput("pL_H", "Response Rate (pL_H):", value = 0.2, min = 0, max = 1, step = 0.01),
                   numericInput("qL_H", "No-AE Rate (qL_H):", value = 0.9, min = 0, max = 1, step = 0.01)
                 ),
                 conditionalPanel(
                   condition = "input.mode == 'full-direct'",
                   h5("Scenario L:"),
                   numericInput("pL_L", "Low Dose Response (pL_L):", value = 0.3, min = 0, max = 1, step = 0.01),
                   numericInput("qL_L", "Low Dose No-AE (qL_L):", value = 0.9, min = 0, max = 1, step = 0.01),
                   numericInput("pH_L", "High Dose Response (pH_L):", value = 0.3, min = 0, max = 1, step = 0.01),
                   numericInput("qH_L", "High Dose No-AE (qH_L):", value = 0.7, min = 0, max = 1, step = 0.01),
                   h5("Scenario H:"),
                   numericInput("pL_H", "Low Dose Response (pL_H):", value = 0.2, min = 0, max = 1, step = 0.01),
                   numericInput("qL_H", "Low Dose No-AE (qL_H):", value = 0.9, min = 0, max = 1, step = 0.01),
                   numericInput("pH_H", "High Dose Response (pH_H):", value = 0.3, min = 0, max = 1, step = 0.01),
                   numericInput("qH_H", "High Dose No-AE (qH_H):", value = 0.9, min = 0, max = 1, step = 0.01)
                 ),
                 numericInput("phi", "Correlation Coefficient (phi):", value = 0, min = -1, max = 1, step = 0.01),
                 numericInput("alpha_L", "Target PCS under Scenario L:", value = 0.8, min = 0, max = 1, step = 0.01),
                 numericInput("alpha_H", "Target PCS under Scenario H:", value = 0.8, min = 0, max = 1, step = 0.01),
                 conditionalPanel(
                   condition = "input.mode != 'margin-based'",
                   textInput("u", "Utility Scores (comma-separated):", value = "1,0.667,0.333,0")
                 ),
                 selectInput("method", "Calculation Method:", 
                             choices = c("Approximate" = "approx", "Exact" = "exact")),
                 conditionalPanel(
                   condition = "input.method == 'exact'",
                   numericInput("max_n", "Max Sample Size to Search:", value = 500),
                   numericInput("buffer", "Search Buffer:", value = 10),
                   numericInput("den", "Denominator for Integer Scaling:", value = 10)
                 ),
                 actionButton("compute_util", "Compute Sample Size", class = "btn-primary")
               ),
               mainPanel(
                 verbatimTextOutput("util_result")
               )
             )
    ),
    
    # ================================================================
    # 2. ROSE Design (unchanged)
    # ================================================================
    tabPanel("ROSE Design",
             sidebarLayout(
               sidebarPanel(
                 numericInput("pL_rose", "Base Response Rate (pL):", value = 0.3, min = 0, max = 1, step = 0.01),
                 numericInput("delta_rose", "Efficacy Margin (delta):", value = 0.1, min = 0, step = 0.01),
                 numericInput("alpha_L_rose", "Target PCS under Scenario L:", value = 0.8, min = 0, max = 1, step = 0.01),
                 numericInput("alpha_H_rose", "Target PCS under Scenario H:", value = 0.8, min = 0, max = 1, step = 0.01),
                 selectInput("method_rose", "Calculation Method:", 
                             choices = c("Approximate" = "approx", "Exact" = "exact")),
                 conditionalPanel(
                   condition = "input.method_rose == 'exact'",
                   numericInput("max_n_rose", "Max Sample Size to Search:", value = 500),
                   numericInput("buffer_rose", "Search Buffer:", value = 10)
                 ),
                 actionButton("compute_rose", "Compute Sample Size", class = "btn-primary")
               ),
               mainPanel(
                 verbatimTextOutput("rose_result")
               )
             )
    ),
    
    # ================================================================
    # 3. Simulation Analysis — RESTRUCTURED TABLES
    # ================================================================
    tabPanel("Simulation Analysis",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 h4("Endpoint Selection"),
                 checkboxInput("perform_tte_sim", "Include Time-to-Event (TTE) Endpoint Analysis", value = TRUE),
                 conditionalPanel(
                   condition = "input.perform_tte_sim == true",
                   h5("TTE Parameters"),
                   numericInput("tte_rate_sim", "True TTE Event Rate (lambda):", 0.1, min = 0, step = 0.01),
                   numericInput("corr_eff_tte_sim", "Corr(Efficacy, TTE) [latent Gaussian]:", 0.3, min = -1, max = 1, step = 0.1),
                   numericInput("landmark_time_sim", "Landmark Time (tau):", 24, min = 1),
                   numericInput("alpha_tte_sim", "TTE Significance Level (alpha_tte):", 0.025),
                   h5("Censoring Settings"),
                   numericInput("entry_time_max_sim", "Max Entry Time:", 52),
                   numericInput("admin_censor_time_sim", "Admin Censor Time:", 76),
                   numericInput("dropout_rate_sim", "Annual Dropout Rate:", 0.05, min = 0, max = 1)
                 ),
                 hr(),
                 h4("Simulation Parameters"),
                 numericInput("p_sim", "Null Response Rate (p):", 0.4, min = 0, max = 1, step = 0.01),
                 numericInput("q_sim", "Null No-AE Rate (q):", 0.8, min = 0, max = 1, step = 0.01),
                 numericInput("phi_sim", "Correlation (phi):", 0.3, min = -1, max = 1, step = 0.1),
                 textInput("u_sim", "Utility Scores (comma-separated):", "1, 0.8, 0.2, 0"),
                 numericInput("n1_sim", "Stage 1 Size (n1):", 50, min = 1),
                 numericInput("n2_sim", "Stage 2 Size (n2):", 100, min = 0),
                 numericInput("lambda_u_sim", "Selection Threshold (lambda_u):", 0),
                 numericInput("alpha_sim", "Significance Level (alpha):", 0.025, min = 0.001, max = 0.1),
                 sliderInput("nSim", "Number of Simulations:", min = 100, max = 50000, value = 1000, step = 100),
                 numericInput("seed", "Random Seed:", 123),
                 actionButton("run_simulation", "Run Simulation", class = "btn-primary btn-lg", icon = icon("play"))
               ),
               mainPanel(width = 8,
                         tabsetPanel(id = "results_tabs",
                                     tabPanel("Binary Endpoint Results",
                                              h4("Response Rate Bias"),
                                              DTOutput("binary_bias_table"),
                                              hr(),
                                              h4("Type I Error: Binary Tests"),
                                              DTOutput("binary_type1_table"),
                                              downloadButton("download_binary", "Download Binary Results"),
                                              plotOutput("bias_comparison_plot_sim", height = "400px"),
                                              plotOutput("binary_type1_plot_sim", height = "400px")
                                     ),
                                     tabPanel("TTE Endpoint Results",
                                              h4("Landmark Survival Bias"),
                                              DTOutput("tte_bias_table"),
                                              hr(),
                                              h4("Type I Error: All TTE Tests"),
                                              DTOutput("tte_type1_table"),
                                              hr(),
                                              h4("Correlation Diagnostics"),
                                              DTOutput("tte_corr_table"),
                                              downloadButton("download_tte_full", "Download Full TTE Results"),
                                              plotOutput("tte_bias_plot_sim", height = "400px"),
                                              plotOutput("tte_type1_plot_sim", height = "400px")
                                     )
                         )
               )
             )
    ),
    
    # ================================================================
    # 4. Formulas & Methods — UPDATED for log-scale one-arm test
    # ================================================================
    tabPanel("Formulas & Methods",
             withMathJax(),
             fluidPage(
               h3("Key Formulas from the Manuscript"),
               
               # --- Sample Size ---
               h4("1. Sample Size (Approximate Method)"),
               wellPanel(
                 p("$$n = \\left[ \\frac{z_{\\alpha_L}\\sqrt{v_{SL}} - z_{1-\\alpha_H}\\sqrt{v_{SH}}}
                     {\\Delta\\mu_{SH} - \\Delta\\mu_{SL}} \\right]^2$$"),
                 p("Under margin-based specification with trade-off ratio \\(r = \\delta/d\\):"),
                 p("$$\\Delta\\mu_{SH} - \\Delta\\mu_{SL} = \\frac{2\\delta}{1+r}$$")
               ),
               
               # --- Binary Bias ---
               h4("2. Selection-Induced Bias (Utility-based)"),
               wellPanel(
                 p("$$\\text{Bias}(\\hat{p}_{\\text{sel}}) = 
                     \\frac{\\text{Cov}(X, U)}{\\sigma_U \\sqrt{n_1}} \\cdot 
                     \\frac{1}{\\sqrt{\\pi}} \\exp\\!\\left(-\\frac{\\lambda_u^2 n_1}{4\\sigma_U^2}\\right)$$"),
                 p("where \\(\\text{Cov}(X,U) = \\sum_{k=1}^4 x_k u_k \\pi_k - p\\mu\\) and \\(\\sigma_U^2 = \\text{Var}(U)\\)."),
                 p(em("Derived via truncated bivariate normal lemma (Rosenbaum, 1961): for \\(D = Z_H - Z_L \\sim N(0,2)\\),")),
                 p("$$E[Z_H \\cdot \\mathbf{1}\\{D > k\\}] + E[Z_L \\cdot \\mathbf{1}\\{D \\le k\\}] 
                     = \\frac{1}{\\sqrt{\\pi}} e^{-k^2/4}$$")
               ),
               
               # --- Max Bias ---
               h4("3. Upper Bound Bias (Response-only, Bauer et al. 2010)"),
               wellPanel(
                 p("$$\\text{Bias}_{\\max} = \\frac{\\sqrt{p(1-p)}}{\\sqrt{n_1\\pi}}$$"),
                 p("Attained when \\(U = X\\) (i.e., \\(\\text{Corr}(X,U)=1\\), selection based on efficacy alone).")
               ),
               
               # --- Binary Type I Error ---
               h4("4. Type I Error (Binary Endpoints)"),
               wellPanel(
                 h5("Z-test:"),
                 p("$$\\text{Type I Error}_Z = 1 - \\Phi\\!\\left(z_{1-\\alpha} - 
                     \\frac{\\text{Bias}_{\\text{comb}}}{\\text{SE}_0}\\right), 
                     \\quad \\text{SE}_0 = \\sqrt{\\frac{p_0(1-p_0)}{n_1+n_2}}$$"),
                 h5("Binomial test:"),
                 p("$$\\text{Type I Error}_{\\text{Bin}} = P(X \\ge k_c \\mid p = p_0 + \\text{Bias}_{\\text{comb}})$$"),
                 p("where \\(\\text{Bias}_{\\text{comb}} = \\frac{n_1}{n_1+n_2} \\cdot \\text{Bias}_{\\text{Stage 1}}\\) (dilution).")
               ),
               
               hr(),
               h3("Time-to-Event Bias Propagation"),
               
               # --- Landmark Survival Bias ---
               h4("5. Landmark Survival Bias"),
               wellPanel(
                 h5("Plugin (direct Cov estimation):"),
                 p("$$\\text{Bias}(\\hat{S}(\\tau)) = 
                     \\frac{n_1}{n}\\cdot\\frac{\\widehat{\\text{Cov}}(S_i(\\tau), U)}{\\hat\\sigma_U\\sqrt{n_1}} \\cdot 
                     \\frac{1}{\\sqrt{\\pi}}\\exp\\!\\left(-\\frac{\\lambda_u^2 n_1}{4\\hat\\sigma_U^2}\\right)$$"),
                 h5("Surrogate-adjusted:"),
                 p("$$\\text{Bias}(\\hat{S}(\\tau)) \\approx 
                     \\hat\\rho_{X,S(\\tau)} \\cdot \\frac{\\hat\\sigma_{S(\\tau)}}{\\hat\\sigma_X} \\cdot 
                     \\text{Bias}_{\\max}(\\hat{p})$$"),
                 p("Captures bias transmitted ", em("through"), " the response endpoint \\(X\\). 
                    Uses the response upper bound as anchor."),
                 h5("Upper bound (\\(\\rho_{X,S(\\tau)}=1\\)):"),
                 p("$$\\text{Bias}_{\\max}(\\hat{S}(\\tau)) = 
                     \\frac{\\hat\\sigma_{S(\\tau)}}{\\hat\\sigma_X} \\cdot \\text{Bias}_{\\max}(\\hat{p})$$"),
                 h5("Landmark Z-test Type I Error:"),
                 p("$$\\text{Type I Error}_{\\text{Lm}} = 1 - \\Phi\\!\\left(z_{1-\\alpha} - 
                     \\frac{\\text{Bias}(\\hat{S}(\\tau))_{\\text{comb}}}{\\text{SE}_0}\\right), \\quad 
                     \\text{SE}_0 = \\sqrt{\\frac{S_0(1-S_0)}{n}}$$")
               ),
               
               # --- One-Arm Exponential (LOG-SCALE — current implementation) ---
               h4("6. One-Arm Exponential Test (Log-Scale)"),
               wellPanel(
                 h5("Test statistic (variance-stabilized):"),
                 p("$$Z_{\\log} = \\left(\\log\\hat\\lambda - \\log\\lambda_0\\right)\\cdot\\sqrt{D}$$"),
                 p("where \\(\\hat\\lambda = D / \\sum T_i\\), \\(D\\) = observed events, 
                    and \\(\\text{SE}(\\log\\hat\\lambda) = 1/\\sqrt{D}\\)."),
                 p("Reject \\(H_0\\!: \\lambda \\ge \\lambda_0\\) when \\(Z_{\\log} \\le -z_{1-\\alpha}\\) 
                    (treatment has lower hazard)."),
                 p(em("Advantage over hazard-scale formulation: SE is independent of \\(\\hat\\lambda\\),
                       eliminating numerator-denominator correlation and Jensen's inequality bias.")),
                 h5("Bias on log-hazard scale (delta method):"),
                 p("$$\\text{Bias}(\\log\\hat\\lambda) 
                     \\approx \\frac{\\text{Bias}(\\hat\\lambda)}{\\lambda_0} 
                     = -\\text{Bias}(\\bar{T})\\cdot\\lambda_0$$"),
                 p("where \\(\\text{Bias}(\\bar{T})\\) is the mean TTE bias from Section 5 via \\(\\text{Cov}(T, U)\\)."),
                 h5("Plugin Type I Error:"),
                 p("$$\\text{Type I Error}_{\\text{Exp}} = \\Phi\\!\\left(-z_{1-\\alpha} 
                     - \\text{Bias}(\\log\\hat\\lambda) \\cdot \\sqrt{E[D]}\\right)$$"),
                 h5("Expected events with staggered entry (Unif\\((0, E_{\\max})\\)):"),
                 p("$$E[D] = n\\left(1 - \\frac{e^{-\\lambda A}(e^{\\lambda E_{\\max}}-1)}{\\lambda E_{\\max}}\\right)$$"),
                 p("where \\(A\\) = admin censor time, \\(E_{\\max}\\) = max entry time.")
               ),
               
               # --- Two-Arm Cox ---
               h4("7. Two-Arm Cox / Log-Rank Test"),
               wellPanel(
                 h5("Four-step delta method chain:"),
                 p("$$\\text{Bias}(\\bar{T}) \\xrightarrow{\\times(-\\lambda_0^2)} 
                     \\text{Bias}(\\hat\\lambda) \\xrightarrow{\\div\\lambda_0}
                     \\text{Bias}(\\log\\text{HR}) \\xrightarrow{\\times\\sqrt{D_{\\text{tot}}/4}}
                     \\text{Bias}(Z_{\\text{cox}})$$"),
                 h5("Type I Error:"),
                 p("$$\\text{Type I Error}_{\\text{Cox}} = 
                     \\Phi\\!\\left(-z_{1-\\alpha} - \\text{Bias}(Z_{\\text{cox}})\\right)$$"),
                 p("Since \\(\\text{Bias}(Z_{\\text{cox}}) < 0\\) (treatment appears better), 
                    the argument exceeds \\(-z_{1-\\alpha}\\), confirming inflation."),
                 h5("Why the event-count approach fails:"),
                 p("Error 1: Landmark-to-event conversion depends on the full censoring distribution, 
                    not just \\(S_0(\\tau)\\)."),
                 p("Error 2: In a two-arm test, perturbing \\(D_{\\text{trt}}\\) by \\(\\Delta\\) shifts 
                    \\(O_{\\text{trt}} - E_{\\text{trt}}\\) by only \\(\\Delta/2\\), not \\(\\Delta\\).")
               ),
               
               # --- Inflation Hierarchy ---
               h4("8. Type I Error Inflation Hierarchy"),
               wellPanel(
                 p("Binary Z-test (\\(\\sim 1.9\\times\\alpha\\)) > 
                    Landmark Z-test (\\(\\sim 1.2\\times\\)) > 
                    One-arm exponential (\\(\\sim 1.1\\text{--}1.3\\times\\)) > 
                    Two-arm Cox (\\(\\sim 1.1\\times\\))"),
                 p("Follows from: \\(|\\text{Corr}(X,U)| > |\\text{Corr}(S(\\tau),U)| > |\\text{Corr}(T,U)|\\) 
                    under imperfect surrogacy.")
               ),
               
               # --- Copula Note ---
               h4("9. Copula Correlation Attenuation"),
               wellPanel(
                 p("The Gaussian copula parameter \\(\\rho\\) is the ", em("latent"), " correlation. 
                    The observed point-biserial between Bernoulli(\\(p\\)) and Exp(\\(\\lambda\\)) 
                    is always attenuated:"),
                 p("$$r_{\\text{obs}} < \\rho \\quad \\text{(typically 70--80\\% of } \\rho \\text{)}$$"),
                 p("Sources: discretization of \\(Z_1\\) into binary; nonlinear quantile transform 
                    \\(T = -\\log(1-\\Phi(Z_2))/\\lambda\\) (Genest & Neslehova 2007; Demirtas & Hedeker 2011).")
               ),
               
               hr(),
               p(em("Formulas from manuscript Sections 2.5--2.7 and Appendices B--D. 
                     Key references: Bauer et al. (2010), Rosenbaum (1961), Prentice (1989)."))
             )
    )
  )
)


# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(simulation = NULL)
  
  # ──────────────────────────────────────────────────────────────
  # Helper: safe column accessor
  # ──────────────────────────────────────────────────────────────
  safe_get <- function(s, name) {
    if (name %in% names(s)) as.numeric(s[[name]]) else NA_real_
  }
  
  # ──────────────────────────────────────────────────────────────
  # Utility Score Design (unchanged)
  # ──────────────────────────────────────────────────────────────
  observeEvent(input$compute_util, {
    tryCatch({
      u_vec <- if (input$mode != "margin-based" && nchar(input$u) > 0) {
        as.numeric(unlist(strsplit(input$u, ",")))
      } else NULL
      
      params <- list(
        pL = if (input$mode %in% c("margin-based", "partial-direct")) input$pL else NULL,
        qL = if (input$mode %in% c("margin-based", "partial-direct")) input$qL else NULL,
        delta = if (input$mode == "margin-based") input$delta else NULL,
        d = if (input$mode == "margin-based") input$d else NULL,
        pL_L = if (input$mode == "full-direct") input$pL_L else NULL,
        qL_L = if (input$mode == "full-direct") input$qL_L else NULL,
        pH_L = if (input$mode == "full-direct") input$pH_L else if (input$mode == "partial-direct") input$pH_L else NULL,
        qH_L = if (input$mode == "full-direct") input$qH_L else if (input$mode == "partial-direct") input$qH_L else NULL,
        pL_H = if (input$mode == "full-direct") input$pL_H else if (input$mode == "partial-direct") input$pL_H else NULL,
        qL_H = if (input$mode == "full-direct") input$qL_H else if (input$mode == "partial-direct") input$qL_H else NULL,
        pH_H = if (input$mode == "full-direct") input$pH_H else NULL,
        qH_H = if (input$mode == "full-direct") input$qH_H else NULL,
        phi = input$phi,
        alpha_L = input$alpha_L, alpha_H = input$alpha_H,
        u = u_vec
      )
      
      if (input$method == "approx") {
        res <- do.call(calc_sample_size_utility_approx, params)
      } else {
        params$max_n <- input$max_n
        params$buffer <- input$buffer
        params$den <- input$den
        res <- do.call(calc_sample_size_utility_exact, params)
      }
      
      output$util_result <- renderPrint({ print(res) })
      
    }, error = function(e) {
      output$util_result <- renderPrint({ cat("Error: ", e$message) })
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # ──────────────────────────────────────────────────────────────
  # ROSE Design (unchanged)
  # ──────────────────────────────────────────────────────────────
  observeEvent(input$compute_rose, {
    tryCatch({
      if (input$method_rose == "approx") {
        res <- calc_sample_size_rose_approx(
          pL = input$pL_rose, delta = input$delta_rose,
          alpha_L = input$alpha_L_rose, alpha_H = input$alpha_H_rose
        )
      } else {
        res <- calc_sample_size_rose_exact(
          pL = input$pL_rose, delta = input$delta_rose,
          alpha_L = input$alpha_L_rose, alpha_H = input$alpha_H_rose,
          max_n = input$max_n_rose, buffer = input$buffer_rose, verbose = FALSE
        )
      }
      output$rose_result <- renderPrint({ print(res) })
    }, error = function(e) {
      output$rose_result <- renderPrint({ cat("Error: ", e$message) })
      showNotification(paste("ROSE Error:", e$message), type = "error")
    })
  })
  
  # ──────────────────────────────────────────────────────────────
  # Simulation Analysis – Main observer
  # ──────────────────────────────────────────────────────────────
  observeEvent(input$run_simulation, {
    withProgress(message = 'Running Simulation...', value = 0, {
      incProgress(0.1, detail = "Gathering parameters...")
      
      p <- input$p_sim; q <- input$q_sim; phi <- input$phi_sim
      u <- tryCatch(as.numeric(unlist(strsplit(input$u_sim, ","))),
                    error = function(e) stop("Invalid utility scores"))
      if (length(u) != 4) stop("Utility must have 4 values")
      
      n1 <- input$n1_sim; n2 <- input$n2_sim
      lambda_u <- input$lambda_u_sim; alpha <- input$alpha_sim
      nSim <- input$nSim; seed <- input$seed
      perform_tte <- input$perform_tte_sim
      
      incProgress(0.1, detail = "Calculating analytical benchmarks...")
      analytical_bias_comp <- calc_analytical_bias(p, q, phi, u, n1, n2, lambda_u)
      analytical_t1e_util <- calc_analytical_type1_error(p, analytical_bias_comp$bias_utility_combined, n1 + n2, alpha)
      analytical_t1e_resp <- calc_analytical_type1_error(p, analytical_bias_comp$bias_response_combined, n1 + n2, alpha)
      
      incProgress(0.3, detail = "Simulating trial stages...")
      sim_results <- compare_bias_methods_fast_enhanced_v2(
        pL = p, pH = p, qL = q, qH = q, phi = phi,
        N1 = n1, N2 = n2, lambda_u = lambda_u, u = u,
        nSim = nSim, ranseed = seed, Alpha = alpha,
        perform_tte_analysis = perform_tte,
        tte_rate = input$tte_rate_sim,
        corr_efficacy_tte = input$corr_eff_tte_sim,
        landmark_time = input$landmark_time_sim,
        entry_time_max = input$entry_time_max_sim,
        admin_censor_time = input$admin_censor_time_sim,
        dropout_rate = input$dropout_rate_sim,
        tte_rate_historical = input$tte_rate_sim,
        Alpha_tte = input$alpha_tte_sim,
        return_raw = FALSE,
        two_arm_survival_flag = TRUE,
        use_internal_parallel = TRUE
      )
      
      incProgress(0.9, detail = "Finalizing results...")
      
      rv$simulation <- list(
        sim = sim_results$Summary_Results,
        analytical = list(
          bias_utility = analytical_bias_comp$bias_utility_combined,
          bias_response = analytical_bias_comp$bias_response_combined,
          bias_U_score = analytical_bias_comp$bias_U_score_combined,
          t1e_z_utility = analytical_t1e_util$type1_z,
          t1e_z_response = analytical_t1e_resp$type1_z,
          t1e_bin_utility = analytical_t1e_util$type1_bin,
          t1e_bin_response = analytical_t1e_resp$type1_bin
        ),
        was_tte_run = perform_tte
      )
      
      cat("Simulation completed.\n")
    })
  })
  
  # Show/hide TTE tab
  observe({
    if (!is.null(rv$simulation) && rv$simulation$was_tte_run) {
      showTab("results_tabs", "TTE Endpoint Results")
    } else {
      hideTab("results_tabs", "TTE Endpoint Results")
    }
  })
  
  # ================================================================
  # BINARY ENDPOINT — Bias Table (long format, Observed appears ONCE)
  # ================================================================
  output$binary_bias_table <- renderDT({
    req(rv$simulation)
    res <- rv$simulation; s <- res$sim
    
    df <- data.frame(
      Method = c(
        "Observed (simulation mean)",
        "Analytical (utility-based, Cov(X,U))",
        "Analytical (max bound, response-only)",
        "Plugin mean (utility-based, Cov(X,U))",
        "Plugin mean (max bound, response-only)",
        "Utility score bias (observed, for U itself)",
        "Utility score bias (analytical, for U itself)",
        "Utility score bias (plugin, for U itself)"
      ),
      Bias = c(
        safe_get(s, "mean_bias_observed_selected_combined"),
        res$analytical$bias_utility,
        res$analytical$bias_response,
        safe_get(s, "mean_bias_utility_plugin_selected_combined"),
        safe_get(s, "mean_bias_response_plugin_selected_combined"),
        safe_get(s, "mean_bias_utility_score_observed_combined"),
        res$analytical$bias_U_score,
        safe_get(s, "mean_bias_utility_score_plugin_combined")
      )
    )
    
    datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 10),
              caption = "Response Rate Bias & Utility Score Bias: Observed is the ground truth; utility-based and max-bound are predictions.") %>%
      formatRound(columns = 2, digits = 6)
  })
  
  # ================================================================
  # BINARY ENDPOINT — Type I Error Table (Observed appears ONCE per test)
  # ================================================================
  output$binary_type1_table <- renderDT({
    req(rv$simulation)
    res <- rv$simulation; s <- res$sim
    alpha_val <- input$alpha_sim
    nSim <- input$nSim
    
    # Helper: MC SE for observed (Bernoulli) = sqrt(p*(1-p)/n)
    mc_se_bern <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)
    # Helper: MC SE for plugin mean = sd / sqrt(n)
    mc_se_mean <- function(sd_key) safe_get(s, sd_key) / sqrt(nSim)
    
    t1e_vals <- c(
      # Z-test
      safe_get(s, "mean_reject_H0_z"),
      res$analytical$t1e_z_utility,
      res$analytical$t1e_z_response,
      safe_get(s, "mean_type1_z_utility_plugin"),
      safe_get(s, "mean_type1_z_response_plugin"),
      # Binomial
      safe_get(s, "mean_reject_H0_bin"),
      res$analytical$t1e_bin_utility,
      res$analytical$t1e_bin_response,
      safe_get(s, "mean_type1_bin_utility_plugin"),
      safe_get(s, "mean_type1_bin_response_plugin")
    )
    
    # MC SE: Bernoulli for observed, 0 for analytical (exact), sd/sqrt(n) for plugin
    mc_se_vals <- c(
      mc_se_bern(t1e_vals[1], nSim),
      NA,  # analytical (exact formula, no MC error)
      NA,
      mc_se_mean("sd_type1_z_utility_plugin"),
      mc_se_mean("sd_type1_z_response_plugin"),
      mc_se_bern(t1e_vals[6], nSim),
      NA,
      NA,
      mc_se_mean("sd_type1_bin_utility_plugin"),
      mc_se_mean("sd_type1_bin_response_plugin")
    )
    
    df <- data.frame(
      Test = c(
        rep("Z-test", 5),
        rep("Binomial", 5)
      ),
      Method = rep(c(
        "Observed (simulation)",
        "Analytical (utility-based)",
        "Analytical (max bound)",
        "Plugin mean (utility-based)",
        "Plugin mean (max bound)"
      ), 2),
      Type_I_Error = t1e_vals,
      MC_SE = mc_se_vals
    )
    
    df$Inflation <- round(df$Type_I_Error / alpha_val, 3)
    
    # Caption with conditional MC warning
    cap <- paste0(
      "Binary Type I Error (nominal alpha = ", alpha_val, "). ",
      "Expected ordering per test: Utility-based <= Observed <= Max bound."
    )
    if (nSim < 10000) {
      cap <- paste0(
        cap,
        " WARNING: nSim = ", nSim, " is low; MC SE (~",
        round(mc_se_bern(0.025, nSim), 4),
        " at alpha) may obscure the true ordering. ",
        "Consider nSim >= 10,000 for reliable comparisons."
      )
    }
    
    datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 12),
              caption = cap) %>%
      formatRound(columns = c("Type_I_Error", "MC_SE"), digits = 5) %>%
      formatRound(columns = "Inflation", digits = 3)
  })
  
  # ================================================================
  # BINARY — Download
  # ================================================================
  output$download_binary <- downloadHandler(
    filename = "binary_results.csv",
    content = function(file) {
      res <- rv$simulation; s <- res$sim
      alpha_val <- input$alpha_sim
      nSim <- input$nSim
      mc_se_bern <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)
      mc_se_mean <- function(sd_key) safe_get(s, sd_key) / sqrt(nSim)
      
      df <- data.frame(
        Test = c(rep("Bias", 5), rep("T1E Z-test", 5), rep("T1E Binomial", 5)),
        Method = rep(c("Observed", "Analytical (utility)", "Analytical (max)",
                       "Plugin (utility)", "Plugin (max)"), 3),
        Value = c(
          safe_get(s, "mean_bias_observed_selected_combined"),
          res$analytical$bias_utility, res$analytical$bias_response,
          safe_get(s, "mean_bias_utility_plugin_selected_combined"),
          safe_get(s, "mean_bias_response_plugin_selected_combined"),
          safe_get(s, "mean_reject_H0_z"),
          res$analytical$t1e_z_utility, res$analytical$t1e_z_response,
          safe_get(s, "mean_type1_z_utility_plugin"),
          safe_get(s, "mean_type1_z_response_plugin"),
          safe_get(s, "mean_reject_H0_bin"),
          res$analytical$t1e_bin_utility, res$analytical$t1e_bin_response,
          safe_get(s, "mean_type1_bin_utility_plugin"),
          safe_get(s, "mean_type1_bin_response_plugin")
        ),
        MC_SE = c(
          # Bias
          mc_se_mean("sd_bias_observed_selected_combined"),
          NA, NA,  # analytical (exact)
          mc_se_mean("sd_bias_utility_plugin_selected_combined"),
          mc_se_mean("sd_bias_response_plugin_selected_combined"),
          # T1E Z-test
          mc_se_bern(safe_get(s, "mean_reject_H0_z"), nSim),
          NA, NA,
          mc_se_mean("sd_type1_z_utility_plugin"),
          mc_se_mean("sd_type1_z_response_plugin"),
          # T1E Binomial
          mc_se_bern(safe_get(s, "mean_reject_H0_bin"), nSim),
          NA, NA,
          mc_se_mean("sd_type1_bin_utility_plugin"),
          mc_se_mean("sd_type1_bin_response_plugin")
        )
      )
      write.csv(df, file, row.names = FALSE, na = "")
    }
  )
  
  # ================================================================
  # BINARY — Bias bar plot
  # ================================================================
  output$bias_comparison_plot_sim <- renderPlot({
    req(rv$simulation)
    res <- rv$simulation; s <- res$sim
    
    vals <- c(
      safe_get(s, "mean_bias_observed_selected_combined"),
      res$analytical$bias_utility,
      res$analytical$bias_response,
      safe_get(s, "mean_bias_utility_plugin_selected_combined"),
      safe_get(s, "mean_bias_response_plugin_selected_combined")
    )
    nms <- c("Observed", "Anal.\n(utility)", "Anal.\n(max)", "Plugin\n(utility)", "Plugin\n(max)")
    cols <- c("#2ca02c", "#1f77b4", "#ff7f0e", "#1f77b4", "#ff7f0e")
    border_cols <- c("black", NA, NA, NA, NA)
    
    bp <- barplot(vals, names.arg = nms,
                  main = "Response Rate Bias Comparison",
                  ylab = "Bias", col = cols, border = border_cols)
    abline(h = 0, lty = 2, col = "gray40")
    
    # Mark observed with a heavier border
    rect(bp[1] - 0.5, 0, bp[1] + 0.5, vals[1], border = "black", lwd = 2)
  })
  
  # ================================================================
  # BINARY — Type I Error bar plot
  # ================================================================
  output$binary_type1_plot_sim <- renderPlot({
    req(rv$simulation)
    res <- rv$simulation; s <- res$sim
    alpha_val <- input$alpha_sim
    
    # Z-test values
    z_vals <- c(
      safe_get(s, "mean_reject_H0_z"),
      res$analytical$t1e_z_utility,
      res$analytical$t1e_z_response,
      safe_get(s, "mean_type1_z_utility_plugin"),
      safe_get(s, "mean_type1_z_response_plugin")
    )
    # Binomial values
    bin_vals <- c(
      safe_get(s, "mean_reject_H0_bin"),
      res$analytical$t1e_bin_utility,
      res$analytical$t1e_bin_response,
      safe_get(s, "mean_type1_bin_utility_plugin"),
      safe_get(s, "mean_type1_bin_response_plugin")
    )
    
    plot_mat <- rbind(`Z-test` = z_vals, Binomial = bin_vals)
    colnames(plot_mat) <- c("Observed", "Anal.(util)", "Anal.(max)", "Plugin(util)", "Plugin(max)")
    
    barplot(plot_mat, beside = TRUE,
            main = "Binary Type I Error by Method",
            ylab = "Type I Error Rate",
            col = c("#1f77b4", "#ff7f0e"),
            legend.text = TRUE,
            args.legend = list(x = "topleft", bty = "n"))
    abline(h = alpha_val, lty = 2, col = "red", lwd = 2)
    text(0.5, alpha_val + 0.003, paste0("Nominal alpha = ", alpha_val),
         adj = 0, col = "red", cex = 0.9)
  })
  
  # ================================================================
  # TTE — Landmark Bias Table
  # ================================================================
  output$tte_bias_table <- renderDT({
    req(rv$simulation, rv$simulation$was_tte_run)
    s <- rv$simulation$sim
    
    df <- data.frame(
      Method = c("Observed (simulation mean)",
                 "Plugin (direct Cov(S(tau),U))",
                 "Surrogate-adjusted",
                 "Upper bound (response-based)",
                 "Upper bound (utility-based)"),
      Landmark_Bias_S_tau = c(
        safe_get(s, "mean_bias_landmark_observed_combined"),
        safe_get(s, "mean_bias_landmark_plugin_combined"),
        safe_get(s, "mean_bias_landmark_surrogate_adjusted_combined"),
        safe_get(s, "mean_bias_landmark_upper_bound_combined"),
        safe_get(s, "mean_bias_landmark_utility_upper_combined")
      ),
      Mean_TTE_Bias_T_bar = c(
        NA_real_,
        safe_get(s, "mean_bias_TTE_plugin"),
        NA_real_,
        NA_real_,
        NA_real_
      )
    )
    datatable(df, rownames = FALSE, options = list(dom = 't'),
              caption = paste0(
                "Landmark S(tau) Bias. Expected ordering: ",
                "Plugin <= Observed <= Upper bound(s). ",
                "Surrogate-adjusted anchors on utility plugin bias via rho_{X,S(tau)}."
              )) %>%
      formatRound(columns = 2:3, digits = 6)
  })
  
  # ================================================================
  # TTE — Type I Error Table (Observed ONCE per test)
  # ================================================================
  output$tte_type1_table <- renderDT({
    req(rv$simulation, rv$simulation$was_tte_run)
    s <- rv$simulation$sim
    alpha_tte <- input$alpha_tte_sim
    nSim <- input$nSim
    
    # Helper: MC SE for observed (Bernoulli) = sqrt(p*(1-p)/n)
    mc_se_bern <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)
    # Helper: MC SE for plugin mean = sd / sqrt(n)
    mc_se_mean <- function(sd_key) safe_get(s, sd_key) / sqrt(nSim)
    
    # Extract Type I error values
    t1e_vals <- c(
      safe_get(s, "mean_reject_H0_landmark_z"),
      safe_get(s, "mean_type1_landmark_z_plugin"),
      safe_get(s, "mean_type1_landmark_z_surrogate"),
      safe_get(s, "mean_type1_landmark_z_upper_bound"),
      safe_get(s, "mean_type1_landmark_z_utility_upper"),
      safe_get(s, "mean_reject_H0_exp_one_arm"),
      safe_get(s, "mean_type1_exp_one_arm_plugin"),
      safe_get(s, "mean_type1_exp_one_arm_upper"),
      safe_get(s, "mean_type1_exp_one_arm_utility_upper"),
      safe_get(s, "mean_reject_H0_logrank_two_arm"),
      safe_get(s, "mean_type1_logrank_two_arm_plugin"),
      safe_get(s, "mean_type1_logrank_two_arm_upper"),
      safe_get(s, "mean_type1_logrank_two_arm_utility_upper")
    )
    
    # MC SE: Bernoulli for observed rows, sd/sqrt(n) for plugin/analytical rows
    mc_se_vals <- c(
      mc_se_bern(t1e_vals[1], nSim),
      mc_se_mean("sd_type1_landmark_z_plugin"),
      mc_se_mean("sd_type1_landmark_z_surrogate"),
      mc_se_mean("sd_type1_landmark_z_upper_bound"),
      mc_se_mean("sd_type1_landmark_z_utility_upper"),
      mc_se_bern(t1e_vals[6], nSim),
      mc_se_mean("sd_type1_exp_one_arm_plugin"),
      mc_se_mean("sd_type1_exp_one_arm_upper"),
      mc_se_mean("sd_type1_exp_one_arm_utility_upper"),
      mc_se_bern(t1e_vals[10], nSim),
      mc_se_mean("sd_type1_logrank_two_arm_plugin"),
      mc_se_mean("sd_type1_logrank_two_arm_upper"),
      mc_se_mean("sd_type1_logrank_two_arm_utility_upper")
    )
    
    df <- data.frame(
      Test = c(
        rep("Landmark Z-test", 5),
        rep("One-Arm Exp (log-scale)", 4),
        rep("Two-Arm Cox/Log-Rank", 4)
      ),
      Method = c(
        "Observed (simulation)", "Plugin", "Surrogate-adj", "Upper (resp-based)", "Upper (utility-based)",
        "Observed (simulation)", "Plugin", "Upper (resp-based)", "Upper (utility-based)",
        "Observed (simulation)", "Plugin", "Upper (resp-based)", "Upper (utility-based)"
      ),
      Type_I_Error = t1e_vals,
      MC_SE = mc_se_vals
    )
    
    df$Inflation <- round(df$Type_I_Error / alpha_tte, 3)
    
    # Caption with conditional MC warning
    cap <- paste0(
      "TTE Type I Error (nominal alpha = ", alpha_tte, "). ",
      "Expected ordering per test: Plugin <= Observed <= Upper bound."
    )
    if (nSim < 10000) {
      cap <- paste0(
        cap,
        " WARNING: nSim = ", nSim, " is low; MC SE (~",
        round(mc_se_bern(0.025, nSim), 4),
        " at alpha) may obscure the true Plugin <= Observed ordering. ",
        "Consider nSim >= 10,000 for reliable comparisons."
      )
    }
    
    datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 15),
              caption = cap) %>%
      formatRound(columns = c("Type_I_Error", "MC_SE"), digits = 5) %>%
      formatRound(columns = "Inflation", digits = 3)
  })
  
  # ================================================================
  # TTE — Correlation Diagnostics (unchanged)
  # ================================================================
  output$tte_corr_table <- renderDT({
    req(rv$simulation, rv$simulation$was_tte_run)
    s <- rv$simulation$sim
    
    df <- data.frame(
      Quantity = c(
        "Latent Gaussian copula rho (input)",
        "Observed Corr(Efficacy, TTE) [point-biserial]",
        "Observed Corr(Response, Landmark S(tau))",
        "Observed Corr(TTE, Utility)",
        "Observed phi(Efficacy, Safety)",
        "Mean Events (Stage 1)",
        "Mean Events (Stage 2)",
        "Mean HR (Two-Arm)"
      ),
      Value = c(
        input$corr_eff_tte_sim,
        safe_get(s, "mean_corr_tte_est"),
        safe_get(s, "mean_corr_response_landmark_plugin"),
        safe_get(s, "mean_corr_TTE_U_plugin"),
        safe_get(s, "mean_phi_est"),
        safe_get(s, "mean_events_stage1"),
        safe_get(s, "mean_events_stage2"),
        safe_get(s, "mean_hr")
      )
    )
    datatable(df, rownames = FALSE, options = list(dom = 't'),
              caption = "Correlation & Diagnostic Summaries") %>%
      formatRound(columns = 2, digits = 4)
  })
  
  # ================================================================
  # TTE — Download (unchanged)
  # ================================================================
  output$download_tte_full <- downloadHandler(
    filename = "tte_full_results.csv",
    content = function(file) {
      s <- rv$simulation$sim
      nSim <- input$nSim
      mc_se_bern <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)
      mc_se_mean <- function(sd_key) safe_get(s, sd_key) / sqrt(nSim)
      
      bias_df <- data.frame(
        Category = "Bias",
        Endpoint = c("Landmark S(tau)", "Landmark S(tau)", "Landmark S(tau)", 
                      "Landmark S(tau)", "Landmark S(tau)", "Mean TTE"),
        Method = c("Observed", "Plugin", "Surrogate-adj", "Upper (resp-based)", 
                   "Upper (utility-based)", "Plugin"),
        Value = c(
          safe_get(s, "mean_bias_landmark_observed_combined"),
          safe_get(s, "mean_bias_landmark_plugin_combined"),
          safe_get(s, "mean_bias_landmark_surrogate_adjusted_combined"),
          safe_get(s, "mean_bias_landmark_upper_bound_combined"),
          safe_get(s, "mean_bias_landmark_utility_upper_combined"),
          safe_get(s, "mean_bias_TTE_plugin")
        ),
        MC_SE = c(
          mc_se_mean("sd_bias_landmark_observed_combined"),
          mc_se_mean("sd_bias_landmark_plugin_combined"),
          mc_se_mean("sd_bias_landmark_surrogate_adjusted_combined"),
          mc_se_mean("sd_bias_landmark_upper_bound_combined"),
          mc_se_mean("sd_bias_landmark_utility_upper_combined"),
          mc_se_mean("sd_bias_TTE_plugin")
        )
      )
      
      t1e_vals <- c(
        safe_get(s, "mean_reject_H0_landmark_z"),
        safe_get(s, "mean_type1_landmark_z_plugin"),
        safe_get(s, "mean_type1_landmark_z_surrogate"),
        safe_get(s, "mean_type1_landmark_z_upper_bound"),
        safe_get(s, "mean_type1_landmark_z_utility_upper"),
        safe_get(s, "mean_reject_H0_exp_one_arm"),
        safe_get(s, "mean_type1_exp_one_arm_plugin"),
        safe_get(s, "mean_type1_exp_one_arm_upper"),
        safe_get(s, "mean_type1_exp_one_arm_utility_upper"),
        safe_get(s, "mean_reject_H0_logrank_two_arm"),
        safe_get(s, "mean_type1_logrank_two_arm_plugin"),
        safe_get(s, "mean_type1_logrank_two_arm_upper"),
        safe_get(s, "mean_type1_logrank_two_arm_utility_upper")
      )
      
      t1e_df <- data.frame(
        Category = "Type I Error",
        Endpoint = c(rep("Landmark Z-test", 5), rep("One-Arm Exp", 4), rep("Two-Arm Cox", 4)),
        Method = c("Observed", "Plugin", "Surrogate-adj", "Upper (resp)", "Upper (util)",
                   "Observed", "Plugin", "Upper (resp)", "Upper (util)",
                   "Observed", "Plugin", "Upper (resp)", "Upper (util)"),
        Value = t1e_vals,
        MC_SE = c(
          mc_se_bern(t1e_vals[1], nSim),
          mc_se_mean("sd_type1_landmark_z_plugin"),
          mc_se_mean("sd_type1_landmark_z_surrogate"),
          mc_se_mean("sd_type1_landmark_z_upper_bound"),
          mc_se_mean("sd_type1_landmark_z_utility_upper"),
          mc_se_bern(t1e_vals[6], nSim),
          mc_se_mean("sd_type1_exp_one_arm_plugin"),
          mc_se_mean("sd_type1_exp_one_arm_upper"),
          mc_se_mean("sd_type1_exp_one_arm_utility_upper"),
          mc_se_bern(t1e_vals[10], nSim),
          mc_se_mean("sd_type1_logrank_two_arm_plugin"),
          mc_se_mean("sd_type1_logrank_two_arm_upper"),
          mc_se_mean("sd_type1_logrank_two_arm_utility_upper")
        )
      )
      
      corr_df <- data.frame(
        Category = "Correlation",
        Endpoint = c("Eff-TTE", "Resp-Landmark", "TTE-Utility", "Eff-Safety"),
        Method = rep("Observed", 4),
        Value = c(
          safe_get(s, "mean_corr_tte_est"),
          safe_get(s, "mean_corr_response_landmark_plugin"),
          safe_get(s, "mean_corr_TTE_U_plugin"),
          safe_get(s, "mean_phi_est")
        ),
        MC_SE = c(
          mc_se_mean("sd_corr_tte_est"),
          mc_se_mean("sd_corr_response_landmark_plugin"),
          mc_se_mean("sd_corr_TTE_U_plugin"),
          mc_se_mean("sd_phi_est")
        )
      )
      
      write.csv(rbind(bias_df, t1e_df, corr_df), file, row.names = FALSE, na = "")
    }
  )
  
  # ================================================================
  # TTE — Bias Plot (unchanged)
  # ================================================================
  output$tte_bias_plot_sim <- renderPlot({
    req(rv$simulation, rv$simulation$was_tte_run)
    s <- rv$simulation$sim
    
    values <- c(
      safe_get(s, "mean_bias_landmark_observed_combined"),
      safe_get(s, "mean_bias_landmark_plugin_combined"),
      safe_get(s, "mean_bias_landmark_surrogate_adjusted_combined"),
      safe_get(s, "mean_bias_landmark_upper_bound_combined"),
      safe_get(s, "mean_bias_landmark_utility_upper_combined")
    )
    nms <- c("Observed", "Plugin", "Surrogate-adj", "Upper (resp)", "Upper (util)")
    cols <- c("#2ca02c", "#1f77b4", "#9467bd", "#ff7f0e", "#d62728")
    
    barplot(values, names.arg = nms,
            main = "Landmark Survival Bias: All Estimation Methods",
            ylab = "Bias(S(tau))", col = cols, border = NA, las = 2)
    abline(h = 0, lty = 2, col = "gray40")
  })
  
  # ================================================================
  # TTE — Type I Error Plot
  # ================================================================
  output$tte_type1_plot_sim <- renderPlot({
    req(rv$simulation, rv$simulation$was_tte_run)
    s <- rv$simulation$sim
    alpha_tte <- input$alpha_tte_sim
    
    obs_vals <- c(
      safe_get(s, "mean_reject_H0_landmark_z"),
      safe_get(s, "mean_reject_H0_exp_one_arm"),
      safe_get(s, "mean_reject_H0_logrank_two_arm")
    )
    plugin_vals <- c(
      safe_get(s, "mean_type1_landmark_z_plugin"),
      safe_get(s, "mean_type1_exp_one_arm_plugin"),
      safe_get(s, "mean_type1_logrank_two_arm_plugin")
    )
    resp_upper_vals <- c(
      safe_get(s, "mean_type1_landmark_z_upper_bound"),
      safe_get(s, "mean_type1_exp_one_arm_upper"),
      safe_get(s, "mean_type1_logrank_two_arm_upper")
    )
    util_upper_vals <- c(
      safe_get(s, "mean_type1_landmark_z_utility_upper"),
      safe_get(s, "mean_type1_exp_one_arm_utility_upper"),
      safe_get(s, "mean_type1_logrank_two_arm_utility_upper")
    )
    
    plot_mat <- rbind(Observed = obs_vals, Plugin = plugin_vals, 
                      `Upper(resp)` = resp_upper_vals, `Upper(util)` = util_upper_vals)
    colnames(plot_mat) <- c("Landmark Z", "One-Arm Exp\n(log-scale)", "Two-Arm Cox")
    
    barplot(plot_mat, beside = TRUE,
            main = "TTE Type I Error: Observed vs Plugin vs Upper Bounds",
            ylab = "Type I Error Rate",
            col = c("#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"),
            legend.text = TRUE,
            args.legend = list(x = "topright", bty = "n"))
    abline(h = alpha_tte, lty = 2, col = "red", lwd = 2)
    text(0.5, alpha_tte + 0.002, paste0("Nominal alpha = ", alpha_tte),
         adj = 0, col = "red", cex = 0.9)
  })
}

# Run the application
shinyApp(ui = ui, server = server)