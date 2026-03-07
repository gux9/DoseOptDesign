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
  
  # ========================================================================
  # ANALYTICAL UTILITY SCORE BIAS (for U itself): Bias(Ū) = σ_U/√(n₁π) × exp(...)
  # Since Cov(U,U) = Var(U) = σ_U², the general formula gives σ_U/√(n₁π) × exp(...)
  # ========================================================================
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
  # Plugin: uses Cov(X,U)-based utility bias from design parameters
  type1_z_utility_analytical_combined <- 1 - pnorm(z_crit - bias_utility_analytical_combined / se_null)
  # Upper: uses analytical utility bias (design-based), replaces response-only max
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
    
    # Observed utility score bias (U itself)
    bias_utility_score_observed_stage1   = NA_real_,
    bias_utility_score_observed_combined = NA_real_,
    
    # Plugin bias estimates - Null case (selected dose method)
    bias_utility_plugin_selected_stage1   = NA_real_,
    bias_utility_plugin_selected_combined = NA_real_,
    bias_response_plugin_selected_stage1  = NA_real_,
    bias_response_plugin_selected_combined= NA_real_,
    
    # Plugin utility score bias (for U itself): σ̂_U/√(n₁π) × exp(...)
    bias_utility_score_plugin_stage1   = NA_real_,
    bias_utility_score_plugin_combined = NA_real_,
    
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
    type1_z_utility_plugin  = NA_real_,    # Cov(X,U)-based from Stage 1
    type1_bin_utility_plugin= NA_real_,    # Cov(X,U)-based from Stage 1
    type1_z_response_plugin = NA_real_,    # Response-only max bound: √(p(1-p)/nπ)
    type1_bin_response_plugin= NA_real_,   # Response-only max bound
    type1_z_utility_upper   = NA_real_,    # Analytical utility-based bias from design params
    type1_bin_utility_upper = NA_real_,    # Analytical utility-based bias from design params
    
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
      
      # Upper bound bias (Method 1: response-based, Section 2.7.4)
      bias_landmark_upper_bound_combined = NA_real_,
      # Upper bound bias (utility-score-based: σ_S/σ_U × Bias(Ū))
      bias_landmark_utility_upper_combined = NA_real_,
      
      # Landmark survival Z-test (Equation 31)
      p_value_landmark_z_test = NA_real_,
      reject_H0_landmark_z    = NA_real_,
      
      # Plugin Type I error for landmark test
      type1_landmark_z_plugin = NA_real_,
      type1_landmark_z_upper_bound = NA_real_,      # response-based upper
      type1_landmark_z_utility_upper = NA_real_,     # utility-score-based upper
      type1_landmark_z_surrogate = NA_real_,
      
      # One-sample exponential test
      p_value_exponential = NA_real_,
      reject_H0_exp_one_arm = NA_real_,
      
      # Plugin Type I error for one-arm exponential test
      # Uses direct Cov(TTE, U) estimation from Stage 1 data
      type1_exp_one_arm_plugin = NA_real_,
      type1_exp_one_arm_upper = NA_real_,         # response-based upper
      type1_exp_one_arm_utility_upper = NA_real_,  # utility-score-based upper
      
      # Direct TTE bias estimation from Stage 1 (Cov(TTE, U))
      cov_TTE_U_plugin = NA_real_,
      corr_TTE_U_plugin = NA_real_,
      bias_TTE_plugin = NA_real_,
      bias_TTE_from_corr = NA_real_,
      
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
      type1_logrank_two_arm_upper = NA_real_,          # response-based upper
      type1_logrank_two_arm_utility_upper = NA_real_,   # utility-score-based upper
      
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
  
  # ========================================================================
  # OBSERVED UTILITY SCORE BIAS (for U itself)
  # Under null: true mean utility is E_U_L = E_U_H
  # ========================================================================
  
  # Compute utility for selected dose in stage 1 and combined
  utility_selected_stage1 <- compute_utility(selected_data_efficacy, selected_data_safety, u)
  utility_stage2 <- compute_utility(stage2_data_efficacy, stage2_data_safety, u)
  utility_selected_combined <- (N1 * utility_selected_stage1 + N2 * utility_stage2) / n_total
  
  # True mean utility (same for both doses under null)
  mu_U_true <- E_U_L  # = E_U_H under null
  
  sim_results$bias_utility_score_observed_stage1 <- utility_selected_stage1 - mu_U_true
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
  
  # ========================================================================
  # PLUGIN UTILITY SCORE BIAS (for U itself): σ̂_U/√(n₁π) × exp(...)
  # Cov(U,U) = Var(U) = σ_U², so Cov(U,U)/σ_U = σ_U
  # ========================================================================
  bias_U_score_plugin_stage1 <- sqrt(mom_sel$sigma2) / sqrt(N1 * base::pi) * exp_term_plugin
  bias_U_score_plugin_combined <- (N1 / n_total) * bias_U_score_plugin_stage1
  
  sim_results$bias_utility_score_plugin_stage1 <- bias_U_score_plugin_stage1
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
  # Plugin: uses Cov(X,U)-based bias estimated from Stage 1
  z_shift_util <- bias_util_plugin_combined / se_null
  sim_results$type1_z_utility_plugin <- 1 - pnorm(z_crit - z_shift_util)
  
  # Response-only max upper bound: sqrt(p(1-p))/sqrt(nπ)
  z_shift_resp <- bias_resp_plugin_combined / se_null
  sim_results$type1_z_response_plugin <- 1 - pnorm(z_crit - z_shift_resp)
  
  # Utility-score-based upper: analytical Cov(X,U)/(σ_U√n₁) from design parameters
  z_shift_util_upper <- bias_utility_analytical_combined / se_null
  sim_results$type1_z_utility_upper <- 1 - pnorm(z_crit - z_shift_util_upper)
  
  # FIXED: Binomial test Type I error = P(X >= k_c | p = p0 + bias)
  # = pbinom(k_c - 1, n, p0 + bias, lower.tail = FALSE)
  # Plugin: uses Cov(X,U) from Stage 1
  sim_results$type1_bin_utility_plugin <- pbinom(k_c - 1, size = n_total, 
                                                  prob = p_null + bias_util_plugin_combined, 
                                                  lower.tail = FALSE)
  # Response-only max upper bound
  sim_results$type1_bin_response_plugin <- pbinom(k_c - 1, size = n_total, 
                                                   prob = p_null + bias_resp_plugin_combined, 
                                                   lower.tail = FALSE)
  # Utility-score-based upper bound
  sim_results$type1_bin_utility_upper <- pbinom(k_c - 1, size = n_total, 
                                                 prob = p_null + bias_utility_analytical_combined, 
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
    
    # Method 2: Surrogate-adjusted estimate (Section 2.7.4, Equation 30)
    # Bias(S(tau)) ≈ rho_{XS} × (sigma_S / sigma_X) × Bias_utility(X̂)
    # Uses utility-based response bias plugin (Cov(X,U)/σ_U path)
    bias_landmark_surrogate <- corr_response_landmark * (sigma_S_stage1 / sigma_X_stage1) * 
                               bias_util_plugin_combined
    bias_landmark_surrogate[!is.finite(bias_landmark_surrogate)] <- 0
    sim_results$bias_landmark_surrogate_adjusted_combined <- bias_landmark_surrogate
    
    # Response-based upper bound (assume Corr(S, X) = 1, uses response-only max bias)
    # Bias_max(S(tau)) = (σ_S / σ_X) × √(p(1-p)/n₁π) × (N1/n_total)
    bias_landmark_resp_upper <- (sigma_S_stage1 / sigma_X_stage1) * bias_resp_plugin_combined
    bias_landmark_resp_upper[!is.finite(bias_landmark_resp_upper)] <- 0
    sim_results$bias_landmark_upper_bound_combined <- bias_landmark_resp_upper
    
    # Utility-score-based upper bound (assume Corr(S, U) = σ_S/σ_U)
    # Bias_max(S(tau)) = (σ_S / σ_U) × Bias(Ū) = σ_S / √(n₁π) × (N1/n_total)
    # Accounts for safety-TTE correlation pathway through U
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
    
    # Using utility-score-based upper bound
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
    # Equivalently: Bias(T̄) = Corr(T,U) × (σ_T/σ_U) × Bias(Ū)
    # where Bias(Ū) = σ_U/√(n₁π) × exp(...)
    bias_TTE_plugin_stage1 <- (cov_TTE_U_direct / sqrt(mom_sel$sigma2 * N1)) * 
                               (1 / sqrt(base::pi)) * exp_term_plugin
    bias_TTE_plugin_combined <- (N1 / n_total) * bias_TTE_plugin_stage1
    
    # Verification: Corr(T,U)-based formulation (should equal direct method)
    # Bias(T̄) = Corr(T,U) × σ_T × (1/√(n₁π)) × exp(...)
    bias_TTE_from_corr_stage1 <- corr_TTE_U * sigma_TTE_stage1 / sqrt(N1 * base::pi) * exp_term_plugin
    bias_TTE_from_corr_combined <- (N1 / n_total) * bias_TTE_from_corr_stage1
    
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
    
    # Response-based upper: uses landmark response-based upper bound → hazard via delta method
    bias_log_lambda_resp_upper <- -sim_results$bias_landmark_upper_bound_combined /
      (landmark_time * exp(-tte_rate_historical * landmark_time) * tte_rate_historical)
    
    # Utility-based upper: uses landmark utility-based upper bound → hazard via delta method
    bias_log_lambda_util_upper <- -sim_results$bias_landmark_utility_upper_combined /
      (landmark_time * exp(-tte_rate_historical * landmark_time) * tte_rate_historical)
    
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
    sim_results$bias_TTE_from_corr <- bias_TTE_from_corr_combined
    
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
        # ---- PARALLEL path: standalone / interactive use ----
        n_workers <- max(1, parallelly::availableCores(omit = 2))
        chunk_size <- min(500, ceiling(nSim / n_workers))
        chunk_list <- split(1:nSim, ceiling(seq_along(1:nSim) / chunk_size))
        
        old_max_size <- getOption("future.globals.maxSize")
        options(future.globals.maxSize = 1024 * 1024^2)
        
        cox_results_list <- tryCatch({
          # Attempt cluster setup with retries
          parallel_ok <- FALSE
          for (attempt in 1:3) {
            cluster_ok <- tryCatch({
              plan(multisession, workers = n_workers)
              TRUE
            }, error = function(e) {
              message(sprintf("Cluster setup attempt %d failed: %s", attempt, conditionMessage(e)))
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
              df <- data.frame(
                time   = tte_chunk[, j],
                status = status_chunk[, j],
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
            })
          }, future.seed = TRUE)
          
          plan(sequential)
          res
          
        }, error = function(e) {
          # Fallback: run sequentially
          try(plan(sequential), silent = TRUE)
          message(sprintf("Parallel Cox failed: %s — running sequentially", conditionMessage(e)))
          
          list(lapply(1:nSim, function(i) {
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
            #p_logrank_one_sided = pnorm(beta / se)
            p_logrank_one_sided = pnorm(sign(beta) * sqrt(s$sctest["test"])) #Updated to use logrank p-value, not Wald p-value (21Feb2026)
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
      
      # Response-based upper: use response-based landmark upper → hazard via delta method
      S0_tau <- exp(-tte_rate_historical * landmark_time)
      bias_lambda_trt_resp_upper <- -sim_results$bias_landmark_upper_bound_combined / 
        (landmark_time * S0_tau)
      
      # Utility-based upper: use utility-based landmark upper → hazard via delta method
      bias_lambda_trt_util_upper <- -sim_results$bias_landmark_utility_upper_combined / 
        (landmark_time * S0_tau)
      
      # --- Step 2: Bias in log hazard ratio ---
      # Under null: λ_ctrl = λ₀ (unbiased historical control)
      #             λ_trt  = λ₀ + Bias(λ̂)
      # log(HR) = log(λ_trt / λ_ctrl) ≈ Bias(λ̂) / λ₀  (first-order Taylor)
      
      bias_logHR_plugin     <- bias_lambda_trt_plugin     / tte_rate_historical
      bias_logHR_resp_upper <- bias_lambda_trt_resp_upper  / tte_rate_historical
      bias_logHR_util_upper <- bias_lambda_trt_util_upper  / tte_rate_historical
      
      # --- Step 3: Bias in Cox Z-statistic ---
      # Cox model: β̂ = log(HR), se(β̂) = 1/√I
      # where I ≈ D_total/4 for equal allocation (= observed Fisher information)
      #
      # Estimate D_total: under null with equal allocation, D_total ≈ 2 × D_trt
      # D = colSums(combined_censor_mat) is the treatment arm events
      # D_total_est <- 2 * D
      # V_cox <- pmax(D_total_est / 4, 1)
      # Current: D from Stage 1+2 combined data
      # D_total_est <- 2 * D   # D = colSums(combined_censor_mat)
      
      # Proposed: project from Stage 1 observed event fraction (updated 19Feb2026)
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
    
    # Utility score (U itself) analytical bias
    bias_utility_score_analytical_stage1   = bias_utility_score_analytical_stage1,
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
