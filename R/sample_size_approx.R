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
#' @export
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
