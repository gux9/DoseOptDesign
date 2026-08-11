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
#' @export
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
