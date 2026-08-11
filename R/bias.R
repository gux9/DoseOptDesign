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
#' @export
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
#' @export
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
