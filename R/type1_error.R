#' Type I Error Inflation for Binary Confirmatory Tests
#'
#' Computes the inflated Type I error rate for one-sided Z-test and exact
#' binomial test when Stage 1 dose optimization data are pooled with Stage 2
#' confirmatory data.
#'
#' @param p0 Numeric in (0,1). Null response rate.
#' @param bias_combined Numeric. Combined-stage bias from
#'   \code{\link{calc_bias}} (\code{bias_combined}).
#' @param n_total Positive integer. Total combined sample size \eqn{n_1 + n_2}.
#' @param alpha Numeric in (0,1). Nominal one-sided significance level.
#'   Default 0.025.
#'
#' @return A named list:
#'   \describe{
#'     \item{type1_z}{Inflated Type I error for the Z-test (Eq. 26).}
#'     \item{type1_bin}{Inflated Type I error for the binomial test (Eq. 28).}
#'     \item{k_c}{Binomial critical value (smallest \eqn{k_c} such that
#'       \eqn{\Pr(X > k_c | n, p_0) \le \alpha}).}
#'     \item{se_null}{Standard error under the null.}
#'   }
#'
#' @details
#' The Z-test inflation is:
#' \deqn{\text{Type I Error}_Z = 1 - \Phi\left(z_{1-\alpha} -
#'   \frac{\Delta p_{\text{comb}}}{\text{SE}_0}\right)}
#'
#' The binomial test inflation accounts for the discrete nature of the
#' rejection region (Eq. 28):
#' \deqn{\text{Type I Error}_{\text{Bin}} = \Pr(X > k_c \mid n, p_0 + \Delta p_{\text{comb}})}
#'
#' The binomial test generally shows lower inflation (simulation mean 1.34×)
#' than the Z-test (1.86×) and its plugin estimator is uniformly conservative.
#'
#' @examples
#' u_vec <- calc_utility(r = 1)
#' b <- calc_bias(p = 0.4, q = 0.8, phi = 0, u = u_vec, n1 = 60, n2 = 140)
#' calc_type1_error(p0 = 0.4, bias_combined = b$bias_combined,
#'                  n_total = 200, alpha = 0.025)
#'
#' # With maximum (conservative) bound
#' calc_type1_error(p0 = 0.4, bias_combined = b$bias_max_combined,
#'                  n_total = 200, alpha = 0.025)
#'
#' @export
calc_type1_error <- function(p0, bias_combined, n_total, alpha = 0.025) {
  if (!is.numeric(p0) || p0 <= 0 || p0 >= 1) stop("'p0' must be in (0,1).")
  if (!is.numeric(n_total) || n_total <= 0) stop("'n_total' must be positive.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0,1).")

  se_null <- sqrt(p0 * (1 - p0) / n_total)
  z_crit  <- qnorm(1 - alpha)

  type1_z <- 1 - pnorm(z_crit - bias_combined / se_null)

  # Correct binomial critical value
  probs <- stats::pbinom(0:n_total, size = n_total, prob = p0, lower.tail = FALSE)
  k_c   <- min(which(probs <= alpha))  # k_c: smallest rejection threshold

  type1_bin <- stats::pbinom(k_c - 1, size = n_total,
                              prob = p0 + bias_combined, lower.tail = FALSE)

  list(type1_z = type1_z, type1_bin = type1_bin, k_c = k_c, se_null = se_null)
}


#' Type I Error Inflation for Landmark Survival Z-Test
#'
#' Computes the inflated Type I error for a one-sided Z-test of the
#' \eqn{\tau}-month landmark survival rate when Stage 1 and Stage 2 TTE data
#' are pooled (Equation 31 of the manuscript).
#'
#' @param S0 Numeric in (0,1). Null landmark survival probability.
#' @param bias_landmark_combined Numeric. Bias in combined landmark survival
#'   rate estimate.
#' @param n_total Positive integer. Total sample size \eqn{n_1 + n_2}.
#' @param alpha Numeric. Nominal one-sided significance level. Default 0.025.
#'
#' @return A named list with \code{type1_landmark} and \code{bias_z_shift}.
#'
#' @examples
#' b <- calc_bias(p = 0.4, q = 0.8, phi = 0, u = c(1,1,0,0), n1 = 60, n2 = 140)
#' # Suppose landmark survival SD = 0.3, rho_TX = 0.5
#' b_tte <- calc_bias_TTE(b$bias_combined, sqrt(0.4*0.6), 0.3, rho_TX = 0.5)
#' calc_type1_error_landmark(S0 = 0.9, bias_landmark_combined = b_tte$bias_TTE_adjusted,
#'                           n_total = 200)
#'
#' @export
calc_type1_error_landmark <- function(S0, bias_landmark_combined, n_total, alpha = 0.025) {
  if (!is.numeric(S0) || S0 <= 0 || S0 >= 1) stop("'S0' must be in (0,1).")
  se_null   <- sqrt(S0 * (1 - S0) / n_total)
  z_crit    <- qnorm(1 - alpha)
  bias_z    <- bias_landmark_combined / se_null
  type1     <- 1 - pnorm(z_crit - bias_z)
  list(type1_landmark = type1, bias_z_shift = bias_z)
}


#' Type I Error Inflation for One-Arm Exponential (Log-Scale) Test
#'
#' Computes the inflated Type I error for a one-sided one-arm exponential
#' test on the log-hazard scale (Equations 32-33 of the manuscript).
#'
#' @param bias_mean_T_combined Numeric. Selection-induced bias in mean survival
#'   time for the combined dataset.
#' @param lambda0 Numeric > 0. Null hazard rate.
#' @param expected_events Numeric. Expected number of events under the null,
#'   computed from study design parameters.
#' @param alpha Numeric. Nominal one-sided significance level. Default 0.025.
#'
#' @return A named list with \code{type1_exp} and the log-hazard bias.
#'
#' @details
#' The log-hazard bias is \eqn{-\lambda_0 \cdot \text{Bias}(\bar{T})_{\text{comb}}},
#' and the test statistic bias is this times \eqn{\sqrt{D}}.
#'
#' @export
calc_type1_error_exp <- function(bias_mean_T_combined, lambda0, expected_events,
                                  alpha = 0.025) {
  if (!is.numeric(lambda0) || lambda0 <= 0) stop("'lambda0' must be positive.")
  bias_log_lambda <- -lambda0 * bias_mean_T_combined
  bias_Z          <- bias_log_lambda * sqrt(expected_events)
  z_crit          <- qnorm(1 - alpha)
  type1           <- pnorm(-z_crit - bias_Z)
  list(type1_exp = type1, bias_log_lambda = bias_log_lambda, bias_Z_exp = bias_Z)
}


#' Type I Error Inflation for Two-Arm Cox Regression Test
#'
#' Computes the inflated Type I error for a two-arm Cox test via the
#' delta method chain from mean survival bias to log hazard ratio bias
#' (Proposition 1 of the manuscript).
#'
#' @param bias_mean_T_combined Numeric. Bias in mean survival for combined data.
#' @param lambda0 Numeric > 0. Null hazard rate.
#' @param total_events Numeric. Total expected events across both arms.
#' @param alpha Numeric. Nominal one-sided significance level. Default 0.025.
#'
#' @return A named list with \code{type1_cox}, log-HR bias, and Cox Z bias.
#'
#' @examples
#' # Illustrative example
#' calc_type1_error_cox(bias_mean_T_combined = 0.05, lambda0 = 0.1,
#'                      total_events = 100, alpha = 0.025)
#'
#' @export
calc_type1_error_cox <- function(bias_mean_T_combined, lambda0, total_events,
                                  alpha = 0.025) {
  if (!is.numeric(lambda0) || lambda0 <= 0) stop("'lambda0' must be positive.")
  # Step 2: hazard bias
  bias_lambda   <- -lambda0^2 * bias_mean_T_combined
  # Step 3: log-HR bias
  bias_log_HR   <- bias_lambda / lambda0  # = -lambda0 * bias_mean_T
  # Step 4: Cox Z bias (Fisher info under equal allocation = D_total/4)
  V_cox         <- pmax(total_events / 4, 1)
  bias_Z_cox    <- bias_log_HR * sqrt(V_cox)
  z_crit        <- qnorm(1 - alpha)
  type1         <- pnorm(-z_crit - bias_Z_cox)
  list(type1_cox = type1, bias_log_HR = bias_log_HR, bias_Z_cox = bias_Z_cox)
}


#' Expected Number of Events for One-Arm Exponential Study
#'
#' Computes the expected number of events under exponential survival with
#' uniform accrual and administrative censoring.
#'
#' @param n_total Positive integer. Total number of patients.
#' @param lambda0 Numeric > 0. Hazard rate.
#' @param entry_time_max Numeric > 0. Maximum enrollment time (weeks).
#' @param admin_censor_time Numeric > 0. Administrative censoring calendar
#'   time (weeks).
#'
#' @return Expected number of events (numeric scalar).
#'
#' @examples
#' calc_expected_events(n_total = 200, lambda0 = 0.1,
#'                      entry_time_max = 52, admin_censor_time = 76)
#'
#' @export
calc_expected_events <- function(n_total, lambda0, entry_time_max, admin_censor_time) {
  lam <- lambda0
  n_total * (1 - exp(-lam * admin_censor_time) *
               (exp(lam * entry_time_max) - 1) / (lam * entry_time_max))
}
