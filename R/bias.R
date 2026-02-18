#' Analytical Selection-Induced Bias in Response Rate
#'
#' Computes the analytical selection-induced bias in response rate estimation
#' arising from a two-arm dose optimization study where the dose with the
#' higher observed mean utility is selected. Implements Equation 21 of the
#' manuscript.
#'
#' @param p Numeric in (0,1). True response rate under the null hypothesis
#'   (equal for both doses).
#' @param q Numeric in (0,1). True no-adverse-event rate under the null.
#' @param phi Numeric in [-1,1]. Pearson correlation between efficacy and
#'   safety endpoints.
#' @param u Numeric vector of length 4. Utility scores.
#' @param n1 Positive integer. Stage 1 (dose optimization) sample size per arm.
#' @param n2 Non-negative integer. Stage 2 (confirmatory) sample size.
#' @param lambda_u Non-negative numeric. Decision threshold. Default 0.
#'
#' @return A named list:
#'   \describe{
#'     \item{bias_stage1}{Bias in the selected dose's Stage 1 response rate
#'       estimate (Eq. 21).}
#'     \item{bias_combined}{Bias in the pooled Stage 1 + Stage 2 estimate
#'       (Eq. 25), diluted by \eqn{n_1/(n_1+n_2)}.}
#'     \item{bias_max_stage1}{Upper bound ignoring utility: \eqn{\sqrt{p(1-p)/n_1\pi}}
#'       (Eq. 22).}
#'     \item{bias_max_combined}{Upper bound for the combined estimator.}
#'     \item{cov_XU}{Covariance between efficacy indicator X and utility U.}
#'     \item{sigma_U}{Standard deviation of the utility score.}
#'     \item{mu_U}{Mean of the utility score.}
#'   }
#'
#' @details
#' Under the null hypothesis \eqn{p_L = p_H = p}, \eqn{q_L = q_H = q},
#' both arms have identical distributions. The bias arises from always
#' selecting the arm with the larger observed utility:
#' \deqn{\text{Bias} = \frac{\text{Cov}(X, U)}{\sigma_U \sqrt{n_1}} \cdot
#'   \frac{1}{\sqrt{\pi}} \exp\!\left(-\frac{\lambda_u^2 n_1}{4\sigma_U^2}\right)}
#'
#' When \eqn{\lambda_u = 0} and \eqn{U = X} (ROSE special case), this
#' reduces to \eqn{\sqrt{p(1-p) / (n_1 \pi)}}, the classical result of
#' Bauer et al. (2010).
#'
#' @references
#' Bauer P, Koenig F, Brannath W, Posch M (2010). Selection and bias—two
#' hostile brothers. \emph{Statistics in Medicine}, 29(1):1--13.
#'
#' @examples
#' u_vec <- calc_utility(r = 1)
#' calc_bias(p = 0.4, q = 0.8, phi = 0, u = u_vec, n1 = 60, n2 = 140)
#'
#' # Maximum bound (efficacy-only selection, lambda_u = 0)
#' calc_bias(p = 0.4, q = 0.8, phi = 0, u = c(1,1,0,0), n1 = 60, n2 = 140)
#'
#' @export
calc_bias <- function(p, q, phi = 0, u, n1, n2, lambda_u = 0) {
  if (!is.numeric(p) || p <= 0 || p >= 1) stop("'p' must be in (0,1).")
  if (!is.numeric(q) || q <= 0 || q >= 1) stop("'q' must be in (0,1).")
  if (!is.numeric(u) || length(u) != 4)   stop("'u' must be a numeric vector of length 4.")
  if (!is.numeric(n1) || n1 <= 0)          stop("'n1' must be a positive number.")
  if (!is.numeric(n2) || n2 < 0)           stop("'n2' must be non-negative.")
  if (!is.numeric(lambda_u) || lambda_u < 0) stop("'lambda_u' must be non-negative.")

  pi_vec <- calc_pi(p, q, phi)
  mom    <- calc_utility_moments(pi_vec, u)
  mu_U   <- mom$mu; sigma2_U <- mom$sigma2; sigma_U <- sqrt(sigma2_U)

  x_vec  <- c(1, 1, 0, 0)  # efficacy indicators per outcome
  E_XU   <- sum(x_vec * u * pi_vec)
  cov_XU <- E_XU - p * mu_U

  exp_term <- exp(-(lambda_u^2 * n1) / (4 * sigma2_U))

  bias_stage1   <- (cov_XU / (sigma_U * sqrt(n1))) * (1 / sqrt(base::pi)) * exp_term
  bias_combined <- (n1 / (n1 + n2)) * bias_stage1

  bias_max_stage1   <- sqrt(p * (1 - p)) / sqrt(n1 * base::pi)
  bias_max_combined <- (n1 / (n1 + n2)) * bias_max_stage1

  list(
    bias_stage1        = bias_stage1,
    bias_combined      = bias_combined,
    bias_max_stage1    = bias_max_stage1,
    bias_max_combined  = bias_max_combined,
    cov_XU             = cov_XU,
    sigma_U            = sigma_U,
    mu_U               = mu_U
  )
}


#' TTE Bias Bound Using Surrogacy Correlation
#'
#' Estimates (or bounds) the selection-induced bias in a time-to-event (TTE)
#' confirmatory endpoint based on its correlation with the binary efficacy
#' endpoint used in the dose optimization stage.
#'
#' @param bias_binary Numeric. Binary bias from \code{\link{calc_bias}}
#'   (\code{bias_combined} or \code{bias_stage1}).
#' @param sigma_X Numeric. Standard deviation of the binary efficacy
#'   indicator, \eqn{\sqrt{p(1-p)}}.
#' @param sigma_T Numeric. Standard deviation of the TTE endpoint or
#'   landmark survival indicator.
#' @param rho_TX Numeric in [-1,1]. Surrogacy correlation between TTE and
#'   binary efficacy. Use 1 for the conservative upper bound.
#'
#' @return A named list:
#'   \describe{
#'     \item{bias_TTE_adjusted}{Surrogate-adjusted TTE bias (Eq. 36).}
#'     \item{bias_TTE_upper}{Conservative upper bound (\eqn{\rho_{TX}=1}, Eq. 34).}
#'   }
#'
#' @details
#' From Theorem 1 of the manuscript:
#' \deqn{\text{Bias}(\bar{T}) \approx \rho_{TX} \cdot \frac{\sigma_T}{\sigma_X}
#'   \cdot \text{Bias}(\hat{p})}
#' with the upper bound obtained by setting \eqn{\rho_{TX} = 1}.
#'
#' @examples
#' b <- calc_bias(p = 0.4, q = 0.8, phi = 0, u = c(1,1,0,0), n1 = 60, n2 = 140)
#' # Landmark survival at 24 weeks, ~0.095 SD assuming S0 = 0.9
#' calc_bias_TTE(bias_binary = b$bias_combined, sigma_X = sqrt(0.4*0.6),
#'               sigma_T = sqrt(0.9*0.1), rho_TX = 0.5)
#'
#' @export
calc_bias_TTE <- function(bias_binary, sigma_X, sigma_T, rho_TX = 1) {
  if (!is.numeric(rho_TX) || abs(rho_TX) > 1) stop("'rho_TX' must be in [-1,1].")
  list(
    bias_TTE_adjusted = rho_TX * (sigma_T / sigma_X) * bias_binary,
    bias_TTE_upper    = (sigma_T / sigma_X) * bias_binary
  )
}
