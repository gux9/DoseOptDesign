# ============================================================================
# Utility Score Core Functions
# ============================================================================

#' Compute Joint Outcome Probabilities
#'
#' @description
#' Computes the joint probability vector \eqn{\pi = (\pi_1, \pi_2, \pi_3, \pi_4)}
#' for four mutually exclusive efficacy-safety outcome combinations from marginal
#' rates and their Pearson correlation.
#'
#' @param p Numeric scalar in (0, 1). Response rate \eqn{p = \Pr(X = 1)}.
#' @param q Numeric scalar in (0, 1). No-AE rate \eqn{q = \Pr(Y = 1)}.
#' @param phi Numeric scalar in [-1, 1]. Pearson correlation between the binary
#'   efficacy indicator \eqn{X} and the binary safety indicator \eqn{Y}.
#'   Default is 0 (independence).
#'
#' @return A named numeric vector of length 4:
#' \describe{
#'   \item{\code{pi1}}{Response AND no AE: \eqn{\pi_1 = pq + \phi\sqrt{p(1-p)q(1-q)}}}
#'   \item{\code{pi2}}{Response AND AE: \eqn{\pi_2 = p - \pi_1}}
#'   \item{\code{pi3}}{No response AND no AE: \eqn{\pi_3 = q - \pi_1}}
#'   \item{\code{pi4}}{No response AND AE: \eqn{\pi_4 = 1 - p - q + \pi_1}}
#' }
#'
#' @details
#' The joint probability vector determines the distribution of the categorical
#' random vector \eqn{\mathbf{W} = (W_1, W_2, W_3, W_4) \sim \text{Categorical}(\pi)}.
#' The parameter \eqn{\phi} must lie within the Frechet-Hoeffding bounds; see
#' \code{\link{phi_bounds}} for the achievable range.
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). A Utility Score Framework for Dose
#' Optimization Studies. Equation (4).
#'
#' @seealso \code{\link{phi_bounds}}, \code{\link{utility_moments}}
#'
#' @examples
#' # Independent efficacy and safety
#' joint_probs(p = 0.4, q = 0.8, phi = 0)
#'
#' # Negative correlation (response tends to accompany AE)
#' joint_probs(p = 0.4, q = 0.8, phi = -0.3)
#'
#' @export
joint_probs <- function(p, q, phi = 0) {
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1)
    stop("'p' must be a single numeric value in (0, 1).")
  if (!is.numeric(q) || length(q) != 1 || q <= 0 || q >= 1)
    stop("'q' must be a single numeric value in (0, 1).")
  if (!is.numeric(phi) || length(phi) != 1 || phi < -1 || phi > 1)
    stop("'phi' must be a single numeric value in [-1, 1].")

  sqrt_term <- sqrt(p * (1 - p) * q * (1 - q))
  pi1 <- p * q + phi * sqrt_term
  pi2 <- p - pi1
  pi3 <- q - pi1
  pi4 <- 1 - p - q + pi1
  pi_vec <- c(pi1 = pi1, pi2 = pi2, pi3 = pi3, pi4 = pi4)

  if (any(pi_vec < -1e-10))
    stop(sprintf(
      "phi = %.4f produces invalid joint probabilities (some pi_k < 0). ",
      phi,
      "Use phi_bounds(p, q) to check the feasible range."
    ))

  pmax(pi_vec, 0)
}


#' Frechet-Hoeffding Bounds for the Correlation Parameter
#'
#' @description
#' Computes the minimum and maximum achievable Pearson correlation \eqn{\phi}
#' between the binary efficacy indicator \eqn{X} and the binary safety indicator
#' \eqn{Y}, given marginal probabilities \eqn{p} and \eqn{q}.
#'
#' @param p Numeric scalar in (0, 1). Response rate.
#' @param q Numeric scalar in (0, 1). No-AE rate.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{\code{phi_min}}{Lower Frechet-Hoeffding bound.}
#'   \item{\code{phi_max}}{Upper Frechet-Hoeffding bound.}
#' }
#'
#' @details
#' The joint probability \eqn{\pi_1 = \Pr(X=1, Y=1)} must satisfy
#' \eqn{\max(0, p+q-1) \le \pi_1 \le \min(p, q)}, which translates to
#' bounds on \eqn{\phi} via the formula in \code{\link{joint_probs}}.
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). Equation (5).
#'
#' @seealso \code{\link{joint_probs}}
#'
#' @examples
#' phi_bounds(p = 0.4, q = 0.8)
#'
#' @export
phi_bounds <- function(p, q) {
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1)
    stop("'p' must be a single numeric value in (0, 1).")
  if (!is.numeric(q) || length(q) != 1 || q <= 0 || q >= 1)
    stop("'q' must be a single numeric value in (0, 1).")

  denom <- sqrt(p * (1 - p) * q * (1 - q))
  phi_min <- (max(0, p + q - 1) - p * q) / denom
  phi_max <- (min(p, q) - p * q) / denom
  c(phi_min = phi_min, phi_max = phi_max)
}


#' Derive Utility Scores from Clinical Margins
#'
#' @description
#' Derives the four utility scores \eqn{(u_1, u_2, u_3, u_4)} from a pair of
#' clinically meaningful margins \eqn{\delta} (efficacy) and \eqn{d} (safety)
#' via the margin-based approach with utility independence.
#'
#' @param delta Positive numeric. Clinically meaningful efficacy margin (response
#'   rate improvement that justifies selecting the higher dose).
#' @param d Positive numeric. Clinically meaningful safety margin (no-AE rate
#'   improvement). Must be greater than 0.
#'
#' @return A named numeric vector of length 4:
#' \describe{
#'   \item{\code{u1}}{Best outcome (Response, No AE): always 1.}
#'   \item{\code{u2}}{Response with AE: \eqn{1/(1+r)}.}
#'   \item{\code{u3}}{No response, No AE: \eqn{r/(1+r)}.}
#'   \item{\code{u4}}{Worst outcome (No response, AE): always 0.}
#' }
#' where \eqn{r = \delta/d}.
#'
#' @details
#' The trade-off ratio \eqn{r = \delta/d} encodes the relative clinical importance
#' of efficacy versus safety:
#' \itemize{
#'   \item \eqn{r < 1} (\eqn{\delta < d}): AEs tolerable; \eqn{u_2 > u_3}.
#'     Typical in serious oncology settings.
#'   \item \eqn{r = 1}: Equal weighting; \eqn{u_2 = u_3 = 0.5}.
#'   \item \eqn{r > 1}: Safety-dominant; \eqn{u_3 > u_2}. Uncommon in oncology.
#' }
#'
#' The utility independence constraint \eqn{u_1 - u_2 - u_3 + u_4 = 0} ensures
#' a constant marginal rate of substitution (MRS = \eqn{r}) across all \eqn{(p, q)}
#' and \eqn{\phi} values (Gu et al., 2025, Section 2.4).
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). Equation (9).
#'
#' @seealso \code{\link{utility_moments}}, \code{\link{ss_utility_approx}}
#'
#' @examples
#' # Efficacy margin 0.10, safety margin 0.15 => r = 2/3
#' utility_scores_from_margins(delta = 0.10, d = 0.15)
#'
#' # Equal margins
#' utility_scores_from_margins(delta = 0.15, d = 0.15)
#'
#' @export
utility_scores_from_margins <- function(delta, d) {
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0)
    stop("'delta' must be a positive numeric scalar.")
  if (!is.numeric(d) || length(d) != 1 || d <= 0)
    stop("'d' must be a positive numeric scalar.")

  r <- delta / d
  u <- c(u1 = 1, u2 = 1 / (1 + r), u3 = r / (1 + r), u4 = 0)
  attr(u, "r") <- r
  u
}


#' Compute Mean and Variance of the Utility Score
#'
#' @description
#' Computes the mean \eqn{\mu = E[U]} and variance \eqn{\sigma^2 = \text{Var}(U)}
#' of the per-patient utility score \eqn{U = \mathbf{u}^T \mathbf{W}}.
#'
#' @param pi Numeric vector of length 4. Joint outcome probabilities
#'   \eqn{(\pi_1, \pi_2, \pi_3, \pi_4)} summing to 1; see \code{\link{joint_probs}}.
#' @param u Numeric vector of length 4. Utility scores satisfying
#'   \eqn{u_1 \ge u_2 \ge u_3 \ge u_4 \ge 0}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{mu}}{Mean utility \eqn{\mu = \sum_k u_k \pi_k}.}
#'   \item{\code{sigma2}}{Variance \eqn{\sigma^2 = \sum_k u_k^2 \pi_k - \mu^2}.
#'     Floored at \eqn{10^{-10}} for numerical stability.}
#'   \item{\code{sigma}}{Standard deviation \eqn{\sigma = \sqrt{\sigma^2}}.}
#' }
#'
#' @details
#' For a sample of \eqn{n} patients, the sample mean utility \eqn{\bar{U}} is
#' approximately \eqn{N(\mu, \sigma^2/n)} by the Central Limit Theorem. This
#' asymptotic normality underlies the approximate sample size formulas in
#' \code{\link{ss_utility_approx}}.
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). Equation (6).
#'
#' @seealso \code{\link{joint_probs}}, \code{\link{utility_scores_from_margins}}
#'
#' @examples
#' pi_vec <- joint_probs(p = 0.4, q = 0.8, phi = 0)
#' u_vec  <- utility_scores_from_margins(delta = 0.10, d = 0.15)
#' utility_moments(pi_vec, u_vec)
#'
#' @export
utility_moments <- function(pi, u) {
  if (!is.numeric(pi) || length(pi) != 4)
    stop("'pi' must be a numeric vector of length 4.")
  if (any(pi < -1e-10) || abs(sum(pi) - 1) > 1e-6)
    stop("'pi' must be non-negative and sum to 1.")
  if (!is.numeric(u) || length(u) != 4)
    stop("'u' must be a numeric vector of length 4.")

  mu <- sum(u * pi)
  sigma2 <- max(sum(u^2 * pi) - mu^2, 1e-10)
  list(mu = mu, sigma2 = sigma2, sigma = sqrt(sigma2))
}


#' Estimate Parameters from Observed 2x2 Count Data
#'
#' @description
#' Estimates the response rate \eqn{p}, no-AE rate \eqn{q}, and efficacy-safety
#' correlation \eqn{\phi} from observed 2x2 count data in a prior trial.
#'
#' @param n11 Non-negative integer. Number of patients with response AND no AE.
#' @param n10 Non-negative integer. Number of patients with response AND AE.
#' @param n01 Non-negative integer. Number of patients with no response AND no AE.
#' @param n00 Non-negative integer. Number of patients with no response AND AE.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{\code{p}}{Estimated response rate \eqn{\hat{p}}.}
#'   \item{\code{q}}{Estimated no-AE rate \eqn{\hat{q}}.}
#'   \item{\code{phi}}{Estimated Pearson correlation \eqn{\hat{\phi}},
#'     truncated to the Frechet-Hoeffding bounds if necessary.}
#'   \item{\code{n}}{Total sample size.}
#' }
#'
#' @details
#' Uses the standard Pearson chi-square-type correlation estimate for 2x2
#' contingency tables (Eq. 3 in Gu et al., 2025). The estimate is truncated
#' to the achievable range \eqn{[\phi_{\min}, \phi_{\max}]} via
#' \code{\link{phi_bounds}}.
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). Equation (3).
#'
#' @seealso \code{\link{phi_bounds}}, \code{\link{joint_probs}}
#'
#' @examples
#' # From a prior trial with 100 patients
#' estimate_params(n11 = 30, n10 = 10, n01 = 45, n00 = 15)
#'
#' @export
estimate_params <- function(n11, n10, n01, n00) {
  for (nm in c("n11", "n10", "n01", "n00")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1 || val < 0 || val != floor(val))
      stop(sprintf("'%s' must be a non-negative integer.", nm))
  }
  n <- n11 + n10 + n01 + n00
  if (n == 0) stop("Total count must be > 0.")

  p_hat <- (n11 + n10) / n
  q_hat <- (n11 + n01) / n

  num <- n * n11 - (n11 + n10) * (n11 + n01)
  den <- sqrt((n11 + n10) * (n11 + n01) * (n10 + n00) * (n01 + n00))
  phi_hat <- if (abs(den) < sqrt(.Machine$double.eps)) 0 else num / den

  # Truncate to achievable range
  if (p_hat > 0 && p_hat < 1 && q_hat > 0 && q_hat < 1) {
    bounds <- phi_bounds(p_hat, q_hat)
    phi_hat <- max(bounds["phi_min"], min(bounds["phi_max"], phi_hat))
  }

  c(p = p_hat, q = q_hat, phi = phi_hat, n = n)
}


# ============================================================================
# Internal helper: compute joint probabilities (vectorized, no input checks)
# ============================================================================
.joint_probs_internal <- function(p, q, phi = 0) {
  sqrt_term <- sqrt(pmax(p * (1 - p) * q * (1 - q), 0))
  pi1 <- p * q + phi * sqrt_term
  pi2 <- p - pi1
  pi3 <- q - pi1
  pi4 <- 1 - p - q + pi1
  c(pi1, pi2, pi3, pi4)
}

# Internal helper: compute utility moments (no checks)
.utility_moments_internal <- function(pi, u) {
  mu <- sum(u * pi)
  sigma2 <- max(sum(u^2 * pi) - mu^2, 1e-10)
  list(mu = mu, sigma2 = sigma2)
}

# Internal helper: validate pi vector
.validate_pi <- function(pi, label = "") {
  if (any(pi < -1e-10) || abs(sum(pi) - 1) > 1e-6)
    stop(sprintf("Invalid joint probabilities%s: pi = (%s), sum = %.6f",
                 if (nchar(label) > 0) paste0(" for ", label) else "",
                 paste(round(pi, 6), collapse = ", "), sum(pi)))
}

# Internal helper: compute utility from r
.utility_from_r <- function(r) c(1, 1 / (1 + r), r / (1 + r), 0)

# Internal helper: parse and validate utility vector
.parse_utility <- function(u, label = "u") {
  if (!is.numeric(u) || length(u) != 4)
    stop(sprintf("'%s' must be a numeric vector of length 4.", label))
  if (any(diff(u) > 1e-10))
    message("Note: utility scores are not in non-increasing order (u1 >= u2 >= u3 >= u4).")
  u
}
