#' @title Internal utility helper functions
#' @description Core computational helpers used throughout the package.
#' @keywords internal
NULL

#' Compute Joint Multinomial Probabilities
#'
#' Computes the four joint probabilities for binary efficacy (X) and binary
#' safety (Y) outcomes given marginal rates and their Pearson correlation.
#'
#' @param p Numeric scalar in (0,1). Response rate \eqn{\Pr(X=1)}.
#' @param q Numeric scalar in (0,1). No-adverse-event rate \eqn{\Pr(Y=1)}.
#' @param phi Numeric scalar in [-1,1]. Pearson correlation between X and Y.
#'
#' @return A numeric vector of length 4: \eqn{(\pi_1, \pi_2, \pi_3, \pi_4)}
#'   corresponding to (Response+NoAE, Response+AE, NoResponse+NoAE,
#'   NoResponse+AE).
#'
#' @details
#' The joint probability is \eqn{\pi_1 = pq + \phi\sqrt{p(1-p)q(1-q)}},
#' with the remaining probabilities derived by subtraction. The achievable
#' range of \code{phi} is constrained by the Frechet-Hoeffding bounds:
#' \deqn{\phi_{\min} = \frac{\max(0, p+q-1) - pq}{\sqrt{p(1-p)q(1-q)}}}
#' \deqn{\phi_{\max} = \frac{\min(p,q) - pq}{\sqrt{p(1-p)q(1-q)}}}
#'
#' @examples
#' calc_pi(p = 0.4, q = 0.8, phi = 0)
#' calc_pi(p = 0.4, q = 0.8, phi = -0.3)
#'
#' @export
calc_pi <- function(p, q, phi = 0) {
  if (!is.numeric(p) || p <= 0 || p >= 1) stop("'p' must be in (0, 1).")
  if (!is.numeric(q) || q <= 0 || q >= 1) stop("'q' must be in (0, 1).")
  if (!is.numeric(phi) || phi < -1 || phi > 1) stop("'phi' must be in [-1, 1].")

  sqrt_term <- sqrt(p * (1 - p) * q * (1 - q))
  pi1 <- p * q + phi * sqrt_term
  pi2 <- p - pi1
  pi3 <- q - pi1
  pi4 <- 1 - p - q + pi1

  pi_vec <- c(pi1, pi2, pi3, pi4)
  if (any(pi_vec < -1e-9)) {
    stop(sprintf(
      "Invalid joint probabilities (phi incompatible with p=%.3f, q=%.3f). ",
      p, q
    ))
  }
  pmax(pi_vec, 0)
}


#' Feasible Correlation Bounds
#'
#' Returns the Frechet-Hoeffding lower and upper bounds on the Pearson
#' correlation \code{phi} between binary efficacy and safety endpoints.
#'
#' @inheritParams calc_pi
#'
#' @return A named numeric vector with elements \code{phi_min} and \code{phi_max}.
#'
#' @examples
#' phi_bounds(p = 0.4, q = 0.8)
#'
#' @export
phi_bounds <- function(p, q) {
  if (!is.numeric(p) || p <= 0 || p >= 1) stop("'p' must be in (0, 1).")
  if (!is.numeric(q) || q <= 0 || q >= 1) stop("'q' must be in (0, 1).")
  denom <- sqrt(p * (1 - p) * q * (1 - q))
  phi_min <- (max(0, p + q - 1) - p * q) / denom
  phi_max <- (min(p, q) - p * q) / denom
  c(phi_min = phi_min, phi_max = phi_max)
}


#' Compute Utility Scores from Trade-off Ratio
#'
#' Derives utility scores \eqn{(u_1, u_2, u_3, u_4)} from the clinical
#' trade-off ratio \eqn{r = \delta / d} using the margin-based approach
#' with utility independence.
#'
#' @param r Positive numeric. Trade-off ratio \eqn{r = \delta/d}, where
#'   \eqn{\delta} is the clinically meaningful efficacy margin and \eqn{d}
#'   is the safety margin.
#'
#' @return A numeric vector of length 4:
#'   \eqn{u_1 = 1, u_2 = 1/(1+r), u_3 = r/(1+r), u_4 = 0}.
#'
#' @details
#' When \eqn{r < 1} (efficacy margin < safety margin), efficacy is weighted
#' more and \eqn{u_2 > u_3}. When \eqn{r = 1}, both endpoints contribute
#' equally. When \eqn{r > 1}, safety dominates; the user should verify
#' that \eqn{u_2 \ge u_3} remains appropriate.
#'
#' @examples
#' # Equal weighting
#' calc_utility(r = 1)
#'
#' # Efficacy-dominant (r < 1)
#' calc_utility(r = 0.5)
#'
#' @export
calc_utility <- function(r) {
  if (!is.numeric(r) || length(r) != 1 || r <= 0) stop("'r' must be a single positive number.")
  c(u1 = 1, u2 = 1 / (1 + r), u3 = r / (1 + r), u4 = 0)
}


#' Compute Mean and Variance of Utility Score
#'
#' Given a probability vector over the four outcomes and a utility score
#' vector, computes the mean and variance of a single patient's utility score.
#'
#' @param pi Numeric vector of length 4. Joint probabilities
#'   \eqn{(\pi_1, \pi_2, \pi_3, \pi_4)} from \code{\link{calc_pi}}.
#' @param u Numeric vector of length 4. Utility scores
#'   \eqn{(u_1, u_2, u_3, u_4)} satisfying \eqn{u_1 \ge u_2 \ge u_3 \ge u_4}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{mu}{Expected utility \eqn{\mu = \sum_k u_k \pi_k}.}
#'     \item{sigma2}{Variance \eqn{\sigma^2 = \sum_k u_k^2 \pi_k - \mu^2}
#'       (floored at 1e-10).}
#'   }
#'
#' @examples
#' pi_vec <- calc_pi(0.4, 0.8, phi = 0)
#' u_vec  <- calc_utility(r = 1)
#' calc_utility_moments(pi_vec, u_vec)
#'
#' @export
calc_utility_moments <- function(pi, u) {
  if (!is.numeric(pi) || length(pi) != 4) stop("'pi' must be a numeric vector of length 4.")
  if (!is.numeric(u)  || length(u)  != 4) stop("'u' must be a numeric vector of length 4.")
  mu     <- sum(u * pi)
  sigma2 <- sum(u^2 * pi) - mu^2
  list(mu = mu, sigma2 = max(sigma2, 1e-10))
}
