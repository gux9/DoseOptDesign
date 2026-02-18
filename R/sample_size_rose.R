#' Sample Size for ROSE Design (Normal Approximation)
#'
#' Computes the per-arm sample size for the one-stage Randomized Optimal
#' Selection (ROSE) design using efficacy-only selection. This is a special
#' case of \code{\link{calc_sample_size_utility_approx}} with
#' \eqn{u = (1, 1, 0, 0)}.
#'
#' @param pL Numeric in (0,1). Response rate for the reference dose.
#' @param delta Numeric > 0. Clinically meaningful efficacy margin.
#' @param alpha_L Numeric in (0,1). Target PCS under Scenario \eqn{S_L}.
#'   Default 0.8.
#' @param alpha_H Numeric in (0,1). Target PCS under Scenario \eqn{S_H}.
#'   Default 0.8.
#'
#' @return An object of class \code{c("rose_design", "list")} with elements:
#'   \describe{
#'     \item{n}{Per-arm sample size.}
#'     \item{lambda_u}{Optimal decision threshold.}
#'     \item{PCS_L, PCS_H}{Achieved PCS values.}
#'     \item{delta}{Efficacy margin.}
#'     \item{method}{"rose_approximate".}
#'   }
#'
#' @references
#' Wang S, Yuan Y, Liu S (2025). Randomized Optimal Selection Design for
#' Dose Optimization. \emph{Biometrics}, 81(4):ujaf124.
#'
#' @examples
#' res <- calc_sample_size_rose_approx(pL = 0.4, delta = 0.15,
#'                                     alpha_L = 0.8, alpha_H = 0.8)
#' print(res)
#'
#' @export
calc_sample_size_rose_approx <- function(pL, delta,
                                          alpha_L = 0.8, alpha_H = 0.8) {
  if (!is.numeric(pL) || pL <= 0 || pL >= 1) stop("'pL' must be in (0,1).")
  if (!is.numeric(delta) || delta <= 0)       stop("'delta' must be positive.")
  if (pL - delta < 0) stop("pL - delta must be non-negative.")

  pH_H <- pL; pL_H <- pL - delta
  v_L  <- 2 * pL * (1 - pL)
  v_H  <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)

  z_aL  <- qnorm(alpha_L); z_1aH <- qnorm(1 - alpha_H)
  n     <- ceiling(((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / delta)^2)
  lam   <- delta + z_1aH * sqrt(v_H / n)
  pcs_L <- pnorm(lam / sqrt(v_L / n))
  pcs_H <- 1 - pnorm((lam - delta) / sqrt(v_H / n))

  structure(
    list(n = n, lambda_u = lam, PCS_L = pcs_L, PCS_H = pcs_H,
         utility = c(1, 1, 0, 0), delta = delta,
         method = "rose_approximate",
         inputs = list(alpha_L = alpha_L, alpha_H = alpha_H)),
    class = c("rose_design", "list")
  )
}


#' Sample Size for ROSE Design (Exact Binomial)
#'
#' Computes the exact minimum per-arm sample size for the ROSE design via
#' grid search over the exact Binomial PMF.
#'
#' @inheritParams calc_sample_size_rose_approx
#' @param max_n Integer. Upper bound for sample size search. Default 500.
#' @param buffer Integer. Search buffer below approximate solution. Default 10.
#'
#' @return An object of class \code{c("rose_design", "list")}.
#'
#' @examples
#' res <- calc_sample_size_rose_exact(pL = 0.4, delta = 0.15,
#'                                    alpha_L = 0.8, alpha_H = 0.8)
#' print(res)
#'
#' @export
calc_sample_size_rose_exact <- function(pL, delta,
                                         alpha_L = 0.8, alpha_H = 0.8,
                                         max_n = 500, buffer = 10) {
  if (!is.numeric(pL) || pL <= 0 || pL >= 1) stop("'pL' must be in (0,1).")
  if (!is.numeric(delta) || delta <= 0)       stop("'delta' must be positive.")
  if (pL - delta < 0) stop("pL - delta must be non-negative.")

  pL_H <- pL - delta; pH_H <- pL
  v_L  <- 2 * pL * (1 - pL)
  v_H  <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)
  z_aL <- qnorm(alpha_L); z_1aH <- qnorm(1 - alpha_H)
  n_approx  <- ceiling(((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / delta)^2)
  min_search <- max(5, n_approx - buffer)

  found <- FALSE
  for (n in min_search:max_n) {
    pmf_L_SL  <- stats::dbinom(0:n, n, pL)
    pmf_H_SL  <- stats::dbinom(0:n, n, pL)
    pmf_L_SH  <- stats::dbinom(0:n, n, pL_H)
    pmf_H_SH  <- stats::dbinom(0:n, n, pH_H)

    pmf_diff_SL <- .binom_diff_pmf(pmf_L_SL, pmf_H_SL, n)
    pmf_diff_SH <- .binom_diff_pmf(pmf_L_SH, pmf_H_SH, n)
    cdf_SL      <- cumsum(pmf_diff_SL)
    cdf_SH      <- cumsum(pmf_diff_SH)
    offset      <- n + 1

    q_SL <- min(which(cdf_SL >= alpha_L))   - offset
    q_SH <- max(which(cdf_SH <= 1 - alpha_H)) - offset

    if (q_SL <= q_SH) {
      found          <- TRUE
      threshold_cnt  <- round((q_SL + q_SH) / 2)
      lambda_u_final <- threshold_cnt / n
      pcs_L_final    <- cdf_SL[threshold_cnt + offset]
      pcs_H_final    <- 1 - cdf_SH[threshold_cnt + offset]
      break
    }
  }
  if (!found) stop(sprintf("No exact solution in [%d, %d]. Increase max_n.", min_search, max_n))

  structure(
    list(n = n, lambda_u = lambda_u_final, PCS_L = pcs_L_final, PCS_H = pcs_H_final,
         utility = c(1, 1, 0, 0), delta = delta, n_approx = n_approx,
         method = "rose_exact",
         inputs = list(alpha_L = alpha_L, alpha_H = alpha_H)),
    class = c("rose_design", "list")
  )
}

#' @keywords internal
.binom_diff_pmf <- function(pmf_L, pmf_H, n) {
  pmf_d  <- rep(0, 2 * n + 1); offset <- n + 1
  for (xL in 0:n) {
    if (pmf_L[xL + 1] == 0) next
    for (xH in 0:n) {
      d <- xH - xL
      pmf_d[d + offset] <- pmf_d[d + offset] + pmf_L[xL + 1] * pmf_H[xH + 1]
    }
  }
  pmf_d
}
