#' Sample Size for Utility-Based Dose Selection (Normal Approximation)
#'
#' Computes the per-arm sample size for a one-stage two-arm dose optimization
#' study using the normal approximation. Three specification modes are supported.
#'
#' @param pL Numeric. Base response rate (Dose L).
#' @param qL Numeric. Base no-adverse-event rate (Dose L).
#' @param delta Numeric. Clinically meaningful efficacy margin \eqn{\delta > 0}.
#' @param d Numeric. Clinically meaningful safety margin \eqn{d \ge 0}.
#' @param pL_L,qL_L,pH_L,qH_L Numeric. Response and no-AE rates for Dose L
#'   and Dose H under Scenario \eqn{S_L} (Mode 2/3).
#' @param pL_H,qL_H,pH_H,qH_H Numeric. Response and no-AE rates for Dose L
#'   and Dose H under Scenario \eqn{S_H} (Mode 2/3).
#' @param phi Numeric in [-1,1]. Pearson correlation between efficacy and
#'   safety endpoints. Default 0.
#' @param alpha_L Numeric in (0,1). Target PCS under Scenario \eqn{S_L}.
#'   Default 0.8.
#' @param alpha_H Numeric in (0,1). Target PCS under Scenario \eqn{S_H}.
#'   Default 0.8.
#' @param u Numeric vector of length 4. Utility scores \eqn{(u_1, u_2, u_3, u_4)}.
#'   If \code{NULL} and Mode 1 is used, derived automatically from \eqn{r = \delta/d}
#'   via \code{\link{calc_utility}}.
#'
#' @details
#' \strong{Mode 1 (Margin-based):} Supply \code{pL}, \code{qL}, \code{delta},
#' \code{d}. Scenarios are constructed as:
#' \itemize{
#'   \item \eqn{S_L}: Dose L = \eqn{(p, q)}, Dose H = \eqn{(p, q-d)}.
#'   \item \eqn{S_H}: Dose L = \eqn{(p-\delta, q)}, Dose H = \eqn{(p, q)}.
#' }
#'
#' \strong{Mode 2 (Partial-direct):} Supply \code{pL}, \code{qL},
#' \code{pH_L}, \code{qH_L}, \code{pL_H}, \code{qL_H}. Dose H under \eqn{S_H}
#' and Dose L under \eqn{S_L} default to \code{(pL, qL)}.
#'
#' \strong{Mode 3 (Full-direct):} Supply all eight rate parameters explicitly,
#' and you must also supply \code{u}.
#'
#' The jointly optimal sample size formula (Equation 13 of the manuscript) is:
#' \deqn{n = \left[\frac{z_{\alpha_L}\sqrt{v_{S_L}} - z_{1-\alpha_H}\sqrt{v_{S_H}}}
#'              {\Delta\mu_{S_H} - \Delta\mu_{S_L}}\right]^2}
#'
#' @return An object of class \code{c("dose_design", "list")} with elements:
#'   \describe{
#'     \item{n}{Required sample size per arm.}
#'     \item{lambda_u}{Optimal decision threshold.}
#'     \item{PCS_L, PCS_H}{Achieved PCS under each scenario.}
#'     \item{scenario_L, scenario_H}{Named lists with scenario parameters.}
#'     \item{utility}{Utility score vector used.}
#'     \item{r, delta, d}{Trade-off ratio and margins (if applicable).}
#'     \item{inputs}{All input parameters.}
#'   }
#'
#' @seealso \code{\link{calc_sample_size_utility_exact}},
#'   \code{\link{calc_sample_size_rose_approx}},
#'   \code{\link{calc_utility}}, \code{\link{calc_pi}}
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2026). A Utility Score Framework for Dose
#' Optimization Studies with Binary Efficacy-Safety Endpoints.
#'
#' Wang S, Yuan Y, Liu S (2025). Randomized Optimal Selection Design for
#' Dose Optimization. \emph{Biometrics}, 81(4):ujaf124.
#'
#' @examples
#' # Mode 1: margin-based
#' res <- calc_sample_size_utility_approx(
#'   pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
#'   phi = 0, alpha_L = 0.8, alpha_H = 0.8
#' )
#' print(res)
#'
#' # ROSE special case (utility based on efficacy only)
#' res_rose <- calc_sample_size_utility_approx(
#'   pL = 0.4, qL = 0.8, delta = 0.15, d = 0.15,
#'   phi = 0, alpha_L = 0.8, alpha_H = 0.8,
#'   u = c(1, 1, 0, 0)
#' )
#' print(res_rose)
#'
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
  if (!is.numeric(phi) || phi < -1 || phi > 1) stop("'phi' must be in [-1, 1].")
  if (!is.numeric(alpha_L) || alpha_L <= 0 || alpha_L >= 1) stop("'alpha_L' must be in (0, 1).")
  if (!is.numeric(alpha_H) || alpha_H <= 0 || alpha_H >= 1) stop("'alpha_H' must be in (0, 1).")

  parsed <- .parse_dose_params(pL, qL, delta, d, pL_L, qL_L, pH_L, qH_L,
                                pL_H, qL_H, pH_H, qH_H, u)
  pL_L <- parsed$pL_L; qL_L <- parsed$qL_L; pH_L <- parsed$pH_L; qH_L <- parsed$qH_L
  pL_H <- parsed$pL_H; qL_H <- parsed$qL_H; pH_H <- parsed$pH_H; qH_H <- parsed$qH_H
  u    <- parsed$u;    r    <- parsed$r

  pi_L_L <- calc_pi(pL_L, qL_L, phi); pi_H_L <- calc_pi(pH_L, qH_L, phi)
  pi_L_H <- calc_pi(pL_H, qL_H, phi); pi_H_H <- calc_pi(pH_H, qH_H, phi)

  mom_L_L <- calc_utility_moments(pi_L_L, u); mom_H_L <- calc_utility_moments(pi_H_L, u)
  mom_L_H <- calc_utility_moments(pi_L_H, u); mom_H_H <- calc_utility_moments(pi_H_H, u)

  delta_mu_L <- mom_H_L$mu - mom_L_L$mu
  v_L        <- mom_H_L$sigma2 + mom_L_L$sigma2
  delta_mu_H <- mom_H_H$mu - mom_L_H$mu
  v_H        <- mom_H_H$sigma2 + mom_L_H$sigma2

  if (delta_mu_L >= 0) warning("Scenario S_L: Dose H is not inferior to Dose L in utility.")
  if (delta_mu_H <= 0) warning("Scenario S_H: Dose H is not superior to Dose L in utility.")
  denom <- delta_mu_H - delta_mu_L
  if (denom <= 0) stop("Utility difference between scenarios is not positive.")

  z_aL  <- qnorm(alpha_L)
  z_1aH <- qnorm(1 - alpha_H)
  n <- ceiling(((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / denom)^2)
  if (n < 2) stop("Calculated sample size < 2. Check input parameters.")

  lambda_u <- delta_mu_H + z_1aH * sqrt(v_H / n)
  pcs_L <- pnorm((lambda_u - delta_mu_L) / sqrt(v_L / n))
  pcs_H <- 1 - pnorm((lambda_u - delta_mu_H) / sqrt(v_H / n))

  structure(
    list(
      n = n, lambda_u = lambda_u, PCS_L = pcs_L, PCS_H = pcs_H,
      scenario_L = list(delta_mu = delta_mu_L, v = v_L,
                        pL = pL_L, qL = qL_L, pH = pH_L, qH = qH_L,
                        pi_L = pi_L_L, pi_H = pi_H_L),
      scenario_H = list(delta_mu = delta_mu_H, v = v_H,
                        pL = pL_H, qL = qL_H, pH = pH_H, qH = qH_H,
                        pi_L = pi_L_H, pi_H = pi_H_H),
      utility = u, r = r, delta = parsed$delta, d = parsed$d,
      method = "approximate",
      inputs = list(mode = parsed$mode_name, phi = phi,
                    alpha_L = alpha_L, alpha_H = alpha_H)
    ),
    class = c("dose_design", "list")
  )
}


#' Sample Size for Utility-Based Dose Selection (Exact Multinomial)
#'
#' Computes the exact minimum per-arm sample size via grid search over
#' the exact Multinomial probability mass function, as described in
#' Equations 16-18 of the manuscript.
#'
#' @inheritParams calc_sample_size_utility_approx
#' @param max_n Integer. Upper bound for sample size search. Default 500.
#' @param buffer Integer. Number of sample sizes below the approximate
#'   solution to begin searching. Default 10.
#' @param den Integer. Denominator for integer-scaling of utility scores
#'   (e.g., \code{den = 10} maps \eqn{u_2 = 0.8} to integer 8). Higher
#'   values give finer resolution but increase computation. Default 10.
#' @param diff_method Character. Method for computing the PMF of the utility
#'   difference: \code{"fft"} (faster, default for large n) or
#'   \code{"nested"} (slower but exact). Default \code{"fft"}.
#'
#' @return An object of class \code{c("dose_design_exact", "dose_design", "list")}
#'   with the same elements as \code{\link{calc_sample_size_utility_approx}},
#'   plus:
#'   \describe{
#'     \item{search_info}{List: approximate n, approximate lambda_u, number
#'       of sample sizes tested, and den.}
#'   }
#'
#' @examples
#' res <- calc_sample_size_utility_exact(
#'   pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
#'   phi = 0, alpha_L = 0.8, alpha_H = 0.8,
#'   max_n = 200
#' )
#' print(res)
#'
#' @seealso \code{\link{calc_sample_size_utility_approx}}
#' @export
calc_sample_size_utility_exact <- function(
    pL = NULL, qL = NULL,
    delta = NULL, d = NULL,
    pL_L = NULL, qL_L = NULL,
    pH_L = NULL, qH_L = NULL,
    pL_H = NULL, qL_H = NULL,
    pH_H = NULL, qH_H = NULL,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8,
    u = NULL,
    max_n = 500, buffer = 10, den = 10,
    diff_method = c("fft", "nested")
) {
  diff_method <- match.arg(diff_method)
  if (!is.numeric(phi) || phi < -1 || phi > 1) stop("'phi' must be in [-1, 1].")

  parsed <- .parse_dose_params(pL, qL, delta, d, pL_L, qL_L, pH_L, qH_L,
                                pL_H, qL_H, pH_H, qH_H, u)
  pL_L <- parsed$pL_L; qL_L <- parsed$qL_L; pH_L <- parsed$pH_L; qH_L <- parsed$qH_L
  pL_H <- parsed$pL_H; qL_H <- parsed$qL_H; pH_H <- parsed$pH_H; qH_H <- parsed$qH_H
  u    <- parsed$u

  pi_L_SL <- calc_pi(pL_L, qL_L, phi); pi_H_SL <- calc_pi(pH_L, qH_L, phi)
  pi_L_SH <- calc_pi(pL_H, qL_H, phi); pi_H_SH <- calc_pi(pH_H, qH_H, phi)

  mom_L_L <- calc_utility_moments(pi_L_SL, u); mom_H_L <- calc_utility_moments(pi_H_SL, u)
  mom_L_H <- calc_utility_moments(pi_L_SH, u); mom_H_H <- calc_utility_moments(pi_H_SH, u)
  delta_mu_L <- mom_H_L$mu - mom_L_L$mu; v_L <- mom_H_L$sigma2 + mom_L_L$sigma2
  delta_mu_H <- mom_H_H$mu - mom_L_H$mu; v_H <- mom_H_H$sigma2 + mom_L_H$sigma2

  denom   <- delta_mu_H - delta_mu_L
  if (denom <= 0) stop("Utility difference between scenarios is not positive.")
  z_aL    <- qnorm(alpha_L); z_1aH <- qnorm(1 - alpha_H)
  n_approx   <- ceiling(((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / denom)^2)
  lambda_approx <- delta_mu_H + z_1aH * sqrt(v_H / n_approx)
  min_search <- max(5, n_approx - buffer)

  # For ROSE special case, exact integer arithmetic is trivial
  if (isTRUE(all.equal(u, c(1, 1, 0, 0)))) den <- 1
  u_int <- as.integer(round(u * den))

  found <- FALSE; n_tested <- 0
  for (n in min_search:max_n) {
    n_tested <- n_tested + 1
    max_S <- n * max(u_int)
    pmf_L_SL <- .compute_pmf_S(pi_L_SL, u_int, n)
    pmf_H_SL <- .compute_pmf_S(pi_H_SL, u_int, n)
    pmf_L_SH <- .compute_pmf_S(pi_L_SH, u_int, n)
    pmf_H_SH <- .compute_pmf_S(pi_H_SH, u_int, n)

    pmf_diff_SL <- .compute_pmf_diff(pmf_L_SL, pmf_H_SL, max_S, diff_method)
    pmf_diff_SH <- .compute_pmf_diff(pmf_L_SH, pmf_H_SH, max_S, diff_method)
    cdf_SL <- cumsum(pmf_diff_SL); cdf_SH <- cumsum(pmf_diff_SH)

    offset  <- max_S + 1
    q_SL_int <- min(which(cdf_SL >= alpha_L)) - offset
    q_SH_int <- max(which(cdf_SH <= 1 - alpha_H)) - offset

    if (q_SL_int <= q_SH_int) {
      found     <- TRUE
      n_final   <- n
      lambda_int <- round((q_SL_int + q_SH_int) / 2)
      lambda_u  <- lambda_int / (den * n)
      pcs_L     <- cdf_SL[lambda_int + offset]
      pcs_H     <- 1 - cdf_SH[lambda_int + offset]
      break
    }
  }
  if (!found) stop(sprintf("No solution found in [%d, %d]. Increase max_n.", min_search, max_n))

  structure(
    list(
      n = n_final, lambda_u = lambda_u, PCS_L = pcs_L, PCS_H = pcs_H,
      scenario_L = list(delta_mu = delta_mu_L, v = v_L,
                        pL = pL_L, qL = qL_L, pH = pH_L, qH = qH_L),
      scenario_H = list(delta_mu = delta_mu_H, v = v_H,
                        pL = pL_H, qL = qL_H, pH = pH_H, qH = qH_H),
      utility = u, r = parsed$r, delta = parsed$delta, d = parsed$d,
      method = sprintf("exact_%s", diff_method),
      inputs = list(mode = parsed$mode_name, phi = phi,
                    alpha_L = alpha_L, alpha_H = alpha_H, den = den),
      search_info = list(n_approx = n_approx, lambda_u_approx = lambda_approx,
                         n_tested = n_tested, den = den)
    ),
    class = c("dose_design_exact", "dose_design", "list")
  )
}


# ---- Internal helpers for exact method ----

#' @keywords internal
.compute_pmf_S <- function(pi, u_int, n) {
  max_S <- n * max(u_int)
  pmf   <- rep(0, max_S + 1); pmf[1] <- 1
  for (i in seq_len(n)) {
    new_pmf <- rep(0, max_S + 1)
    for (j in seq_along(u_int)) {
      if (pi[j] == 0) next
      s <- u_int[j]
      if (s == 0) {
        new_pmf <- new_pmf + pi[j] * pmf
      } else {
        new_pmf[(s + 1):(max_S + 1)] <-
          new_pmf[(s + 1):(max_S + 1)] + pi[j] * pmf[1:(max_S + 1 - s)]
      }
    }
    pmf <- new_pmf
  }
  pmf
}

#' @keywords internal
.compute_pmf_diff <- function(pmf_L, pmf_H, max_S, method = "fft") {
  if (method == "nested") {
    pmf_d  <- rep(0, 2 * max_S + 1); offset <- max_S + 1
    idx_L  <- which(pmf_L > 0) - 1; idx_H  <- which(pmf_H > 0) - 1
    for (sL in idx_L) {
      pL_val <- pmf_L[sL + 1]
      for (sH in idx_H) {
        d <- sH - sL
        pmf_d[d + offset] <- pmf_d[d + offset] + pL_val * pmf_H[sH + 1]
      }
    }
    pmf_d
  } else {
    pmf_L_rev     <- rev(pmf_L)
    target_length <- 2 * max_S + 1
    fft_length    <- 2^ceiling(log2(target_length + max(length(pmf_L_rev), length(pmf_H)) - 1))
    pmf_L_padded  <- c(pmf_L_rev, rep(0, fft_length - length(pmf_L_rev)))
    pmf_H_padded  <- c(pmf_H,     rep(0, fft_length - length(pmf_H)))
    conv_result   <- Re(stats::fft(
      stats::fft(pmf_L_padded) * stats::fft(pmf_H_padded), inverse = TRUE
    )) / fft_length
    start_idx <- length(pmf_L_rev)
    pmf_d     <- conv_result[start_idx:(start_idx + target_length - 1)]
    pmf_d[abs(pmf_d) < 1e-15] <- 0
    pmax(pmf_d, 0)
  }
}


# ---- Internal parameter parser ----
#' @keywords internal
.parse_dose_params <- function(pL, qL, delta, d,
                                pL_L, qL_L, pH_L, qH_L,
                                pL_H, qL_H, pH_H, qH_H, u) {
  margin_based  <- !is.null(delta) || !is.null(d)
  partial_direct <- !is.null(pL) && !is.null(qL) &&
    (!is.null(pH_L) || !is.null(qH_L) || !is.null(pL_H) || !is.null(qL_H))
  full_direct   <- !is.null(pL_L) && !is.null(qL_L) && !is.null(pL_H) &&
    !is.null(qL_H) && !is.null(pH_L) && !is.null(qH_L) &&
    !is.null(pH_H) && !is.null(qH_H)

  n_modes <- sum(margin_based, partial_direct, full_direct)
  if (n_modes > 1) stop("Cannot mix parameters from different specification modes.")
  if (n_modes == 0) stop("Must specify parameters for Mode 1, 2, or 3.")

  r <- NA; mode_name <- ""
  delta_out <- NA; d_out <- NA

  if (margin_based) {
    mode_name <- "margin-based"
    if (is.null(pL) || is.null(qL) || is.null(delta) || is.null(d))
      stop("Mode 1 requires: pL, qL, delta, d.")
    if (pL + delta > 1) stop("pL + delta must not exceed 1.")
    if (qL - d < 0)     stop("qL - d must be non-negative.")
    pL_L <- pL; qL_L <- qL; pH_L <- pL;          qH_L <- qL - d
    pL_H <- pL - delta; qL_H <- qL; pH_H <- pL; qH_H <- qL
    delta_out <- delta; d_out <- d
    if (is.null(u)) {
      if (d == 0) stop("'d' cannot be zero when deriving utility automatically.")
      r <- delta / d; u <- calc_utility(r)
    }
  }
  if (partial_direct) {
    mode_name <- "partial-direct"
    if (is.null(pH_L) || is.null(qH_L) || is.null(pL_H) || is.null(qL_H))
      stop("Mode 2 requires: pL, qL, pH_L, qH_L, pL_H, qL_H.")
    pL_L <- pL; qL_L <- qL; pH_H <- pL; qH_H <- qL
    delta_out <- pL - pL_H; d_out <- qL - qH_L
    if (is.null(u)) {
      if (d_out == 0) stop("Implied d = qL - qH_L cannot be zero.")
      r <- delta_out / d_out; u <- calc_utility(r)
    }
  }
  if (full_direct) {
    mode_name <- "full-direct"
    if (is.null(u)) stop("Mode 3 (full-direct) requires utility scores 'u'.")
    delta_out <- mean(c(pH_L - pL_L, pH_H - pL_H))
    d_out     <- mean(c(qL_L - qH_L, qL_H - qH_H))
  }

  if (!is.numeric(u) || length(u) != 4) stop("'u' must be a numeric vector of length 4.")

  list(pL_L = pL_L, qL_L = qL_L, pH_L = pH_L, qH_L = qH_L,
       pL_H = pL_H, qL_H = qL_H, pH_H = pH_H, qH_H = qH_H,
       u = u, r = r, delta = delta_out, d = d_out, mode_name = mode_name)
}
