# ============================================================================
# Sample Size Determination: Normal Approximation (Modes 1-3, ROSE)
# ============================================================================

#' Sample Size for Utility-Based Dose Optimization (Normal Approximation)
#'
#' @description
#' Computes the minimum per-arm sample size \eqn{n} to achieve prespecified
#' Probabilities of Correct Selection (PCS) under two clinical scenarios,
#' using the normal approximation. Supports three specification modes.
#'
#' @param pL Numeric. Base response rate for the reference dose. Required for
#'   Modes 1 and 2.
#' @param qL Numeric. Base no-AE rate for the reference dose. Required for
#'   Modes 1 and 2.
#' @param delta Numeric. Clinically meaningful efficacy margin (> 0). Required
#'   for Mode 1.
#' @param d Numeric. Clinically meaningful safety margin (>= 0). Required for
#'   Mode 1.
#' @param pH_L Numeric. High-dose response rate under \eqn{\mathcal{S}_L}
#'   (Mode 2 only).
#' @param qH_L Numeric. High-dose no-AE rate under \eqn{\mathcal{S}_L}
#'   (Mode 2 only).
#' @param pL_H Numeric. Low-dose response rate under \eqn{\mathcal{S}_H}
#'   (Mode 2 only).
#' @param qL_H Numeric. Low-dose no-AE rate under \eqn{\mathcal{S}_H}
#'   (Mode 2 only).
#' @param pL_L,qL_L,pH_L_fd,qH_L_fd,pL_H_fd,qL_H_fd,pH_H,qH_H Numeric.
#'   Full rate specification for Mode 3 (8 parameters, all dose-scenario
#'   combinations).
#' @param phi Numeric in [-1, 1]. Pearson correlation between efficacy and
#'   safety binary indicators. Default 0.
#' @param alpha_L Numeric in (0, 1). Target PCS under Scenario L (probability
#'   of correctly selecting Dose L when it is superior). Default 0.8.
#' @param alpha_H Numeric in (0, 1). Target PCS under Scenario H (probability
#'   of correctly selecting Dose H when it is superior). Default 0.8.
#' @param u Numeric vector of length 4. User-supplied utility scores
#'   \eqn{(u_1, u_2, u_3, u_4)} with \eqn{u_1 \ge u_2 \ge u_3 \ge u_4}.
#'   If \code{NULL} (default), derived automatically from \eqn{r = \delta/d}
#'   (Modes 1-2) or required (Mode 3).
#'
#' @return An object of class \code{"dose_design"} with elements:
#' \describe{
#'   \item{\code{n}}{Required per-arm sample size.}
#'   \item{\code{lambda_u}}{Optimal decision threshold \eqn{\lambda_u^*}.}
#'   \item{\code{PCS_L}}{Achieved PCS under Scenario L.}
#'   \item{\code{PCS_H}}{Achieved PCS under Scenario H.}
#'   \item{\code{utility}}{Utility score vector \eqn{(u_1, u_2, u_3, u_4)}.}
#'   \item{\code{delta}, \code{d}}{Efficacy and safety margins used.}
#'   \item{\code{scenario_L}, \code{scenario_H}}{Lists with rates and variance
#'     under each scenario.}
#'   \item{\code{inputs}}{List of all input parameters.}
#' }
#'
#' @details
#' \strong{Sample size formula} (Equation 12 in Gu et al., 2025):
#' \deqn{n = \left[\frac{z_{\alpha_L}\sqrt{v(S_L)} - z_{1-\alpha_H}\sqrt{v(S_H)}}
#'   {\Delta\mu(S_H) - \Delta\mu(S_L)}\right]^2}
#'
#' \strong{Mode 1 (Margin-based, recommended):} Specify \code{pL}, \code{qL},
#' \code{delta}, \code{d}. Utility scores are derived via
#' \code{\link{utility_scores_from_margins}} and scenarios are constructed
#' automatically. This ensures symmetric utility differences across scenarios.
#'
#' \strong{Mode 2 (Partial direct):} Specify \code{pL}, \code{qL} plus all
#' four scenario-specific rates. Utility can be specified or derived.
#'
#' \strong{Mode 3 (Full direct):} Specify all 8 rates across both scenarios
#' plus user-supplied utility scores.
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). Equations (10)-(14).
#'
#' Wang S, Yuan Y, Liu S (2025). Biometrics, 81(4), ujaf124.
#'
#' @seealso \code{\link{ss_utility_exact}}, \code{\link{ss_rose_approx}},
#'   \code{\link{ss_direct_approx}}
#'
#' @examples
#' # Mode 1: margin-based (recommended)
#' res <- ss_utility_approx(pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
#'                          alpha_L = 0.8, alpha_H = 0.8)
#' print(res)
#'
#' # Mode 1 with negative efficacy-safety correlation
#' res2 <- ss_utility_approx(pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
#'                           phi = -0.3, alpha_L = 0.8, alpha_H = 0.8)
#'
#' # Mode 2: partial direct specification
#' res3 <- ss_utility_approx(pL = 0.3, qL = 0.8,
#'                           pH_L = 0.3, qH_L = 0.65,
#'                           pL_H = 0.2, qL_H = 0.8,
#'                           alpha_L = 0.8, alpha_H = 0.8)
#'
#' @export
ss_utility_approx <- function(
    pL = NULL, qL = NULL,
    delta = NULL, d = NULL,
    pH_L = NULL, qH_L = NULL,
    pL_H = NULL, qL_H = NULL,
    pL_L = NULL, qL_L = NULL,
    pH_L_fd = NULL, qH_L_fd = NULL,
    pL_H_fd = NULL, qL_H_fd = NULL,
    pH_H = NULL, qH_H = NULL,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8,
    u = NULL
) {
  # --- Detect mode ---
  mode_out <- .detect_mode(pL, qL, delta, d, pH_L, qH_L, pL_H, qL_H,
                            pL_L, qL_L, pH_L_fd, qH_L_fd,
                            pL_H_fd, qL_H_fd, pH_H, qH_H, u)

  pL_L <- mode_out$pL_L; qL_L <- mode_out$qL_L
  pH_L_use <- mode_out$pH_L; qH_L_use <- mode_out$qH_L
  pL_H_use <- mode_out$pL_H; qL_H_use <- mode_out$qL_H
  pH_H_use <- mode_out$pH_H; qH_H_use <- mode_out$qH_H
  u_vec <- mode_out$u
  delta_out <- mode_out$delta; d_out <- mode_out$d

  # Validate phi and alpha
  .check_scalar_range(phi, "phi", -1, 1)
  .check_scalar_range(alpha_L, "alpha_L", 0, 1, open = TRUE)
  .check_scalar_range(alpha_H, "alpha_H", 0, 1, open = TRUE)

  # Compute moments under each scenario
  pi_LL <- .joint_probs_internal(pL_L, qL_L, phi)
  pi_HL <- .joint_probs_internal(pH_L_use, qH_L_use, phi)
  pi_LH <- .joint_probs_internal(pL_H_use, qL_H_use, phi)
  pi_HH <- .joint_probs_internal(pH_H_use, qH_H_use, phi)

  .validate_pi(pi_LL, "Scenario L, Dose L")
  .validate_pi(pi_HL, "Scenario L, Dose H")
  .validate_pi(pi_LH, "Scenario H, Dose L")
  .validate_pi(pi_HH, "Scenario H, Dose H")

  mom_LL <- .utility_moments_internal(pi_LL, u_vec)
  mom_HL <- .utility_moments_internal(pi_HL, u_vec)
  mom_LH <- .utility_moments_internal(pi_LH, u_vec)
  mom_HH <- .utility_moments_internal(pi_HH, u_vec)

  delta_mu_L <- mom_HL$mu - mom_LL$mu
  delta_mu_H <- mom_HH$mu - mom_LH$mu
  v_L <- mom_LL$sigma2 + mom_HL$sigma2
  v_H <- mom_LH$sigma2 + mom_HH$sigma2

  if (delta_mu_L >= 0)
    warning("Under Scenario L, Dose H is not inferior to Dose L (delta_mu_L >= 0).")
  if (delta_mu_H <= 0)
    warning("Under Scenario H, Dose H is not superior to Dose L (delta_mu_H <= 0).")
  if ((delta_mu_H - delta_mu_L) <= 0)
    stop("Utility difference (delta_mu_H - delta_mu_L) is not positive. Check inputs.")

  z_aL  <- stats::qnorm(alpha_L)
  z_1aH <- stats::qnorm(1 - alpha_H)
  n_cont <- ((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / (delta_mu_H - delta_mu_L))^2
  n <- ceiling(n_cont)
  if (n < 2) stop("Calculated sample size is less than 2.")

  lambda_u <- delta_mu_H + z_1aH * sqrt(v_H / n)
  pcs_L <- stats::pnorm((lambda_u - delta_mu_L) / sqrt(v_L / n))
  pcs_H <- 1 - stats::pnorm((lambda_u - delta_mu_H) / sqrt(v_H / n))

  structure(
    list(
      n = n, lambda_u = lambda_u, PCS_L = pcs_L, PCS_H = pcs_H,
      utility = u_vec, delta = delta_out, d = d_out,
      scenario_L = list(
        delta_mu = delta_mu_L, variance = v_L,
        pL = pL_L, qL = qL_L, pH = pH_L_use, qH = qH_L_use,
        pi_L = pi_LL, pi_H = pi_HL
      ),
      scenario_H = list(
        delta_mu = delta_mu_H, variance = v_H,
        pL = pL_H_use, qL = qL_H_use, pH = pH_H_use, qH = qH_H_use,
        pi_L = pi_LH, pi_H = pi_HH
      ),
      method = "approximate",
      inputs = list(
        mode = mode_out$mode, phi = phi,
        alpha_L = alpha_L, alpha_H = alpha_H
      )
    ),
    class = c("dose_design", "list")
  )
}


#' Sample Size for ROSE Design (Normal Approximation)
#'
#' @description
#' Computes the minimum per-arm sample size for the Randomized Optimal
#' Selection (ROSE) design using the normal approximation. ROSE is a
#' special case of the utility framework with \eqn{u = (1, 1, 0, 0)}
#' (selection based on efficacy only).
#'
#' @param pL Numeric in (0, 1). Response rate for both doses under
#'   Scenario L (and Dose H under Scenario H).
#' @param delta Numeric > 0. Efficacy margin (response rate advantage of
#'   Dose H under Scenario H): \eqn{p_H - p_L = \delta}.
#' @param alpha_L Numeric in (0, 1). Target PCS under Scenario L. Default 0.8.
#' @param alpha_H Numeric in (0, 1). Target PCS under Scenario H. Default 0.8.
#'
#' @return An object of class \code{"dose_design"} with the same structure
#'   as \code{\link{ss_utility_approx}}.
#'
#' @details
#' The ROSE design corresponds to utility scores \eqn{u_1 = u_2 = 1},
#' \eqn{u_3 = u_4 = 0}, so \eqn{\bar{U}_j = \hat{p}_j}. Scenario L has
#' \eqn{\Delta\mu(S_L) = 0}, making the optimal threshold
#' \eqn{\lambda_u^* > 0} necessary for adequate PCS under both scenarios
#' (Wang et al., 2025, Table 1).
#'
#' @references
#' Wang S, Yuan Y, Liu S (2025). Biometrics, 81(4), ujaf124.
#'
#' Gu X, Xu C, Xu L, Yuan Y (2025). Section 3.5.
#'
#' @seealso \code{\link{ss_rose_exact}}, \code{\link{ss_utility_approx}}
#'
#' @examples
#' # Reproduce ROSE design from Wang et al. (2025) Table 1
#' res <- ss_rose_approx(pL = 0.4, delta = 0.15, alpha_L = 0.8, alpha_H = 0.8)
#' print(res)  # n = 58, lambda_u ~ 0.077
#'
#' @export
ss_rose_approx <- function(pL, delta, alpha_L = 0.8, alpha_H = 0.8) {
  .check_scalar_range(pL, "pL", 0, 1, open = TRUE)
  .check_scalar_range(delta, "delta", 0, Inf, left_open = TRUE)
  if (pL - delta < 0) stop("pL - delta must be >= 0.")
  .check_scalar_range(alpha_L, "alpha_L", 0, 1, open = TRUE)
  .check_scalar_range(alpha_H, "alpha_H", 0, 1, open = TRUE)

  pL_H <- pL - delta; pH_H <- pL
  v_L <- 2 * pL * (1 - pL)
  v_H <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)

  z_aL  <- stats::qnorm(alpha_L)
  z_1aH <- stats::qnorm(1 - alpha_H)
  n_cont <- ((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / delta)^2
  n <- ceiling(n_cont)

  lambda_u <- delta + z_1aH * sqrt(v_H / n)
  pcs_L <- stats::pnorm(lambda_u / sqrt(v_L / n))
  pcs_H <- 1 - stats::pnorm((lambda_u - delta) / sqrt(v_H / n))

  structure(
    list(
      n = n, lambda_u = lambda_u, PCS_L = pcs_L, PCS_H = pcs_H,
      utility = c(u1 = 1, u2 = 1, u3 = 0, u4 = 0), delta = delta, d = NA,
      scenario_L = list(delta_mu = 0, variance = v_L, pL = pL, pH = pL),
      scenario_H = list(delta_mu = delta, variance = v_H, pL = pL_H, pH = pH_H),
      method = "rose_approximate",
      inputs = list(pL = pL, delta = delta, alpha_L = alpha_L, alpha_H = alpha_H)
    ),
    class = c("dose_design", "list")
  )
}


# ============================================================================
# Internal helpers
# ============================================================================

# Detect specification mode and standardize all rate parameters
.detect_mode <- function(pL, qL, delta, d, pH_L, qH_L, pL_H, qL_H,
                          pL_L, qL_L, pH_L_fd, qH_L_fd,
                          pL_H_fd, qL_H_fd, pH_H, qH_H, u) {

  mode1 <- !is.null(delta) || !is.null(d)
  mode2 <- !is.null(pL) && !is.null(qL) &&
    (!is.null(pH_L) || !is.null(qH_L) || !is.null(pL_H) || !is.null(qL_H))
  mode3 <- !is.null(pL_L) && !is.null(qL_L) && !is.null(pH_H) && !is.null(qH_H)

  if (sum(mode1, mode2, mode3) > 1)
    stop("Cannot mix parameters from different modes. Use one of Mode 1, 2, or 3.")
  if (sum(mode1, mode2, mode3) == 0)
    stop("Must specify parameters for Mode 1 (delta, d), Mode 2, or Mode 3.")

  delta_out <- NA; d_out <- NA; mode_name <- ""
  u_vec <- u

  if (mode1) {
    mode_name <- "margin-based"
    if (is.null(pL) || is.null(qL) || is.null(delta) || is.null(d))
      stop("Mode 1 requires: pL, qL, delta, d.")
    .check_scalar_range(pL, "pL", 0, 1, open = TRUE)
    .check_scalar_range(qL, "qL", 0, 1, open = TRUE)
    .check_scalar_range(delta, "delta", 0, Inf, left_open = TRUE)
    .check_scalar_range(d, "d", 0, Inf, left_open = TRUE)
    if (pL + delta > 1) stop("pL + delta must not exceed 1.")
    if (qL - d < 0) stop("qL - d must be >= 0.")
    pL_L_out <- pL; qL_L_out <- qL
    pH_L_out <- pL;  qH_L_out <- qL - d
    pL_H_out <- pL - delta; qL_H_out <- qL
    pH_H_out <- pL;  qH_H_out <- qL
    delta_out <- delta; d_out <- d
    if (is.null(u_vec)) {
      u_vec <- .utility_from_r(delta / d)
    }
  } else if (mode2) {
    mode_name <- "partial-direct"
    if (is.null(pH_L) || is.null(qH_L) || is.null(pL_H) || is.null(qL_H))
      stop("Mode 2 requires: pL, qL, pH_L, qH_L, pL_H, qL_H.")
    pL_L_out <- pL; qL_L_out <- qL
    pH_L_out <- pH_L; qH_L_out <- qH_L
    pL_H_out <- pL_H; qL_H_out <- qL_H
    pH_H_out <- pL; qH_H_out <- qL
    delta_out <- pL - pL_H; d_out <- qL - qH_L
    if (is.null(u_vec)) {
      if (abs(d_out) < 1e-10) stop("Safety margin (qL - qH_L) is zero; supply u manually.")
      u_vec <- .utility_from_r(delta_out / d_out)
    }
  } else {
    mode_name <- "full-direct"
    if (is.null(pL_L) || is.null(qL_L) || is.null(pH_L_fd) || is.null(qH_L_fd) ||
        is.null(pL_H_fd) || is.null(qL_H_fd) || is.null(pH_H) || is.null(qH_H))
      stop("Mode 3 requires all 8 rate parameters.")
    if (is.null(u_vec)) stop("Mode 3 requires user-supplied utility scores u.")
    pL_L_out <- pL_L; qL_L_out <- qL_L
    pH_L_out <- pH_L_fd; qH_L_out <- qH_L_fd
    pL_H_out <- pL_H_fd; qL_H_out <- qL_H_fd
    pH_H_out <- pH_H; qH_H_out <- qH_H
  }

  u_vec <- .parse_utility(u_vec)

  list(
    mode = mode_name,
    pL_L = pL_L_out, qL_L = qL_L_out,
    pH_L = pH_L_out, qH_L = qH_L_out,
    pL_H = pL_H_out, qL_H = qL_H_out,
    pH_H = pH_H_out, qH_H = qH_H_out,
    u = u_vec, delta = delta_out, d = d_out
  )
}

# Check scalar in range
.check_scalar_range <- function(x, name, lo, hi, open = FALSE,
                                 left_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x))
    stop(sprintf("'%s' must be a single non-missing numeric value.", name))
  lo_ok <- if (open || left_open) x > lo else x >= lo
  hi_ok <- if (open) x < hi else x <= hi
  if (!lo_ok || !hi_ok) {
    bracket_l <- if (open || left_open) "(" else "["
    bracket_r <- if (open) ")" else "]"
    stop(sprintf("'%s' must be in %s%.4g, %.4g%s.",
                 name, bracket_l, lo, hi, bracket_r))
  }
}
