# ============================================================================
# Sample Size Determination: Exact Multinomial Method (Modes 1-3, ROSE)
# ============================================================================

#' Sample Size for Utility-Based Dose Optimization (Exact Multinomial Method)
#'
#' @description
#' Computes the minimum per-arm sample size to guarantee prespecified
#' Probabilities of Correct Selection (PCS) using exact Multinomial PMF
#' calculations via dynamic programming. Supports Modes 1-3.
#'
#' @inheritParams ss_utility_approx
#' @param max_n Positive integer. Maximum per-arm sample size to search.
#'   Default 500.
#' @param buffer Non-negative integer. Number of sample sizes below the
#'   approximate solution to begin the exact search. Default 10.
#' @param den Positive integer. Denominator for scaling utility scores to
#'   integers: \code{u_int = round(u * den)}. Larger values increase
#'   accuracy at a computational cost. Default 10.
#' @param diff_method Character. Method for computing the PMF of the utility
#'   difference: \code{"nested"} (exact, slow for large n) or \code{"fft"}
#'   (fast Fourier transform, faster). Default \code{"nested"}.
#'
#' @return An object of class \code{"dose_design_exact"} (inherits from
#'   \code{"dose_design"}) with all elements from \code{\link{ss_utility_approx}}
#'   plus:
#' \describe{
#'   \item{\code{search_info}}{List with \code{n_approx} (approximate solution),
#'     \code{lambda_u_approx}, \code{n_tested}, and \code{den}.}
#' }
#'
#' @details
#' The exact method directly computes:
#' \deqn{\text{PCS}_L = \sum_{\mathbf{n}_L, \mathbf{n}_H}
#'   \mathbf{1}(\bar{U}_H - \bar{U}_L \le \lambda_u)
#'   \cdot \text{Mult}(\mathbf{n}_L|\boldsymbol{\pi}_L^{S_L})
#'   \cdot \text{Mult}(\mathbf{n}_H|\boldsymbol{\pi}_H^{S_L})}
#' via a dynamic programming PMF convolution, then searches for the minimum
#' \eqn{n} satisfying both PCS constraints. The normal approximation is used
#' to initialize the search.
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2025). Equations (17)-(19).
#'
#' @seealso \code{\link{ss_utility_approx}}, \code{\link{ss_rose_exact}}
#'
#' @examples
#' res <- ss_utility_exact(pL = 0.3, qL = 0.7, delta = 0.10, d = 0.15,
#'                         alpha_L = 0.8, alpha_H = 0.8)
#' print(res)
#'
#' @export
ss_utility_exact <- function(
    pL = NULL, qL = NULL,
    delta = NULL, d = NULL,
    pH_L = NULL, qH_L = NULL,
    pL_H = NULL, qL_H = NULL,
    pL_L = NULL, qL_L = NULL,
    pH_L_fd = NULL, qH_L_fd = NULL,
    pL_H_fd = NULL, qL_H_fd = NULL,
    pH_H = NULL, qH_H = NULL,
    phi = 0, alpha_L = 0.8, alpha_H = 0.8,
    u = NULL,
    max_n = 500, buffer = 10, den = 10,
    diff_method = c("nested", "fft")
) {
  diff_method <- match.arg(diff_method)

  # Use approx to initialize and validate
  approx_res <- ss_utility_approx(
    pL = pL, qL = qL, delta = delta, d = d,
    pH_L = pH_L, qH_L = qH_L, pL_H = pL_H, qL_H = qL_H,
    pL_L = pL_L, qL_L = qL_L, pH_L_fd = pH_L_fd, qH_L_fd = qH_L_fd,
    pL_H_fd = pL_H_fd, qL_H_fd = qL_H_fd, pH_H = pH_H, qH_H = qH_H,
    phi = phi, alpha_L = alpha_L, alpha_H = alpha_H, u = u
  )

  u_vec    <- approx_res$utility
  n_approx <- approx_res$n
  lambda_u_approx <- approx_res$lambda_u

  sL <- approx_res$scenario_L
  sH <- approx_res$scenario_H
  pi_LL <- sL$pi_L; pi_HL <- sL$pi_H
  pi_LH <- sH$pi_L; pi_HH <- sH$pi_H

  # Scale utilities to integers
  if (all(u_vec == c(1, 1, 0, 0))) den <- 1
  u_int <- as.integer(round(u_vec * den))
  if (all(u_int == 0)) stop("Integer-scaled utilities all zero. Increase 'den'.")

  min_search <- max(2, n_approx - buffer)
  found <- FALSE
  n_tested <- 0

  for (n in min_search:max_n) {
    n_tested <- n_tested + 1
    max_S <- n * max(u_int)

    pmf_LL <- .pmf_multinomial_dp(pi_LL, u_int, n)
    pmf_HL <- .pmf_multinomial_dp(pi_HL, u_int, n)
    pmf_LH <- .pmf_multinomial_dp(pi_LH, u_int, n)
    pmf_HH <- .pmf_multinomial_dp(pi_HH, u_int, n)

    pmf_diff_SL <- .pmf_diff(pmf_LL, pmf_HL, max_S, diff_method)
    pmf_diff_SH <- .pmf_diff(pmf_LH, pmf_HH, max_S, diff_method)

    cdf_SL <- cumsum(pmf_diff_SL)
    cdf_SH <- cumsum(pmf_diff_SH)
    offset  <- max_S + 1

    q_SL_int <- min(which(cdf_SL >= alpha_L)) - offset
    q_SH_int <- max(which(cdf_SH <= 1 - alpha_H)) - offset

    if (q_SL_int <= q_SH_int) {
      found <- TRUE
      lambda_int <- as.integer(round((q_SL_int + q_SH_int) / 2))
      lambda_u_final <- lambda_int / (den * n)
      pcs_L_final <- cdf_SL[lambda_int + offset]
      pcs_H_final <- 1 - cdf_SH[lambda_int + offset]
      break
    }
  }

  if (!found)
    stop(sprintf(
      "No exact solution found in sample sizes [%d, %d]. Increase 'max_n'.",
      min_search, max_n
    ))

  result <- approx_res
  result$n         <- n
  result$lambda_u  <- lambda_u_final
  result$PCS_L     <- pcs_L_final
  result$PCS_H     <- pcs_H_final
  result$method    <- sprintf("exact_%s", diff_method)
  result$search_info <- list(
    n_approx        = n_approx,
    lambda_u_approx = lambda_u_approx,
    n_tested        = n_tested,
    den             = den
  )
  class(result) <- c("dose_design_exact", "dose_design", "list")
  result
}


#' Sample Size for ROSE Design (Exact Binomial Method)
#'
#' @description
#' Computes the exact minimum per-arm sample size for the ROSE design using
#' exact binomial probability calculations.
#'
#' @param pL Numeric in (0, 1). Common response rate under Scenario L.
#' @param delta Numeric > 0. Efficacy margin.
#' @param alpha_L Numeric in (0, 1). Target PCS under Scenario L. Default 0.8.
#' @param alpha_H Numeric in (0, 1). Target PCS under Scenario H. Default 0.8.
#' @param max_n Positive integer. Maximum sample size to search. Default 500.
#' @param buffer Non-negative integer. Search buffer below approximate solution.
#'   Default 10.
#'
#' @return An object of class \code{"dose_design"}.
#'
#' @references
#' Wang S, Yuan Y, Liu S (2025). Biometrics, 81(4), ujaf124.
#'
#' @seealso \code{\link{ss_rose_approx}}, \code{\link{ss_utility_exact}}
#'
#' @examples
#' res <- ss_rose_exact(pL = 0.4, delta = 0.15, alpha_L = 0.8, alpha_H = 0.8)
#' print(res)  # Matches Wang et al. (2025) Table 1
#'
#' @export
ss_rose_exact <- function(pL, delta, alpha_L = 0.8, alpha_H = 0.8,
                           max_n = 500, buffer = 10) {
  .check_scalar_range(pL, "pL", 0, 1, open = TRUE)
  .check_scalar_range(delta, "delta", 0, Inf, left_open = TRUE)
  if (pL - delta < 0) stop("pL - delta must be >= 0.")
  .check_scalar_range(alpha_L, "alpha_L", 0, 1, open = TRUE)
  .check_scalar_range(alpha_H, "alpha_H", 0, 1, open = TRUE)

  pL_H <- pL - delta; pH_H <- pL
  v_L  <- 2 * pL * (1 - pL)
  v_H  <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)
  z_aL <- stats::qnorm(alpha_L); z_1aH <- stats::qnorm(1 - alpha_H)
  n_approx <- ceiling(((z_aL * sqrt(v_L) - z_1aH * sqrt(v_H)) / delta)^2)
  min_search <- max(2, n_approx - buffer)
  found <- FALSE

  for (n in min_search:max_n) {
    pmf_LL <- stats::dbinom(0:n, n, pL)
    pmf_HL <- stats::dbinom(0:n, n, pL)
    pmf_LH <- stats::dbinom(0:n, n, pL_H)
    pmf_HH <- stats::dbinom(0:n, n, pH_H)

    pmf_diff_SL <- .pmf_diff_binomial(pmf_LL, pmf_HL, n)
    pmf_diff_SH <- .pmf_diff_binomial(pmf_LH, pmf_HH, n)

    cdf_SL <- cumsum(pmf_diff_SL)
    cdf_SH <- cumsum(pmf_diff_SH)
    offset  <- n + 1

    q_SL <- min(which(cdf_SL >= alpha_L)) - offset
    q_SH <- max(which(cdf_SH <= 1 - alpha_H)) - offset

    if (q_SL <= q_SH) {
      found <- TRUE
      threshold_count  <- as.integer(round((q_SL + q_SH) / 2))
      lambda_u_final   <- threshold_count / n
      pcs_L_final      <- cdf_SL[threshold_count + offset]
      pcs_H_final      <- 1 - cdf_SH[threshold_count + offset]
      break
    }
  }
  if (!found)
    stop(sprintf("No ROSE exact solution in [%d, %d]. Increase 'max_n'.", min_search, max_n))

  structure(
    list(
      n = n, lambda_u = lambda_u_final, PCS_L = pcs_L_final, PCS_H = pcs_H_final,
      utility = c(u1 = 1, u2 = 1, u3 = 0, u4 = 0), delta = delta, d = NA,
      scenario_L = list(delta_mu = 0, variance = v_L, pL = pL, pH = pL),
      scenario_H = list(delta_mu = delta, variance = v_H, pL = pL_H, pH = pH_H),
      method = "rose_exact",
      search_info = list(n_approx = n_approx),
      inputs = list(pL = pL, delta = delta, alpha_L = alpha_L, alpha_H = alpha_H)
    ),
    class = c("dose_design", "list")
  )
}


# ============================================================================
# Internal: Dynamic programming for Multinomial PMF
# ============================================================================

# Compute PMF of n-patient sum of integer utilities via DP
.pmf_multinomial_dp <- function(pi, u_int, n) {
  max_S <- n * max(u_int)
  pmf <- numeric(max_S + 1)
  pmf[1] <- 1.0
  for (i in seq_len(n)) {
    new_pmf <- numeric(max_S + 1)
    for (j in 1:4) {
      if (pi[j] < 1e-15) next
      s <- u_int[j]
      if (s == 0L) {
        new_pmf <- new_pmf + pi[j] * pmf
      } else {
        idx_to   <- (s + 1):(max_S + 1)
        idx_from <- 1:(max_S + 1 - s)
        new_pmf[idx_to] <- new_pmf[idx_to] + pi[j] * pmf[idx_from]
      }
    }
    pmf <- new_pmf
  }
  pmf
}

# Compute PMF of difference (score_H - score_L) by nested convolution
.pmf_diff_nested <- function(pmf_L, pmf_H, max_S) {
  len  <- 2 * max_S + 1
  offset <- max_S + 1
  pmf_d  <- numeric(len)
  idx_L  <- which(pmf_L > 1e-15) - 1
  idx_H  <- which(pmf_H > 1e-15) - 1
  for (sL in idx_L) {
    pL_val <- pmf_L[sL + 1]
    for (sH in idx_H) {
      d_idx <- sH - sL + offset
      pmf_d[d_idx] <- pmf_d[d_idx] + pL_val * pmf_H[sH + 1]
    }
  }
  pmf_d
}

# Compute PMF of difference using FFT (faster for large n)
.pmf_diff_fft <- function(pmf_L, pmf_H, max_S) {
  target_len  <- 2 * max_S + 1
  pmf_L_rev   <- rev(pmf_L)
  fft_len     <- 2^ceiling(log2(target_len + length(pmf_L_rev) - 1))
  fft_L <- stats::fft(c(pmf_L_rev, rep(0, fft_len - length(pmf_L_rev))))
  fft_H <- stats::fft(c(pmf_H,     rep(0, fft_len - length(pmf_H))))
  conv  <- Re(stats::fft(fft_L * fft_H, inverse = TRUE)) / fft_len
  start <- length(pmf_L_rev)
  pmf_d <- conv[start:(start + target_len - 1)]
  pmf_d[pmf_d < 0] <- 0
  pmf_d
}

.pmf_diff <- function(pmf_L, pmf_H, max_S, method = "nested") {
  if (method == "fft") .pmf_diff_fft(pmf_L, pmf_H, max_S)
  else .pmf_diff_nested(pmf_L, pmf_H, max_S)
}

# Binomial difference PMF for ROSE
.pmf_diff_binomial <- function(pmf_L, pmf_H, n) {
  len <- 2 * n + 1; offset <- n + 1
  pmf_d <- numeric(len)
  for (xL in 0:n) {
    if (pmf_L[xL + 1] < 1e-15) next
    for (xH in 0:n) {
      d_idx <- xH - xL + offset
      pmf_d[d_idx] <- pmf_d[d_idx] + pmf_L[xL + 1] * pmf_H[xH + 1]
    }
  }
  pmf_d
}
