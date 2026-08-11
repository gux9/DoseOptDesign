#' Attained Probability of Correct Selection for a Given Design
#'
#' Computes the exact probability of correct selection actually attained by a
#' design `(n, lambda_u)`, by direct multinomial enumeration of the distribution
#' of `Ubar_H - Ubar_L` under each scenario.
#'
#' This is the quantity a design delivers in practice. It differs from the
#' normal-approximation PCS returned by [calc_sample_size_utility_approx()] and
#' [calc_sample_size_rose_approx()], which equals the nominal target by
#' construction because the sample size was solved from the same normal
#' equations. At the small sample sizes typical of dose optimization studies the
#' two can differ materially: the utility score lives on a lattice of spacing
#' `1/(den * n)`, and when the optimal threshold falls inside the first lattice
#' cell the design behaves as if `lambda_u = 0`.
#'
#' @param pi_L_SL,pi_H_SL Outcome probability vectors (length 4) for the low and
#'   high dose under Scenario L.
#' @param pi_L_SH,pi_H_SH Outcome probability vectors (length 4) for the low and
#'   high dose under Scenario H.
#' @param u Utility scores (length 4).
#' @param n Sample size per arm.
#' @param lambda_u Selection threshold; Dose H is selected when
#'   `Ubar_H - Ubar_L > lambda_u`.
#' @param den Denominator used to place the utility scores on an integer
#'   lattice, as in [calc_sample_size_utility_exact()]. Default 10.
#'
#' @return A list with `PCS_L` and `PCS_H`, the exact attained probabilities of
#'   correct selection under Scenarios L and H. Both are `NA` (with a warning)
#'   if `u * den` is not integral, since the lattice representation is then
#'   inexact.
#'
#' @examples
#' pi_L_SL <- calc_pi(0.3, 0.5, 0)
#' pi_H_SL <- calc_pi(0.3, 0.35, 0)
#' pi_L_SH <- calc_pi(0.2, 0.5, 0)
#' pi_H_SH <- calc_pi(0.3, 0.5, 0)
#' calc_attained_pcs(pi_L_SL, pi_H_SL, pi_L_SH, pi_H_SH,
#'                   u = c(1, 0.6, 0.4, 0), n = 44, lambda_u = 0.00142)
#'
#' @export
calc_attained_pcs <- function(pi_L_SL, pi_H_SL, pi_L_SH, pi_H_SH,
                              u, n, lambda_u, den = 10) {

  if (length(u) != 4) stop("u must be length 4")
  if (!is.numeric(n) || n < 1) stop("n must be a positive integer")
  n <- as.integer(n)

  # ROSE (u = (1, 1, 0, 0)) needs no scaling; keeps the lattice small.
  if (isTRUE(all.equal(as.numeric(u), c(1, 1, 0, 0)))) den <- 1

  u_int <- round(u * den)
  if (max(abs(u * den - u_int)) > 1e-9) {
    warning("u * den is not integral; the attained PCS cannot be computed on ",
            "this lattice. Increase 'den' (e.g. den = 100) or supply utility ",
            "scores with at most log10(den) decimal places.")
    return(list(PCS_L = NA_real_, PCS_H = NA_real_))
  }

  max_S <- n * max(u_int)

  pmf_diff_SL <- compute_pmf_diff_unified(
    compute_pmf_S(pi_L_SL, u_int, n), compute_pmf_S(pi_H_SL, u_int, n),
    max_S = max_S, method = "nested", validate = FALSE)
  pmf_diff_SH <- compute_pmf_diff_unified(
    compute_pmf_S(pi_L_SH, u_int, n), compute_pmf_S(pi_H_SH, u_int, n),
    max_S = max_S, method = "nested", validate = FALSE)

  cdf_SL <- cumsum(pmf_diff_SL)
  cdf_SH <- cumsum(pmf_diff_SH)
  offset <- max_S + 1

  # The difference S_H - S_L is supported on the integers -max_S .. max_S.
  # Selecting H when Ubar_H - Ubar_L > lambda_u is the event
  # S_H - S_L > lambda_u * den * n, i.e. S_H - S_L > floor(lambda_u * den * n).
  lambda_int <- floor(lambda_u * den * n + 1e-9)

  if (lambda_int < -max_S) {
    return(list(PCS_L = 0, PCS_H = 1))
  }
  if (lambda_int >= max_S) {
    return(list(PCS_L = 1, PCS_H = 0))
  }

  idx <- lambda_int + offset
  list(PCS_L = unname(cdf_SL[idx]), PCS_H = unname(1 - cdf_SH[idx]))
}
