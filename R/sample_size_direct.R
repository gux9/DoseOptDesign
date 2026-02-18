#' Sample Size for Utility-Based Design — Direct Approach (Approximate)
#'
#' Computes per-arm sample sizes for multiple \eqn{(p, q)} scenarios
#' simultaneously (Mode 4 / Direct Approach), as in Equations 11-12 of the
#' manuscript. Returns the maximum across all scenarios as the design size.
#'
#' @param p Numeric vector. Response rates for each scenario.
#' @param q Numeric vector. No-adverse-event rates for each scenario.
#' @param delta Numeric. Common efficacy margin.
#' @param d Numeric. Common safety margin.
#' @param phi Numeric scalar or vector. Correlation(s) (recycled). Default 0.
#' @param alpha_L Numeric. Target PCS under \eqn{S_L}. Default 0.8.
#' @param alpha_H Numeric. Target PCS under \eqn{S_H}. Default 0.8.
#' @param lambda_u Numeric or \code{NULL}. Decision threshold. When \code{NULL}
#'   (default), the jointly optimal threshold is computed per scenario.
#' @param r Numeric or \code{NULL}. Trade-off ratio; overrides \code{delta/d}
#'   for utility computation.
#' @param u Numeric vector of length 4 or \code{NULL}. Custom utility scores.
#'
#' @return An object of class \code{"ss_direct"} with elements:
#'   \describe{
#'     \item{n_design}{Design sample size (per arm) = maximum across scenarios.}
#'     \item{lambda_u}{Decision threshold vector (one per scenario).}
#'     \item{scenarios}{Data frame with per-scenario details.}
#'     \item{utility, alpha_L, alpha_H, method}{Design metadata.}
#'   }
#'
#' @examples
#' res <- sample_size_direct_approx(
#'   p = c(0.3, 0.4, 0.5), q = c(0.5, 0.5, 0.5),
#'   delta = 0.10, d = 0.10,
#'   phi = 0, alpha_L = 0.8, alpha_H = 0.8
#' )
#' print(res)
#'
#' @seealso \code{\link{sample_size_direct_exact}}
#' @export
sample_size_direct_approx <- function(p, q, delta, d,
                                       phi = 0,
                                       alpha_L = 0.8, alpha_H = 0.8,
                                       lambda_u = NULL,
                                       r = NULL, u = NULL) {
  n_scen <- max(length(p), length(q))
  p   <- rep_len(p,   n_scen); q   <- rep_len(q,   n_scen)
  phi <- rep_len(phi, n_scen)

  if (is.null(u)) {
    if (is.null(r)) r <- delta / d
    u <- calc_utility(r)
  }
  if (!is.numeric(u) || length(u) != 4) stop("'u' must be length-4 numeric.")

  z_aL  <- qnorm(alpha_L); z_aH <- qnorm(alpha_H)
  optimize_lambda <- is.null(lambda_u)

  results <- data.frame(
    p = p, q = q, delta = delta, d = d, phi = phi,
    Delta_mu_L = NA_real_, v_L = NA_real_, n_L = NA_integer_,
    Delta_mu_H = NA_real_, v_H = NA_real_, n_H = NA_integer_,
    n_scenario = NA_integer_, binding = NA_character_,
    lambda_u = NA_real_, PCS_L = NA_real_, PCS_H = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_scen)) {
    pi_L_SL <- calc_pi(p[i], q[i], phi[i])
    pi_H_SL <- calc_pi(p[i], q[i] - d, phi[i])
    pi_L_SH <- calc_pi(p[i] - delta, q[i], phi[i])
    pi_H_SH <- calc_pi(p[i], q[i], phi[i])

    mom_L_SL <- calc_utility_moments(pi_L_SL, u)
    mom_H_SL <- calc_utility_moments(pi_H_SL, u)
    mom_L_SH <- calc_utility_moments(pi_L_SH, u)
    mom_H_SH <- calc_utility_moments(pi_H_SH, u)

    DL <- mom_H_SL$mu - mom_L_SL$mu; vL <- mom_L_SL$sigma2 + mom_H_SL$sigma2
    DH <- mom_H_SH$mu - mom_L_SH$mu; vH <- mom_L_SH$sigma2 + mom_H_SH$sigma2

    lam_i <- if (optimize_lambda) {
      a <- z_aL * sqrt(vL); b <- -qnorm(1 - alpha_H) * sqrt(vH)
      (DL * b + DH * a) / (a + b)
    } else lambda_u

    if (lam_i <= DL || lam_i >= DH)
      stop(sprintf("Scenario %d: lambda_u=%.4f not in feasible range (%.4f, %.4f).",
                   i, lam_i, DL, DH))

    nL_i <- ceiling(z_aL^2 * vL / (lam_i - DL)^2)
    nH_i <- ceiling(z_aH^2 * vH / (DH - lam_i)^2)
    n_i  <- max(nL_i, nH_i)

    results$Delta_mu_L[i] <- DL; results$v_L[i]    <- vL; results$n_L[i] <- nL_i
    results$Delta_mu_H[i] <- DH; results$v_H[i]    <- vH; results$n_H[i] <- nH_i
    results$n_scenario[i] <- n_i
    results$binding[i]    <- if (nL_i > nH_i) "S_L" else if (nH_i > nL_i) "S_H" else "both"
    results$lambda_u[i]   <- lam_i
    results$PCS_L[i]      <- pnorm((lam_i - DL) / sqrt(vL / n_i))
    results$PCS_H[i]      <- 1 - pnorm((lam_i - DH) / sqrt(vH / n_i))
  }

  n_design <- max(results$n_scenario)
  for (i in seq_len(n_scen)) {
    lam_i <- results$lambda_u[i]
    results$PCS_L[i] <- pnorm((lam_i - results$Delta_mu_L[i]) / sqrt(results$v_L[i] / n_design))
    results$PCS_H[i] <- 1 - pnorm((lam_i - results$Delta_mu_H[i]) / sqrt(results$v_H[i] / n_design))
  }

  structure(
    list(n_design = n_design,
         lambda_u = results$lambda_u,
         scenarios = results, utility = u,
         alpha_L = alpha_L, alpha_H = alpha_H, method = "approximate"),
    class = "ss_direct"
  )
}


#' Sample Size for Utility-Based Design — Direct Approach (Exact)
#'
#' Exact version of \code{\link{sample_size_direct_approx}} using dynamic
#' programming on the Multinomial PMF.
#'
#' @inheritParams sample_size_direct_approx
#' @param den Integer. Denominator for integer-scaling of utility scores.
#'   Default 10.
#' @param buffer Integer. Search buffer below approximate solution. Default 10.
#' @param max_n Integer. Upper bound for sample size search. Default 500.
#'
#' @return An object of class \code{"ss_direct"} with the same structure as
#'   \code{\link{sample_size_direct_approx}}, plus exact PCS values and
#'   the integer-scaling metadata.
#'
#' @examples
#' res <- sample_size_direct_exact(
#'   p = c(0.3, 0.4), q = c(0.5, 0.5),
#'   delta = 0.10, d = 0.10,
#'   phi = 0, alpha_L = 0.8, alpha_H = 0.8,
#'   max_n = 200
#' )
#' print(res)
#'
#' @export
sample_size_direct_exact <- function(p, q, delta, d,
                                      phi = 0,
                                      alpha_L = 0.8, alpha_H = 0.8,
                                      lambda_u = NULL,
                                      r = NULL, u = NULL,
                                      den = 10, buffer = 10, max_n = 500) {
  approx <- sample_size_direct_approx(
    p = p, q = q, delta = delta, d = d,
    phi = phi, alpha_L = alpha_L, alpha_H = alpha_H,
    lambda_u = lambda_u, r = r, u = u
  )
  u_used <- approx$utility
  u_int  <- as.integer(round(u_used * den))
  if (all(u_int == 0)) stop("Integer-scaled utilities are all zero. Increase 'den'.")

  df    <- approx$scenarios
  n_scen <- nrow(df)
  df$n_approx  <- df$n_scenario
  df$n_exact   <- NA_integer_
  df$PCS_L_exact <- NA_real_
  df$PCS_H_exact <- NA_real_
  phi_vec <- df$phi; p_vec <- df$p; q_vec <- df$q

  for (i in seq_len(n_scen)) {
    lam_i   <- df$lambda_u[i]
    n_start <- max(5, df$n_approx[i] - buffer)

    pi_L_SL <- calc_pi(p_vec[i], q_vec[i], phi_vec[i])
    pi_H_SL <- calc_pi(p_vec[i], q_vec[i] - d, phi_vec[i])
    pi_L_SH <- calc_pi(p_vec[i] - delta, q_vec[i], phi_vec[i])
    pi_H_SH <- calc_pi(p_vec[i], q_vec[i], phi_vec[i])

    found <- FALSE
    for (n_try in n_start:max_n) {
      thr_int <- round(n_try * lam_i * den)
      pcs_L   <- .exact_pcs(pi_L_SL, pi_H_SL, u_int, n_try, thr_int, "L")
      pcs_H   <- .exact_pcs(pi_L_SH, pi_H_SH, u_int, n_try, thr_int, "H")
      if (pcs_L >= alpha_L && pcs_H >= alpha_H) {
        df$n_exact[i]    <- n_try
        df$PCS_L_exact[i] <- pcs_L
        df$PCS_H_exact[i] <- pcs_H
        found <- TRUE; break
      }
    }
    if (!found) {
      warning(sprintf("Scenario %d: no exact solution for n <= %d.", i, max_n))
      df$n_exact[i] <- max_n
    }
  }

  n_design <- max(df$n_exact)
  for (i in seq_len(n_scen)) {
    pi_L_SL <- calc_pi(p_vec[i], q_vec[i], phi_vec[i])
    pi_H_SL <- calc_pi(p_vec[i], q_vec[i] - d, phi_vec[i])
    pi_L_SH <- calc_pi(p_vec[i] - delta, q_vec[i], phi_vec[i])
    pi_H_SH <- calc_pi(p_vec[i], q_vec[i], phi_vec[i])
    thr_int <- round(n_design * df$lambda_u[i] * den)
    df$PCS_L_exact[i] <- .exact_pcs(pi_L_SL, pi_H_SL, u_int, n_design, thr_int, "L")
    df$PCS_H_exact[i] <- .exact_pcs(pi_L_SH, pi_H_SH, u_int, n_design, thr_int, "H")
  }

  df$binding <- ifelse(df$PCS_L_exact - alpha_L < df$PCS_H_exact - alpha_H, "S_L",
                       ifelse(df$PCS_H_exact - alpha_H < df$PCS_L_exact - alpha_L, "S_H", "both"))

  structure(
    list(n_design = n_design,
         lambda_u = if (length(unique(df$lambda_u)) == 1) df$lambda_u[1] else df$lambda_u,
         scenarios = df, utility = u_used, u_int = u_int, den = den,
         alpha_L = alpha_L, alpha_H = alpha_H, method = "exact"),
    class = "ss_direct"
  )
}

#' @keywords internal
.exact_pcs <- function(pi_L, pi_H, u_int, n, threshold_int, direction = "L") {
  pmf_L <- .compute_pmf_S(pi_L, u_int, n)
  pmf_H <- .compute_pmf_S(pi_H, u_int, n)
  max_S  <- n * max(u_int); offset <- max_S + 1
  pmf_d  <- .compute_pmf_diff(pmf_L, pmf_H, max_S, "nested")
  idx    <- threshold_int + offset
  if (idx < 1)                  return(if (direction == "L") 0 else 1)
  if (idx >= length(pmf_d) + 1) return(if (direction == "L") 1 else 0)
  cdf_at <- sum(pmf_d[seq_len(idx)])
  if (direction == "L") cdf_at else 1 - cdf_at
}
