#' @title Print Methods for dosopt Design Objects
#' @description S3 print methods for output classes returned by sample size
#'   functions.
#' @name print_dosopt
NULL

#' @rdname print_dosopt
#' @param x A \code{dose_design} or \code{dose_design_exact} object.
#' @param ... Ignored.
#' @export
print.dose_design <- function(x, ...) {
  is_exact <- inherits(x, "dose_design_exact")
  label    <- if (is_exact) "EXACT (Multinomial PMF)" else "APPROXIMATE (Normal)"
  cat("\nUtility-Based Dose Selection Design —", label, "\n")
  cat(strrep("=", 55), "\n")
  cat(sprintf("  Sample size per arm : %d  (Total: %d)\n", x$n, 2L * x$n))
  cat(sprintf("  Decision threshold  : lambda_u = %.6f\n", x$lambda_u))
  cat(sprintf("  Achieved PCS_L      : %.4f  (target: %.2f)\n", x$PCS_L, x$inputs$alpha_L))
  cat(sprintf("  Achieved PCS_H      : %.4f  (target: %.2f)\n", x$PCS_H, x$inputs$alpha_H))
  cat(sprintf("  Utility scores      : u = (%.4f, %.4f, %.4f, %.4f)\n",
              x$utility[1], x$utility[2], x$utility[3], x$utility[4]))
  cat(sprintf("  Correlation (phi)   : %.2f\n", x$inputs$phi))
  if (!is.null(x$delta) && !is.na(x$delta))
    cat(sprintf("  Margins             : delta = %.3f, d = %.3f\n", x$delta, x$d))
  if (is_exact)
    cat(sprintf("  Approx. n (normal)  : %d\n", x$search_info$n_approx))
  cat("\n  Scenarios:\n")
  cat(sprintf("    S_L: DoseL=(%.3f,%.3f), DoseH=(%.3f,%.3f), Delta_mu=%.4f\n",
              x$scenario_L$pL, x$scenario_L$qL,
              x$scenario_L$pH, x$scenario_L$qH,
              x$scenario_L$delta_mu))
  cat(sprintf("    S_H: DoseL=(%.3f,%.3f), DoseH=(%.3f,%.3f), Delta_mu=%.4f\n",
              x$scenario_H$pL, x$scenario_H$qL,
              x$scenario_H$pH, x$scenario_H$qH,
              x$scenario_H$delta_mu))
  cat("\n")
  invisible(x)
}

#' @rdname print_dosopt
#' @export
print.rose_design <- function(x, ...) {
  label <- if (x$method == "rose_exact") "EXACT (Binomial PMF)" else "APPROXIMATE (Normal)"
  cat("\nROSE Design (Efficacy-Only Selection) —", label, "\n")
  cat(strrep("=", 55), "\n")
  cat(sprintf("  Sample size per arm : %d  (Total: %d)\n", x$n, 2L * x$n))
  cat(sprintf("  Decision threshold  : lambda_u = %.6f\n", x$lambda_u))
  cat(sprintf("  Achieved PCS_L      : %.4f  (target: %.2f)\n", x$PCS_L, x$inputs$alpha_L))
  cat(sprintf("  Achieved PCS_H      : %.4f  (target: %.2f)\n", x$PCS_H, x$inputs$alpha_H))
  cat(sprintf("  Efficacy margin     : delta = %.3f\n", x$delta))
  if (!is.null(x$n_approx))
    cat(sprintf("  Approx. n (normal)  : %d\n", x$n_approx))
  cat("\n")
  invisible(x)
}

#' @rdname print_dosopt
#' @export
print.ss_direct <- function(x, ...) {
  cat(sprintf("\nDose Optimization — Direct Approach (%s method)\n", x$method))
  cat(strrep("=", 55), "\n")
  cat(sprintf("  Design size (per arm) : n = %d\n", x$n_design))
  lam <- if (length(x$lambda_u) == 1) x$lambda_u else x$scenarios$lambda_u[1]
  cat(sprintf("  Decision threshold    : lambda_u = %.4f%s\n", lam,
              if (length(x$lambda_u) > 1) " (first scenario)" else ""))
  cat(sprintf("  Utility scores        : u = (%s)\n",
              paste(round(x$utility, 4), collapse = ", ")))
  cat(sprintf("  Target PCS            : alpha_L = %.2f, alpha_H = %.2f\n",
              x$alpha_L, x$alpha_H))
  if (x$method == "exact")
    cat(sprintf("  Integer scaling       : den = %d, u_int = (%s)\n",
                x$den, paste(x$u_int, collapse = ", ")))

  df <- x$scenarios
  n_scen <- nrow(df)
  cat(sprintf("\n  Per-Scenario Results (%d scenario%s):\n",
              n_scen, if (n_scen > 1) "s" else ""))

  for (i in seq_len(n_scen)) {
    s <- df[i, ]
    cat(sprintf("  [%d] p=%.2f, q=%.2f, delta=%.2f, d=%.2f, phi=%.2f\n",
                i, s$p, s$q, s$delta, s$d, s$phi))
    if (x$method == "approximate") {
      cat(sprintf("      n_L=%d, n_H=%d => n=%d [%s]  PCS_L=%.4f, PCS_H=%.4f\n",
                  s$n_L, s$n_H, s$n_scenario, s$binding, s$PCS_L, s$PCS_H))
    } else {
      cat(sprintf("      n_approx=%d, n_exact=%d [%s]  PCS_L=%.4f, PCS_H=%.4f\n",
                  s$n_approx, s$n_exact, s$binding, s$PCS_L_exact, s$PCS_H_exact))
    }
  }
  cat("\n")
  invisible(x)
}
