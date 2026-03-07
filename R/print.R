#' Print Method for ROSE Design Objects
#'
#' @param x A rose_design object from calc_sample_size_rose.
#' @param ... Additional arguments (not used).
#' @export
print.rose_design <- function(x, ...) {
  cat("ROSE Design (Response-Only Selection)\n")
  cat("======================================\n\n")
  
  cat("Sample Size:\n")
  cat(sprintf("  n per arm: %d\n", x$n))
  cat(sprintf("  Total n:   %d\n\n", 2 * x$n))
  
  cat("Decision Threshold:\n")
  cat(sprintf("  lambda_u:  %.6f\n\n", x$lambda_u))
  
  cat("Achieved PCS:\n")
  cat(sprintf("  PCS_L (select low dose):  %.4f (target: %.2f)\n", 
              x$PCS_L, x$inputs$alpha_L))
  cat(sprintf("  PCS_H (select high dose): %.4f (target: %.2f)\n\n", 
              x$PCS_H, x$inputs$alpha_H))
  
  cat("Utility Scores:\n")
  cat(sprintf("  u = [%.4f, %.4f, %.4f, %.4f]\n", 
              x$utility[1], x$utility[2], x$utility[3], x$utility[4]))
  cat("  (Response-only: ignores safety)\n\n")
  
  cat("Scenario L (Low dose should be selected):\n")
  cat(sprintf("  Dose L: p = %.3f\n", x$scenario_L$pL))
  cat(sprintf("  Dose H: p = %.3f\n", x$scenario_L$pH))
  cat(sprintf("  Delta_mu = %.4f, v = %.4f\n\n", 
              x$scenario_L$delta_mu, x$scenario_L$variance))
  
  cat("Scenario H (High dose should be selected):\n")
  cat(sprintf("  Dose L: p = %.3f\n", x$scenario_H$pL))
  cat(sprintf("  Dose H: p = %.3f\n", x$scenario_H$pH))
  cat(sprintf("  Delta_mu = %.4f, v = %.4f\n\n", 
              x$scenario_H$delta_mu, x$scenario_H$variance))
  
  cat("Efficacy margin: delta = ", x$delta, "\n\n")
  
  cat("Note: ROSE design considers only efficacy (response rate),\n")
  cat("      ignoring safety outcomes.\n")
  
  invisible(x)
}

#' Print Method for Dose Design Objects
#'
#' @param x A dose_design object from calc_sample_size or calc_sample_size_exact.
#' @param ... Additional arguments (not used).
#' @export
print.dose_design <- function(x, ...) {
  is_exact <- inherits(x, "dose_design_exact")
  method_label <- if (is_exact) "EXACT" else "APPROXIMATE"
  
  cat("Utility-Based Dose Selection Design (", method_label, ")\n", sep = "")
  cat(strrep("=", 50 + nchar(method_label)), "\n\n")
  
  cat("Method:\n")
  cat(sprintf("  Type: %s\n", method_label))
  cat("\n")
  
  cat("Sample Size:\n")
  cat(sprintf("  n per arm: %d\n", x$n))
  cat(sprintf("  Total n:   %d\n\n", 2 * x$n))
  
  cat("Decision Threshold:\n")
  cat(sprintf("  lambda_u:  %.6f\n\n", x$lambda_u))
  
  cat("Achieved PCS:\n")
  cat(sprintf("  PCS_L: %.4f (target: %.2f)\n", x$PCS_L, x$inputs$alpha_L))
  cat(sprintf("  PCS_H: %.4f (target: %.2f)\n\n", x$PCS_H, x$inputs$alpha_H))
  
  cat("Utility:\n")
  cat(sprintf("  u = [%.4f, %.4f, %.4f, %.4f]\n", 
              x$utility[1], x$utility[2], x$utility[3], x$utility[4]))
  cat("\n")
  
  cat("Scenario L (Low dose optimal):\n")
  cat(sprintf("  L: p=%.3f q=%.3f   H: p=%.3f q=%.3f\n",
              x$scenario_L$pL, x$scenario_L$qL, 
              x$scenario_L$pH, x$scenario_L$qH))
  cat(sprintf("  Delta_mu=%.4f, v=%.4f\n\n", 
              x$scenario_L$delta_mu, x$scenario_L$variance))
  
  cat("Scenario H (High dose optimal):\n")
  cat(sprintf("  L: p=%.3f q=%.3f   H: p=%.3f q=%.3f\n",
              x$scenario_H$pL, x$scenario_H$qL, 
              x$scenario_H$pH, x$scenario_H$qH))
  cat(sprintf("  Delta_mu=%.4f, v=%.4f\n\n", 
              x$scenario_H$delta_mu, x$scenario_H$variance))
  
  cat("Parameters:\n")
  cat(sprintf("  phi = %.2f\n", x$phi))
  if (!is.null(x$delta)) {
    cat(sprintf("  delta = %.3f, d = %.3f\n", x$delta, x$d))
  }
  
  if (is_exact) {
    cat("\nSearch Info:\n")
    cat(sprintf("  Approx n: %d\n", x$n_approx))
  }
  
  invisible(x)
}
