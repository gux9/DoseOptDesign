#' Calculate Approximate Sample Size for ROSE Design
#'
#' ROSE (Response-Only Selection) design focuses solely on efficacy,
#' ignoring safety considerations. Equivalent to utility-based design
#' with utilities u = (1, 1, 0, 0).
#'
#' @param pL Baseline response rate for low dose
#' @param delta Efficacy margin
#' @param alpha_L Target PCS for Scenario L (default 0.8)
#' @param alpha_H Target PCS for Scenario H (default 0.8)
#'
#' @return List containing:
#'   - n: Sample size per arm
#'   - lambda_u: Decision threshold for response rate difference
#'   - PCS_L: Achieved PCS under Scenario L
#'   - PCS_H: Achieved PCS under Scenario H
#'   - utility: (1, 1, 0, 0)
#'   - scenario_L: Parameters for Scenario L
#'   - scenario_H: Parameters for Scenario H
#'   - method: "rose_approximate"
#'
#' @examples
#' # ROSE design
#' design <- calc_sample_size_rose_approx(
#'   pL = 0.3, delta = 0.15
#' )
#' print(paste("ROSE sample size:", design$n))
#' @export
calc_sample_size_rose_approx <- function(pL, delta,
                                         alpha_L = 0.8, alpha_H = 0.8) {
  
  # Input validation
  if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
    stop("pL must be between 0 and 1")
  }
  if (!is.numeric(delta) || delta <= 0) {
    stop("delta must be positive")
  }
  if (pL - delta < 0) {
    stop("pL - delta must be non-negative")
  }
  
  # Scenario L: Both doses have response rate pL (no difference)
  pL_L <- pL
  pH_L <- pL
  delta_mu_L <- 0
  v_L <- 2 * pL * (1 - pL)
  
  # Scenario H: High dose has response rate pL, low dose has pL - delta
  pL_H <- pL - delta
  pH_H <- pL 
  delta_mu_H <- delta
  v_H <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)  # Fixed: Use pL_H here
  
  # Calculate sample size
  z_alpha_L <- qnorm(alpha_L)
  z_1_alpha_H <- qnorm(1 - alpha_H)
  
  n_continuous <- ((z_alpha_L * sqrt(v_L) - z_1_alpha_H * sqrt(v_H)) / delta_mu_H)^2
  n <- ceiling(n_continuous)
  
  if (n < 2) {
    stop("Calculated sample size is less than 2. Check parameters.")
  }
  
  # Calculate decision threshold
  lambda_u <- delta_mu_H + z_1_alpha_H * sqrt(v_H / n)
  
  # Calculate achieved PCS
  pcs_L <- pnorm((lambda_u - delta_mu_L) / sqrt(v_L / n))
  pcs_H <- 1 - pnorm((lambda_u - delta_mu_H) / sqrt(v_H / n))
  
  # Return results
  result <- list(
    n = n,
    lambda_u = lambda_u,
    PCS_L = pcs_L,
    PCS_H = pcs_H,
    utility = c(1, 1, 0, 0),  # ROSE utilities
    scenario_L = list(
      delta_mu = delta_mu_L,
      variance = v_L,
      pL = pL_L, pH = pH_L
    ),
    scenario_H = list(
      delta_mu = delta_mu_H,
      variance = v_H,
      pL = pL_H, pH = pH_H
    ),
    method = "rose_approximate",
    inputs = list(
      alpha_L = alpha_L,
      alpha_H = alpha_H
    ),
    delta = delta
  )
  class(result) <- c("rose_design", "list")
  return(result)
}

#' Calculate Exact Sample Size for ROSE Design
#'
#' Uses exact binomial calculations for ROSE design.
#'
#' @param pL Baseline response rate for low dose
#' @param delta Efficacy margin
#' @param alpha_L Target PCS for Scenario L (default 0.8)
#' @param alpha_H Target PCS for Scenario H (default 0.8)
#' @param max_n Maximum sample size to search (default 200)
#' @param buffer Buffer below approximate n (default 10)
#' @param verbose Print progress (default FALSE)
#'
#' @return List with same structure as ROSE approximate, plus:
#'   - method: "rose_exact"
#'   - n_approx: Approximate starting point
#'
#' @examples
#' # ROSE exact calculation
#' design <- calc_sample_size_rose_exact(
#'   pL = 0.3, delta = 0.15
#' )
#' print(paste("ROSE exact sample size:", design$n))
#' @export
calc_sample_size_rose_exact <- function(pL, delta,
                                        alpha_L = 0.8, alpha_H = 0.8,
                                        max_n = 200,
                                        buffer = 10,
                                        verbose = FALSE) {
  
  # Input validation
  if (!is.numeric(pL) || pL <= 0 || pL >= 1) {
    stop("pL must be between 0 and 1")
  }
  if (!is.numeric(delta) || delta <= 0) {
    stop("delta must be positive")
  }
  if (pL - delta < 0) {
    stop("pL - delta must be non-negative")
  }
  
  # Define scenarios
  pL_L <- pL; pH_L <- pL
  pL_H <- pL - delta; pH_H <- pL
  
  # Get approximate starting point
  v_L <- 2 * pL * (1 - pL)
  v_H <- pL_H * (1 - pL_H) + pH_H * (1 - pH_H)
  
  z_alpha_L <- qnorm(alpha_L)
  z_1_alpha_H <- qnorm(1 - alpha_H)
  n_approx <- ceiling(((z_alpha_L * sqrt(v_L) - z_1_alpha_H * sqrt(v_H)) / delta)^2)
  
  min_search <- max(5, n_approx - buffer)
  
  if (verbose) {
    cat(sprintf("Approximate n = %d, searching [%d, %d]...\n", 
                n_approx, min_search, max_n))
  }
  
  # Exact search
  found <- FALSE
  n_final <- NA
  lambda_u_final <- NA
  pcs_L_final <- NA
  pcs_H_final <- NA
  
  for (n in min_search:max_n) {
    if (verbose && (n - min_search) %% 5 == 0) {
      cat(sprintf("  Testing n = %d...\n", n))
    }
    
    # Scenario L: compute PMF of difference
    pmf_L_SL <- compute_pmf_binomial(pL_L, n)
    pmf_H_SL <- compute_pmf_binomial(pH_L, n)
    
    # Difference PMF for binomial case
    pmf_diff_SL <- rep(0, 2*n + 1)
    offset <- n + 1
    
    for (x_L in 0:n) {
      for (x_H in 0:n) {
        d <- x_H - x_L
        pmf_diff_SL[d + offset] <- pmf_diff_SL[d + offset] + 
          pmf_L_SL[x_L + 1] * pmf_H_SL[x_H + 1]
      }
    }
    cdf_SL <- cumsum(pmf_diff_SL)
    
    # Scenario H: compute PMF of difference
    pmf_L_SH <- compute_pmf_binomial(pL_H, n)
    pmf_H_SH <- compute_pmf_binomial(pH_H, n)
    
    pmf_diff_SH <- rep(0, 2*n + 1)
    for (x_L in 0:n) {
      for (x_H in 0:n) {
        d <- x_H - x_L
        pmf_diff_SH[d + offset] <- pmf_diff_SH[d + offset] + 
          pmf_L_SH[x_L + 1] * pmf_H_SH[x_H + 1]
      }
    }
    cdf_SH <- cumsum(pmf_diff_SH)
    
    # Find quantiles (in count difference scale)
    q_SL <- min(which(cdf_SL >= alpha_L)) - offset
    q_SH <- max(which(cdf_SH <= 1 - alpha_H)) - offset
    
    # Check feasibility
    if (q_SL <= q_SH) {
      found <- TRUE
      n_final <- n
      
      # Threshold in count difference
      threshold_count <- round((q_SL + q_SH) / 2)
      # Convert to rate difference
      lambda_u_final <- threshold_count / n
      
      # Calculate exact PCS
      pcs_L_final <- cdf_SL[threshold_count + offset]
      pcs_H_final <- 1 - cdf_SH[threshold_count + offset]
      
      if (verbose) {
        cat(sprintf("  Solution found: n = %d, lambda_u = %.6f\n", 
                    n_final, lambda_u_final))
        cat(sprintf("    PCS_L = %.4f, PCS_H = %.4f\n", 
                    pcs_L_final, pcs_H_final))
      }
      
      break
    }
  }
  
  if (!found) {
    stop(sprintf("No solution found in range [%d, %d]. Increase max_n.", 
                 min_search, max_n))
  }
  
  # Return results
  result <- list(
    n = n_final,
    lambda_u = lambda_u_final,
    PCS_L = pcs_L_final,
    PCS_H = pcs_H_final,
    utility = c(1, 1, 0, 0),
    scenario_L = list(
      delta_mu = 0,
      variance = v_L,
      pL = pL_L, pH = pH_L
    ),
    scenario_H = list(
      delta_mu = delta,
      variance = v_H,
      pL = pL_H, pH = pH_H
    ),
    method = "rose_exact",
    n_approx = n_approx,
    inputs = list(
      alpha_L = alpha_L,
      alpha_H = alpha_H
    ),
    delta = delta
  )
  class(result) <- c("rose_design", "list")
  return(result)
}
