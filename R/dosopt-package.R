#' dosopt: Utility Score Framework for Dose Optimization Studies
#'
#' @description
#' The \pkg{dosopt} package implements a frequentist framework for one-stage
#' randomized dose optimization studies, directly addressing the FDA's
#' \emph{Project Optimus} initiative for patient-centered dosing in oncology.
#'
#' @section Key functions:
#' \describe{
#'   \item{Utility specification}{
#'     \code{\link{calc_pi}}, \code{\link{calc_utility}},
#'     \code{\link{calc_utility_moments}}, \code{\link{phi_bounds}}
#'   }
#'   \item{Sample size — Utility-based design}{
#'     \code{\link{calc_sample_size_utility_approx}},
#'     \code{\link{calc_sample_size_utility_exact}}
#'   }
#'   \item{Sample size — ROSE design (efficacy-only)}{
#'     \code{\link{calc_sample_size_rose_approx}},
#'     \code{\link{calc_sample_size_rose_exact}}
#'   }
#'   \item{Sample size — Multi-scenario direct approach}{
#'     \code{\link{sample_size_direct_approx}},
#'     \code{\link{sample_size_direct_exact}}
#'   }
#'   \item{Bias characterization}{
#'     \code{\link{calc_bias}}, \code{\link{calc_bias_TTE}}
#'   }
#'   \item{Type I error}{
#'     \code{\link{calc_type1_error}},
#'     \code{\link{calc_type1_error_landmark}},
#'     \code{\link{calc_type1_error_exp}},
#'     \code{\link{calc_type1_error_cox}},
#'     \code{\link{calc_expected_events}}
#'   }
#'   \item{Monte Carlo simulation}{
#'     \code{\link{simulate_dose_opt}}
#'   }
#'   \item{Interactive Shiny app}{
#'     \code{\link{run_dosopt_app}}
#'   }
#' }
#'
#' @references
#' Gu X, Xu C, Xu L, Yuan Y (2026). A Utility Score Framework for Dose
#' Optimization Studies with Binary Efficacy-Safety Endpoints.
#'
#' Wang S, Yuan Y, Liu S (2025). Randomized Optimal Selection Design for
#' Dose Optimization. \emph{Biometrics}, 81(4):ujaf124.
#'
#' Shah M, et al. (2021). The drug-dosing conundrum in oncology.
#' \emph{New England Journal of Medicine}, 385(16):1445--1447.
#'
#' @keywords internal
"_PACKAGE"
