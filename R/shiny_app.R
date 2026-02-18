#' Launch the Dose Optimization Shiny Application
#'
#' Opens an interactive browser-based application for exploring utility-based
#' dose optimization designs, ROSE designs, bias characterization, and
#' Monte Carlo simulation validation.
#'
#' @param launch.browser Logical. Open in browser? Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @details
#' The Shiny app provides four tabs:
#' \describe{
#'   \item{General Utility Score Design}{Sample size calculation with four
#'     specification modes: margin-based, partial-direct, full-direct, and
#'     direct approach (multi-scenario).}
#'   \item{ROSE Design}{Efficacy-only selection sample size.}
#'   \item{Simulation Analysis}{Monte Carlo validation of bias and Type I
#'     error formulas, with optional TTE analysis.}
#'   \item{Formulas & Methods}{Mathematical reference for all key equations.}
#' }
#'
#' Requires the Suggests packages: \pkg{shiny}, \pkg{shinythemes},
#' \pkg{DT}, \pkg{shinyBS}.
#'
#' @examples
#' \dontrun{
#' run_dosopt_app()
#' }
#'
#' @export
run_dosopt_app <- function(launch.browser = TRUE, ...) {
  pkg_needed <- c("shiny", "shinythemes", "DT", "shinyBS")
  missing_pkg <- pkg_needed[!vapply(pkg_needed, requireNamespace,
                                     logical(1), quietly = TRUE)]
  if (length(missing_pkg) > 0) {
    stop(sprintf(
      "The following packages are needed to run the Shiny app:\n%s\n",
      paste(sprintf("  install.packages('%s')", missing_pkg), collapse = "\n")
    ))
  }
  app_dir <- system.file("shiny", package = "dosopt")
  if (!nzchar(app_dir) || !file.exists(app_dir)) {
    stop("Shiny app directory not found in the package installation.")
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
