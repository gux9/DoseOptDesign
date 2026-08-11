#' Launch the DoseOptDesign Shiny Application
#'
#' Opens an interactive Shiny app for utility score-based dose optimization
#' design, including sample size calculation, bias analysis, and simulation.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @examples
#' \dontrun{
#' run_app()
#' }
#'
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install it with install.packages('shiny').")
  }
  app_dir <- system.file("shiny", package = "DoseOptDesign")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try reinstalling 'DoseOptDesign'.")
  }
  shiny::runApp(app_dir, ...)
}
