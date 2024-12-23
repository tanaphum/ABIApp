#' Run the standalone version
#'
#' @param options Named options that should be passed to the runApp call.
#' @return
#' Open browser
#' @export run_app_ABI
#'
#' @examples
#' run_app_ABI()

run_app_ABI <- function(options = list()) {
  app_dir <- system.file("ABI", package = "ABI")
  shiny::shinyAppDir(app_dir, options = options)
}
