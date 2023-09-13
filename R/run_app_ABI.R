#' Run the standalone version
#'
#' @param options Named options that should be passed to the runApp call.
#' @return
#' Open browser
#' @export
#'
#' @import shiny shinythemes stringr ggplot2 markdown

run_app_ABI <- function(options = list()) {
  app_dir <- system.file("ABI", package = "ABI")
  shiny::shinyAppDir(app_dir, options = options)
}
