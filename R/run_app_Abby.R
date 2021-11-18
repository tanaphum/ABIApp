#' Run the standalone version
#'
#' @param options Named options that should be passed to the runApp call.
#' @return
#' Open browser
#' @export
#'
#' @import shiny shinythemes stringr ggplot2 DT markdown

run_app_Abby <- function(options = list()) {
  app_dir <- system.file("Abbyapp", package = "Abby")
  shiny::shinyAppDir(app_dir, options = options)
}
