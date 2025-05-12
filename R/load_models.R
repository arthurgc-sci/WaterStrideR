#' Load LDA model for package
#'
#' @param model_name a string corresponding to the name of the model
#' 
loadLDA <- function(model_name) {
  model_path <- system.file("models", paste0(model_name, ".rds"), package = "WaterStrideR")
  readRDS(model_path)
}