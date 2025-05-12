#' Load LDA model for package
#'
#' @param model_name a string corresponding to the name of the model
#' @return A fitted LDA model (object of class \code{lda})
loadLDA <- function(model_name) {
  model_path <- system.file("models", paste0(model_name, ".rds"), package = "WaterStrideR")
  readRDS(model_path)
}