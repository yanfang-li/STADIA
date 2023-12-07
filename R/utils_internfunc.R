# Initialize subtype parameter c in function ParameterInitialize()
#' @importFrom mclust Mclust mclustBIC
.InitialC <- function(Y, q = 7) {
  init <- list()
  clust_result <- mclust::Mclust(data = Y, G = q, modelNames = "EEE", verbose = FALSE)
  init$c <- clust_result$classification
  init$mu <- clust_result$parameters$mean
  init$Lambda <- clust_result$parameters$variance$Sigma

  return(init)
}

# batch vector to covariable matrix M
batchVec2M <- function(batch_vec) {
  n_b <- length(unique(batch_vec))
  n <- length(batch_vec)
  result <- matrix(0, nrow = n, ncol = n_b)
  for (i in 1:n_b) {
    result[batch_vec == i, i] <- 1
  }
  return(result)
}

#' # calculate smoothing parameter in potts model
#' #' @importFrom lctools moransI
#' #' @importFrom dplyr  %>%
#' #' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
#' calEta <- function(object) {
#'   position <- cbind(object$row, object$col)
#'
#'   object <- object %>%
#'     Seurat::NormalizeData(verbose = FALSE) %>%
#'     Seurat::FindVariableFeatures(verbose = FALSE) %>%
#'     Seurat::ScaleData(verbose = FALSE) %>%
#'     Seurat::RunPCA(dims = 1, verbose = FALSE)
#'
#'   m.I <- lctools::moransI(
#'     position, 20,
#'     object@reductions$pca@cell.embeddings[, 1]
#'   )
#'
#'   eta <- (m.I$Morans.I * 2)^{
#'     1.5
#'   }
#'   return(eta)
#' }
