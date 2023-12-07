############################################
# Functions to prepare data (specifically, generate the PCs) from a list or a single batch.
#' @importFrom BiocSingular ExactParam
#' @importFrom BiocParallel SerialParam
#' @importFrom irlba irlba
#' @importFrom BiocNeighbors findMutualNN
#' @importFrom methods as
#'
#' @title Find mutual nearest neighbors pairs
#' @param batch.list list of gene expression.
#' @param d dimension of PCA.
#' @param BSPARAM parallel.
#' @param BPPARAM parallel.
#' @return a sparse matrix - adjacent matrix
#' @export
mnn_adjacent <- function(batch.list, d=50, BSPARAM=ExactParam(), BPPARAM=SerialParam())
{
    nbatches <- length(batch.list)
    if (nbatches < 2L) {
        stop("at least two batches must be specified")
    }

    pc.mat <- lapply(batch.list, FUN = function(x) {
        svd_out <- irlba::irlba(x, nv = d)
        return(svd_out$u %*% diag(svd_out$d))
    })

    adjacent.mat <- matrix(0,
                           nrow = sum(sapply(batch.list, nrow)),
                           ncol = sum(sapply(batch.list, nrow)))

    dataPair <- data.frame(first.data = NULL, second.data = NULL)
    for(i in 1:(nbatches-1)){
        for(j in (i+1):nbatches) {
            dataPair <- rbind(dataPair, data.frame(first.data = i, second.data = j))
        }
    }

    batchNo <- c(0, cumsum(sapply(batch.list, nrow)))
    for(i in 1:nrow(dataPair)){
        mnnResult <- findMutualNN(data1 = pc.mat[[dataPair[i,1]]],
                                  data2 = pc.mat[[dataPair[i,2]]],
                                  k1 = 2)
        adjacent.mat[cbind(batchNo[dataPair[i,1]] + mnnResult$first,
                           batchNo[dataPair[i,2]] + mnnResult$second)] <- 1
        adjacent.mat[cbind(batchNo[dataPair[i,2]] + mnnResult$second,
                           batchNo[dataPair[i,1]] + mnnResult$first)] <- 1
    }

    adjacent.mat <- as(adjacent.mat, "sparseMatrix")
}

