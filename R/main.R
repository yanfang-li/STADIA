#' @title Core function to conduct STADIA algorithm
#' @description core function to conduct STADIA algorithm
#'
#' @importFrom methods slot
#'
#' @param object.list a list of seurat object.
#' @param hyper a list of all hyperParameters.
#' @param dim dimension of hidden space.
#' @param n_cluster number of spatial domains.
#' @param platform string vector to idicate platform generating data
#' @param adj.cutoff cutoff when calculating adjacent matrix.
#' @param icm.maxiter number of iteration of icm.
#' @param em.maxiter number of iteration of em.
#' @param min.features filtering low-quality data.
#' @param min.spots filtering low-quality data.
#' @param nfeatures number of HVGs.
#' @param verbose print information or not.
#' @param verbose.in print information or not.
#' @return a list of estimated parameters.
#'
#' @export
stadia <- function(object.list, hyper, dim = 35, n_cluster = 7,
                  platform = c("visium", "st", "others"),
                  adj.cutoff = 50, icm.maxiter = 10, em.maxiter = 30,
                  min.features = 200, min.spots = 20, nfeatures = 2000,
                  verbose = TRUE, verbose.in = TRUE) {
    # dimensions
    d <- dim
    K <- n_cluster

    ##### preprocess data
    cat("Preprocess data...\n")
    object.list <- PreprocessData(object.list, min.features, min.spots, nfeatures)

    cat("Set initialization...\n")
    init <- InitiParameters(object.list, K, d)

    ##### utilities
    batch_vec <- rep(seq_along(object.list), times = sapply(object.list, ncol))
    ncores <- parallel::detectCores()

    ######  main
    ## extract coordinates and expression matrix
    position <- Reduce(rbind, lapply(object.list, function(x) {
        cbind(x$row, x$col)
    })) # nx2
    X.list <- lapply(object.list, function(x) {
        Seurat::GetAssayData(x, slot = "scale.data")
    })
    batch.list <- lapply(X.list, t)
    X <- do.call(cbind, X.list) # pxn

    ## run model
    if (!is.null(platform) && substr(tolower(platform),1,1) %in% c("v", "s")) {
        platform <- tolower(platform)
        platform <- match.arg(platform)
    } else {
        platform <- "others"
    }

    # adjacent matrix for mnn
    adj_mat_mnn <- mnn_adjacent(batch.list)
    out <- stadia_EM_SP(
        X, position, batch_vec, adj_mat_mnn,
        hyper, init, platform,
        d, K, adj.cutoff, icm.maxiter, em.maxiter, verbose,
        verbose.in, ncores
    )
    rownames(out$L) <- rownames(X)
    colnames(out$factors) <- colnames(X)

    return(out)
}
