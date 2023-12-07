#' @include utils_internfunc.R
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Functions for preprocessing data
#'
#' @importFrom stats lm qnorm rbinom varimax
#' @importFrom Seurat NormalizeData FindVariableFeatures
#' SelectIntegrationFeatures ScaleData GetAssayData DefaultAssay
#' @importFrom Matrix rowSums
#' @importFrom irlba irlba
#' @importFrom mombf qmom
#' @importFrom progress progress_bar
#' @param object.list A list of Seurat objects.
#' @param min.features Filtering spots no more than min.features expression gene.
#' @param min.spots Filter genes no more than min.spots expression spots.
#' @param nfeatures Number of high variable features.
#' @return a list of Seurat objects
#' @export
PreprocessData <- function(object.list, min.features = 200, min.spots = 20,
                           nfeatures = 2000) {
    pb <- progress_bar$new(total = length(object.list)*4, format = "[:bar] :percent :elapsed")
    ##### filtering low-quality genes and spots
    object.list <- lapply(object.list, FUN = function(x) {
        pb$tick()
        suppressWarnings(
            # filtering spots that less than min.features expressed therein
            x <- x[, x@meta.data[, paste0("nFeature_", DefaultAssay(x))] > min.features],
            classes = 'validationWarning'
        )
        suppressWarnings(
            # filtering genes that less than min.spots expressed therein
            x <- x[Matrix::rowSums(x[[DefaultAssay(x)]]@counts > 0) > min.spots, ],
            classes = 'validationWarning'
        )
        return(x)
    })

    #### select features
    ## select HVGs on each batch
    object.list <- lapply(X = object.list, FUN = function(x) {
        pb$tick()
        x <- Seurat::FindVariableFeatures(x, nfeatures = nfeatures, verbose = FALSE)
        return(x)
    })
    ## select joint HVGs among all batches
    features <- Seurat::SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeatures)

    #### library-size normalization and log normalization: cell-dimension to @data
    object.list <- lapply(X = object.list, FUN = function(x) {
        pb$tick()
        x <- Seurat::NormalizeData(x,
                                   verbose = FALSE,
                                   normalization.method = "LogNormalize"
        )
        return(x)
    })

    ##### center without scale data: gene-dimension to @scale.data
    object.list <- lapply(X = object.list, FUN = function(x) {
        pb$tick()
        ## center data on joint HVGs
        x <- Seurat::ScaleData(x,
                               features = features, verbose = FALSE,
                               do.scale = FALSE, do.center = TRUE
        )
        suppressWarnings(
            x <- x[features, ],
            classes = 'validationWarning'
        )
        return(x)
    })

    return(object.list)
}


#' @title Set hyperparamters used in the model
#'
#' @param object.list list of seurat object
#' @param d dimension of latent factor.
#' @param eta smoothing parameters in Potts model.
#' @param mu_mu mean vector of gaussian distribution (mu).
#' @param Sigma_mu shared covariance matrix of gaussian distribution (mu).
#' @param nu_omega shape and rate parameters of gamma distribution (omega).
#' @param n_Lambda degree of freedom of wishart distribution (Lambda).
#' @param Sigma_Lambda pds location matrix (Lambda).
#' @param nu_tau shape and rate parameters of gamma distribution (tau).
#' @param lambda0 variance of gaussian distribution (L).
#' @param lambda1 variance of gaussian distribution (L).
#' @param alpha_p shape1 parameter of beta distribution (p).
#' @param beta_p shape2 parameter of beta distribution (p).
#' @return result: a list of hyperparameters
#' @export
HyperParameters <- function(object.list, d = 35, eta = NULL, mu_mu = NULL, Sigma_mu = NULL,
                            nu_omega = NULL, n_Lambda = NULL, Sigma_Lambda = NULL,
                            nu_tau = NULL, lambda0 = NULL, lambda1 = NULL,
                            alpha_p = NULL, beta_p = NULL) {
    result <- list()

    # parameters in Potts model
    batchNo <- length(object.list)
    if (length(eta) == 1) {
        result$eta <- rep(eta, batchNo)
    } else if (length(eta) == batchNo) {
        result$eta <- eta
    } else {
        warning("length of eta given is not correct, only the first is used!")
        result$eta <- rep(eta[1], batchNo)
    }
    # result$eta <- ifelse(!is.null(eta), eta, 1.5)

    # hyperparameters in layer 2
    result$mu_mu <- if (!is.null(mu_mu)) mu_mu else rep(0, d)
    result$Sigma_mu <- if (!is.null(Sigma_mu)) Sigma_mu else 100 * diag(d) ## check
    result$nu_omega <- ifelse(!is.null(nu_omega), nu_omega, 2)
    result$n_Lambda <- ifelse(!is.null(n_Lambda), n_Lambda, d)
    result$Sigma_Lambda <- if (!is.null(Sigma_Lambda)) Sigma_Lambda else 100 * diag(d) ## check

    # hyperparameters in layer 1
    result$nu_tau <- ifelse(!is.null(nu_tau), nu_tau, 1)
    aus_p <- 0.99 ## check
    # aus_cutoff <- 1e-2
    aus_cutoff <- 0.1
    result$lambda0 <- ifelse(!is.null(lambda0), lambda0, aus_cutoff / qnorm((1 - aus_p) / 2)^2)
    # library(mombf)
    result$lambda1 <- ifelse(!is.null(lambda1), lambda1, aus_cutoff / mombf::qmom(aus_p / 2)^2)
    result$alpha_p <- ifelse(!is.null(alpha_p), alpha_p, 1)
    result$beta_p <- ifelse(!is.null(beta_p), beta_p, 1)

    return(result)
}


#' @title Initializing parameters used in EM algorithm later
#'
#' @importFrom progress progress_bar
#' @param object.list a list of interested seurat object (preprocessed).
#' @param k number of subtypes.
#' @param d dimension of latent space.
#' @return a list of initialized parameters.
#' @export
InitiParameters <- function(object.list, k = 7, d = 35) {
    pb <- progress_bar$new(total = 5, format = "[:bar] :percent :elapsed")
    ###### allocate space for result
    result <- list()

    ## batch vector
    batch_vec <- rep(seq_along(object.list), times = sapply(object.list, ncol))

    ###### dimensions
    n <- sum(sapply(object.list, ncol))
    if (length(unique(sapply(object.list, nrow))) != 1) {
        stop("ParameterInitialize should be conducted after PreprocessData()")
    }
    p <- sapply(object.list, nrow)[1]
    b <- length(object.list)

    ##### combined all (normalized) data
    if (b == 1) {
        data2use <- t(GetAssayData(object.list[[1]], slot = "scale.data"))
    } else {
        data_use_list <- sapply(object.list, function(object) {
            Seurat::GetAssayData(object, slot = "scale.data")
        })
        data2use <- t(Reduce(cbind, data_use_list)) # nxp
    }
    pb$tick()

    ###### initialization for L, F, B, T in layer 1
    ## initialize for B using regression
    M <- batchVec2M(batch_vec) # nxn_b
    fm <- lm(data2use ~ M + 0)
    result$B <- t(fm$coefficients) # pxn_b
    data2use <- fm$residuals
    pb$tick()

    ## initialize for L, F using ppca
    # library(irlba)
    pca_result <- irlba::irlba(data2use, nv = d)
    pb$tick()
    result$L <- pca_result$v %*% diag(pca_result$d) # p x d
    result$F <- pca_result$u # n x d
    varimax_L <- varimax(result$L)
    result$L <- varimax_L$loadings
    result$F <- result$F %*% varimax_L$rotmat
    data2use <- data2use - result$F %*% t(result$L) # n x p
    result$T <- sapply(1:b, function(i) {
        apply(data2use[batch_vec == i, ] * data2use[batch_vec == i, ], 2, sum) / (nrow(object.list[[i]]) - 1)
    })
    result$T <- 1 / result$T # p x b
    pb$tick()

    ###### initialization for c using GMM
    ## initialize c using GMM
    clust_result <- .InitialC(result$F, q = k)
    result$c <- clust_result$c
    result$mu <- clust_result$mu # d x k
    result$Lambda <- solve(clust_result$Lambda) # d x d
    result$omega <- rep(1, n)

    ##### initialize p and gamma
    result$p <- rep(0.5, d)
    result$gamma <- matrix(rbinom(p * d, 1, 0.5), ncol = d) # p x d
    pb$tick()

  # dense0 <- dnorm(result$L, 0, sqrt(lambda0))
  # dense1pMoM <- (result$L^2*dnorm(result$L, 0, sqrt(lambda1)))/lambda1
  # if(verbose) message("Initializing p ...")
  # result$p <- dense1pMoM/(dense1pMoM + dense0)
  # if(verbose) message("Initializing gamma ...")
  # result$gamma <- ifelse(result$p >= 0.5, 1, 0)

  return(result)
}
