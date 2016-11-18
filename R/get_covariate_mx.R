#' Function that returns a matrix of IC coefficients that can be used as covariates.
#'
#' Returning a matrix of IC coefficients that are associated with covariates or have multiple clusters.
#' @param input_list An input ICAobject or PCAobject.
#' @param idx Component indexes that the user wishes to recover.
#' @return Matrix with component coefficients
#' @keywords keywords
#'
#' @export
get_covariate_mx <- function(input_list = NULL, idx = NULL){

    if(is.null(input_list)){
        stop("Input is missing, please specify an ICA or PCA object")
    }

    if(is.null(idx)){
        if(attr(input_list, 'covar_cor') == "yes"){
            cov_idx <- which(sapply(input_list$comp_cov, function(x) length(x) > 1))
        } else {
            cov_idx <- NULL
        }

        if(attr(input_list, 'clustering') == "yes"){
            clust_idx <- which(sapply(input_list$mclust_result, function(x) x$G > 2))
        } else {
            clust_idx <- NULL
        }

        comp_idx <- unique(c(cov_idx, clust_idx))

        if(is.null(comp_idx)){
            stop("No components detected by covariate association and clustering")
        }

    } else {

        comp_idx <- idx

    }

    if(class(input_list) == "ICAobject"){
        cov_mx = t(input_list$A[comp_idx,])
    } else if (class(input_list) == "PCAobject"){
        cov_mx = input_list$x[,comp_idx]
    }

    return(cov_mx)

}
