#' Function to estimate the number of clusters in IC coefficients or PC projections.
#'
#' Performing ICA on a dataset and create a list object with results.
#' @param input_list An input ICAobject or PCAobject.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @import mclust
#' @export
sample_cluster_check <- function(input_list = NULL){

    if(is.null(input_list)){
        stop("Input is missing, please specify an ICA or PCA object")
    }


    if(class(input_list) == "ICAobject"){
        comp_coeff_mx <- input_list$A
    } else if (class(input_list) == "PCAobject"){
        comp_coeff_mx <- t(input_list$x)
    }

    # Turn off warnings to prevent unneccessary panic
    options(warn=-1)
    mclust_result <- apply(comp_coeff_mx, 1, function(x) Mclust(x))
    # Turn warnings back on to prevent disasters
    options(warn=0)


    input_list$n_clust <- sapply(mclust_result, function(x) x$G)
    input_list$mclust_result <- mclust_result

    attr(input_list, 'clustering') <- "yes"

    return(input_list)
}

