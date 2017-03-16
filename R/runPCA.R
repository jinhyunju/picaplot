#' Custom PCA function for analyzing gene expression data.
#'
#' Performing PCA on a dataset and create a list object with results.
#'
#' @param pheno_mx Phenotype matrix with diemnsions g x N, where g is the number of genes and N is the number of samples.
#' @param se_obj The input for the function can also be a \code{SummarizedExperiment} or \code{RangedSummarizedExperiment} object created through the
#'        \code{SummarizedExperiment} package on bioconductor.
#' @param assay_idx The assay index to be used to retrieved a single assay from the SummarizedExperiment object.
#' @param scale_pheno Logical value specifying the scaling of row of the pheno_mx. Default is set to FALSE.
#' @return
#' The following entries will be generated in the output list \code{pca_result} after running the example above. \cr
#' \code{rotation} : Matrix of principal component gene weights where each column represents a single component. (standard \code{prcomp()} output)\cr
#' \code{x} : Matrix of the projections of the original data onto principal componets. Each column holds a projection. (standard \code{prcomp()} output)\cr
#' \code{sdev}: The standard deviation (square root of the eigen values) of each principal components. (standard \code{prcomp()} output)\cr
#' \code{percent_var} : The percent variance each principal component is explaining. Calculated based on \code{sdev}.\cr
#' \code{peaks} : Indicating which gene has a gene weight larger than 2 standard deviations of its component gene weights. \cr
#' \code{center} : The mean values for each gene used to center the data. (standard \code{prcomp()} output)\cr
#' \code{scale} : TRUE or FALSE value indicating whether the data was scaled. (standard \code{prcomp()} output)\cr
#' Three attributes are set within the list object. "PCAobject" for \code{class}, "pca" for \code{method} and "no" for \code{covar_cor}.\cr
#'
#'
#' @examples
#'
#' data(expr_data)
#'
#' pca_result <- runPCA(expr_data)
#'
#' @export
#'
runPCA <- function(pheno_mx = NULL,
                    se_obj = NULL,
                    assay_idx = NULL,
                    scale_pheno = FALSE){

    if(is.null(pheno_mx) & is.null(se_obj)){
        stop("Error: Phenotype matrix is missing \n")
    } else if (!is.null(se_obj) & is.null(pheno_mx)){

        if(!is.null(assay_idx)){
            pheno_mx <- SummarizedExperiment::assay(se_obj, assay_idx)
        } else {
            message("Assay number not specified, using the first entry as default")
            pheno_mx <- SummarizedExperiment::assay(se_obj, 1)
        }

    } else if (!is.null(pheno_mx) & is.null(se_obj)){
            message("Using input pheno_mx for expression values")

    } else {
        stop("Both pheno_mx and se_obj are defined, please only supply one
             phenotype matrix as input")
    }

    pheno_nrow <- nrow(pheno_mx)
    pheno_ncol <- ncol(pheno_mx)
    message("Original dimensions of <pheno_mx> = ",
            pheno_nrow , " x ", pheno_ncol, "\n")

    if (pheno_nrow < pheno_ncol){
        message("[Caution] Number of samples exceeding number of measured
                features, please check rows and columns of <pheno_mx> \n")
        message("* If you are from the future and have more samples than
                measured features, disregard the above message and please
                proceed. \n")
    }



    message("Pre Processing Data \n")
    # removing NA values, centering, scaling (optional)
    pheno_mx <- preProcessData(pheno_mx, scale_pheno = scale_pheno)

    message("Running PCA \n")
    pca_result <- stats::prcomp(t(pheno_mx))

    # calculating percent variance each component explains
    pca_result$percent_var <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100


    pca_result$peaks <- apply(pca_result$rotation, 2, peakDetection)

    # label result as pca for report2me() function
    class(pca_result) <- "PCAobject"
    attr(pca_result, 'method') <- "pca"
    attr(pca_result, 'covar_cor') <- "no"
    attr(pca_result, 'clustering') <- "no"

    return(pca_result)
}
