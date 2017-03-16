#' Custom ICA function for analyzing gene expression data.
#'
#' Performing ICA on a dataset and create a list object with results.
#'
#' @param pheno_mx Phenotype matrix with diemnsions g x N, where g is the number of genes and N is the number of samples.
#' @param se_obj The input for the function can also be a \code{SummarizedExperiment} or \code{RangedSummarizedExperiment} object created through the
#'        \code{SummarizedExperiment} package on bioconductor.
#' @param assay_idx The assay index to be used to retrieved a single assay from the SummarizedExperiment object.
#' @param k_est Number of components to be estimated or method to estimate it.
#' @param var_cutoff Percent variance threshold to use when <k_est> is not supplied.
#' @param n_runs Number of times to run ICA. If this is set to a number larger than 1,
#' only ICs that replicate between runs are going to be returned
#' @param n_cores Number of cores to use when n_runs is larger than 1.
#' @param scale_pheno If set to TRUE the pre-processing step of the data will include a scaling step.
#' This will divide all phenotypes by their standard deviation.
#' @param h_clust_cutoff is the cutoff value used in hierarchical clustering. Default is set to 0.3.
#' @param max_iter Maximum iterations for estimating k for each run. Default value is set to 10.
#' @param similarity_measure How to measure the similarity between ICs.
#' If set to "peaks" only gene weights that are greater than 1 sd are used to calculate similarity.
#' @return The following entries will be generated in the output list \code{ica_object} after running the example above. \cr
#' \code{A} : The IC coefficient matrix, with each row representing coefficients for the corresponding independent component. (standard \code{fastICA()} output) \cr
#' \code{S} : Matrix of gene weights for each independent component. Each column holds a single component. (standard \code{fastICA()} output)\cr
#' \code{percent_var} : The percent variance each independent component is explaining.\cr
#' \code{peaks} : Indicating which gene has a gene weight larger than 2 standard deviations of its component gene weights.\cr
#' \code{order} : The order of independent components based on the variance that they explain.\cr
#' \code{X, K, W} : Standard outputs of \code{fastICA()}. \code{X} is the pre-processed data matrix, \code{K} is the pre-whitening matrix projecting the data onto the first n principal components, and `W` is the estimated unmixing matrix.\cr
#' Three attributes are set within the list object. "ICAobject" for \code{class}, "ica" for \code{method} and "no" for \code{covar_cor}.
#'
#' @examples
#'
#' data(expr_data)
#'
#' ica_result <- runICA(expr_data)
#'
#'
#' @export
runICA <- function(pheno_mx = NULL,
                    se_obj = NULL,
                    assay_idx = NULL,
                    k_est = NULL,
                    var_cutoff = 99,
                    n_runs = 1,
                    n_cores = NULL,
                    scale_pheno = FALSE,
                    h_clust_cutoff = 0.3,
                    max_iter = 10,
                    similarity_measure = "peaks"){

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


    if(is.null(colnames(pheno_mx))){
        message("<pheno_mx> is missing column names, setting to default (sample#). \n")
        colnames(pheno_mx) <- paste("sample",c(1:pheno_ncol), sep = "_")
    }

    if(is.null(rownames(pheno_mx))){
        message("<pheno_mx> is missing row names, set to default. \n")
        rownames(pheno_mx) <- paste("feature",c(1:pheno_nrow), sep = "_")
    }

    if(is.null(n_cores)){
        message("<n_cores> not specified.
                Using single core for calculations. \n")
        n_cores = 1
    }

    message("------ Pre-processing Data ------- \n")
    # removing 0 variance genes and scaling and centering the phenotype matirx
    pheno_mx <- preProcessData(pheno_mx, scale_pheno = scale_pheno)

    if(is.null(k_est)){
        message("Missing input for <k_est>, using the number of principal
                components explaining 99% of total variance \n")
        # running PCA on data
        pca_pheno <- stats::prcomp(t(pheno_mx))
        percent <- (cumsum(pca_pheno$sdev^2) /sum(pca_pheno$sdev^2)) * 100
        k_est <- which(percent > var_cutoff)[1]
        message(k_est," components needed to explain more than ",
                var_cutoff,"% of the variance \n")
        message("<k_est> set to ", k_est, "\n")

        if(k_est == 1){
            stop("1 component explains more than", var_cutoff,
                 "% of the variance, check your data or
                 set <k_est> to a number bigger than 1 \n")
        }
    }

    ica_list <- list()
    ica_result <- list()


    message("------ Running ICA ------- \n")
    message("* This may take some time depending on the size of your dataset \n")
    if(n_runs > 1){
        message("<n_runs> set to ", n_runs,
                " / initiating multiple ICA runs \n")
        message("- Generating ", n_runs,
                " replicates using ", n_cores, " core(s) \n")
        message("- Estimated components = ", k_est, "\n")

        ica_list <- parallel::mclapply(1:n_runs,
                                       function(x)
                                           fastICAgeneExpr(pheno_mx, k_est,
                                                             fun = "logcosh",
                                                             alpha = 1,
                                                             scale_pheno = FALSE,
                                                             maxit = 500,
                                                             tol = 0.0001,
                                                             verbose = FALSE),
                                       mc.cores = n_cores)

        for(i in 1:length(ica_list)){

            ica_list[[i]]$peak_mx <- apply(ica_list[[i]]$S,
                                           2,
                                           function(x)
                                               1*(abs(x) > 2*stats::sd(x)))

        }

        # After running ICA sevaral times
        # combined.A <- do.call(rbind, lapply(ica_list, function(x) x$A))
        # combine all components into a single matrix
        combined_S <- do.call(cbind, lapply(ica_list, function(x) x$S))

        # combine all peak position matrices as well
        peak_matrix <- do.call(cbind, lapply(ica_list, function(x) x$peak_mx))

        if(similarity_measure == "peaks"){
            if(0 %in% apply(peak_matrix, 2, sum)){
                cor.mx <- stats::cor(combined_S)
            } else{
                # do element-wise multiplication to only save values for peaks
                peak.component <- peak_matrix * combined_S
                # calculate correlation between components
                # (only with their peak values)
                cor.mx <- stats::cor(peak.component)
            }
        } else {
            # if similarity measure is not based on "peaks"
            # use all elements for calculating the similarity between estimates
            cor.mx <- cor(combined_S)
        }

        # clustering of the components
        # create dissimilarity matrix
        dissimilarity <- 1 - abs(cor.mx)

        # convert into distance matrix format for clustering
        cor.dist <- stats::as.dist(dissimilarity)

        # run hierarchical clustering
        h_clust <- stats::hclust(cor.dist)

        # cut tree at h_clust_cutoff (absolute correlation > 1 - h_clust_cutoff)
        groups <- stats::cutree(h_clust, h=h_clust_cutoff)

        # count member components for each group
        group_table <- table(groups)

        # get groups with more than n_runs * 0.9 members
        multi_component_group <- which(group_table >= (n_runs * 0.9))

        k_update <- length(multi_component_group)
        message("- ", k_update," Replicating Components Estimated \n")

        if(k_update < 1){
            stop('None of the ICs replicated.
                 You could try to the increase sample size or <n_runs>. \n')
        }

        # Averaging the components estimated in multiple runs
        Avg_S <- matrix(0,nrow = dim(combined_S)[1],ncol = k_update)

        # for each group calculate the average component
        for(i in 1:length(multi_component_group)){

            # get component indexes for groups with multiple components
            group_members <- which(groups %in% multi_component_group[i])

            # subset matrix for those components
            sub_group_S <- combined_S[,group_members]

            # calculate correlation between them
            group_cor_mx <- cor(sub_group_S)

            # in order to average the components
            # signs need to be matched (positive vs negative)
            match_sign <- ifelse(group_cor_mx[,1] < 0, -1,1 )

            # calculate the mean component
            avg_component <- (sub_group_S %*% match_sign) / length(match_sign)

            Avg_S[,i] <- avg_component

        }

        hclust_list <- list()
        hclust_list$plot <- h_clust
        hclust_list$cutoff <- h_clust_cutoff
        hclust_list$dist.mx <- dissimilarity
        ica_result$hclust <- hclust_list
        ica_result$S <- Avg_S
        ica_result$A <- solve(t(Avg_S) %*% Avg_S) %*% t(Avg_S) %*% pheno_mx

        rm(Avg_S)

    } else if (n_runs == 1){
        message("<n_runs> set to ", n_runs, " / initiating single ICA run \n")
        message("- Estimated components = ", k_est, "\n")

        ica_result <- fastICAgeneExpr(pheno_mx, k_est,
                                        fun = "logcosh",
                                        alpha = 1, scale_pheno = FALSE,
                                        maxit=500, tol = 0.0001, verbose = FALSE)
        k_update <- k_est
    }

    # Setting appropriate names for signals and mixing matrix
    rownames(ica_result$S) <- rownames(pheno_mx)
    colnames(ica_result$S) <- paste("IC",c(1:dim(ica_result$S)[2]),sep="")
    colnames(ica_result$A) <- colnames(pheno_mx)
    rownames(ica_result$A) <- paste("IC",c(1:dim(ica_result$A)[1]),sep="")

    # Attaching the sample info dataframe to the ica list


    message("------ Post-ICA Processing ------- \n")

    message("- Labeling peak genes in each IC \n")
    # peaks are defined as gene contributions that are larger than 2 standard deviations

    ica_result$peaks <- apply(ica_result$S, 2, peakDetection)

    message("- Calculating variance explained by each IC \n")
    # get the total variance by the sums of squares of the scaled pheno_mx
    total_var <- sum(pheno_mx^2)

    # applying IC component-wise variance calculations
    var_IC <- sapply(1:dim(ica_result$A)[1],
                     function (x) ICvarianceCalc(ica_result$S[,x],
                                                   ica_result$A[x,]))

    # % variance explained by each IC
    percent_var <- (var_IC / total_var) * 100

    message("- Sanity Check : Total % of variance explained by ",
            k_update," ICs = ", sum(percent_var), "\n")

    message("- Creating index based on Variance explained \n")
    # ordering the ICs based on the amount of variance they explain
    ica_result$order <- order(percent_var,decreasing = TRUE)
    ica_result$percent_var <- percent_var


    class(ica_result) <- "ICAobject"
    attr(ica_result, 'method') <- "ica"
    attr(ica_result, 'covar_cor') <- "no"
    attr(ica_result, 'clustering') <- "no"

    message("------ Process completed without any interuptions -------  \n")

    return(ica_result)
}


# Labeling peaks for each independent component
#
# The function takes in the A matrix from the ica.result object
# with a dataframe of measured covariates to test the association
# between them.
#
# @param ica.input.s A single column of the S matrix with dimensions 1 x g \code{s}
# @return A vector that contains the gene weights of the signficant peaks for a single independent component sorted
# in the decreasing order of absolute magnitude.
#
peakDetection <- function(ica.input.s){
    peaks.idx <- which(abs(ica.input.s) > (2 * stats::sd(ica.input.s)))
    # get the peak indexes
    peaks <- ica.input.s[names(sort(abs(ica.input.s[peaks.idx]),
                        decreasing = TRUE))]
    # sort them in decreasing order of absolute magnitude
    return(peaks)
}

# Calculating the variation explained by each component
#
# The function takes in a single row of the ICA result A matrix and a
# single column of the ICA result S matrix and calculates the sums of squares.
#
# @param s A single column of the S matrix with dimensions 1 x g  \code{s}
# @param a A single row of the A matrix with dimensions N x 1 \code{a}
# @return A single number of how much variantion is explained by a single IC.
#
ICvarianceCalc <- function(s, a){
    var.IC <- sum(  (as.matrix(s) %*% t(as.matrix(a) ) )^2)
    return(var.IC)
}
