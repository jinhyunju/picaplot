#' Custom ICA function for analyzing gene expression data.
#'
#' Performing ICA on a dataset and create a list object with results.
#' @param input_list An input ICAobject or PCAobject.
#' @param covars A matrix of covariates with each row for a sample and each column for a covaraite
#' @param cor_threshold Threshold that is going to be used for identifying significant associations.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @export
covar_association_check <- function(input_list = NULL,
                                    covars = NULL,
                                    cor_threshold = 0.05){

    if(is.null(input_list)){
        stop("Input is missing, please specify an ICA or PCA object")
    }

    if(is.null(covars)){
        stop("Covariate matrix is missing, please specify a covariate matrix")
    }

    if(class(input_list) == "ICAobject"){
        comp_coeff_mx <- input_list$A
    } else if (class(input_list) == "PCAobject"){
        comp_coeff_mx <- t(input_list$x)
    }

    message("Checking dimensions and column/row names \n")

    if(nrow(covars) == ncol(comp_coeff_mx)){
        message("[Good] Number of samples in <covars> and <input_list> match \n")
    } else {
        stop("Error: Number of samples in <covars> and <input_list> don't match.")
    }

    matching.names <- sum(rownames(covars) %in% colnames(comp_coeff_mx))

    if(matching.names == nrow(covars)){
        message("[Good] All samples in <input_list> are accounted for in <covars> \n")
    } else {
        stop("Sample names do not match between <input_list> and <covars> : Check sample names in <covars> and <input_list>")
    }


    n_unique <- apply(covars, 2, function (x) length(table(x, useNA = 'no')))
    use_info <- n_unique > 1 & n_unique < nrow(covars)

    warning(sum(!use_info), " columns excluded from <covars>, due to uniqueness issues.")

    covars <- covars[, use.info]
    
    message("- Checking associations between components and covariates \n")
    # Anova analysis for covariates vs ICA weights (A matrix)
    covar_pvals <- ic_covariate_association_test(comp_coeff_mx, covars)

    if(class(input_list) == "ICAobject"){
        rownames(covar_pvals) <- paste0("IC",1:nrow(covar_pvals))
    } else if (class(input_list) == "PCAobject"){
        rownames(covar_pvals) <- paste0("PC",1:nrow(covar_pvals))
    }

    bon_sig <- cor_threshold / (nrow(covar_pvals) * ncol(covar_pvals))
    comp_cov <- list()

    for(i in 1:nrow(covar_pvals)){

        comp_name <- rownames(covar_pvals)[i]
        comp_cov[[comp_name]] <- covar_pvals[i,covar_pvals[i,] < bon_sig, drop = FALSE]

    }

    input_list$covar_pvals <- covar_pvals
    input_list$comp_cov <- comp_cov
    input_list$covars <- covars

    attr(input_list, 'covar_cor') <- "yes"

    return(input_list)
}

#' Testing the association between independent component coefficients and known covariates.
#'
#' The function takes in the A matrix from the ica.result object
#' with a dataframe of measured covariates to test the association
#' between them.
#'
#' @param input.A A matrix of an ica_result list
#' @param info.input A dataframe that holds the measured covariates for each sample
#' @return output A matrix holding the p-values for each indepenedent component and covariate pair.
#' @keywords keywords
#'
ic_covariate_association_test <- function(input.A, info.input){

    n.components <- dim(input.A)[1]
    covar.names <- colnames(info.input)
    n.covariates <- length(covar.names)

    pval.mx <- matrix(NA, nrow = n.components, ncol = n.covariates)
    colnames(pval.mx) <- covar.names

    covar.df <- data.frame(info.input[colnames(input.A),covar.names])

    for( i in 1: dim(input.A)[1]){
        for (j in 1:length(covar.names)){
            analysis.df <- data.frame("IC"=input.A[i,],"Covar" = covar.df[,j])
            anova.fit <- aov(IC ~ Covar, data=analysis.df)
            model.fstat <- summary.lm(anova.fit)$fstatistic
            p.val <- pf(model.fstat[1],model.fstat[2],model.fstat[3], lower.tail = FALSE, log.p = FALSE)
            pval.mx[i,j] <- p.val
        }
    }

    return(pval.mx)
}
