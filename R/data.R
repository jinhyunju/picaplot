#' A subset of gene expression measurements with 2642 genes and 47 samples.
#'
#'
#' A dataset containing the gene expression measurements obtained from the
#' Gene Expression Omnibus (GEO) accession number GSE60028.
#'
#' The subset was created by randomly sampling 10% of the original 26419 genes.
#'
#' The original dataset is described in the following paper:
#' Dhingra N, Shemer A, Correa da Rosa J, Rozenblit M et al.
#' Molecular profiling of contact dermatitis skin identifies allergen-dependent differences in immune response.
#' J Allergy Clin Immunol 2014 Aug;134(2):362-72. PMID: 24768652
#'
#' @format A matrix with 2642 rows and 47 columns:
#' \describe{
#'   \item{rows}{Each row is a gene expression measurement}
#'   \item{cols}{Each column is a single sample}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60028}
"expr_data"

#' Sample covariate information for 47 samples in \code{expr_data}
#'
#' A dataset containing the sample covariates for \code{expr_data}.
#'
#' The original dataset is described in the following paper:
#' Dhingra N, Shemer A, Correa da Rosa J, Rozenblit M et al.
#' Molecular profiling of contact dermatitis skin identifies allergen-dependent differences in immune response.
#' J Allergy Clin Immunol 2014 Aug;134(2):362-72. PMID: 24768652
#'
#' @format A matrix with 47 rows and 4 columns:
#' \describe{
#'   \item{gender}{Gender of the sample}
#'   \item{allergen}{Allergen that the sample was exposed to}
#'   \item{group}{Group of allergen}
#'   \item{level_of_reaction}{Level of reaction after allergen exposure}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60028}
"sample_info"

#' Positional information for 2642 probes used in \code{expr_data}
#'
#'
#' The original dataset is described in the following paper:
#' Dhingra N, Shemer A, Correa da Rosa J, Rozenblit M et al.
#' Molecular profiling of contact dermatitis skin identifies allergen-dependent differences in immune response.
#' J Allergy Clin Immunol 2014 Aug;134(2):362-72. PMID: 24768652
#'
#' @format A matrix with 2642 rows and 4 columns:
#' \describe{
#'   \item{phenotype}{The probe id used in \code{expr_data}}
#'   \item{pheno_chr}{Chromosome of measured probe.}
#'   \item{pheno_start}{Starting position of measurement probe}
#'   \item{pheno_end}{End position of measurement probe}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60028}
"probe_info"
