#' Generating ICA or PCA summary reports.
#'
#' Generating a HTML report from a ICA or PCA list object.
#'
#' @param input_list ICA or PCA result list created by either \code{run_ica()} or \code{run_pca()}.
#' @param n_comps Number of components to plot.
#'        Default is set to plot every component.
#' @param prefix Output filename prefix. The output file will be named
#'        "prefix_PICA_summary.html".
#' @param geneinfo_df Dataframe that contains positions of the genes. Column names should be
#'        "pheno_chr" for chromosomes, "pheno_start" for starting position and "pheno_end" for
#'        ending positions.
#' @param output_path Directory path for generating the output HTML file.
#'        default is set to current working directory.
#' @param file.ext File extension to be used for saved plots.
#'        Default is set to png for html reports.
#'        Note that if you use pdf plots for html files they will not show up in your report.
#' @return The output of this function is a HTML report containing component-wise plots and information.
#' @keywords keywords
#'
#' @import ggplot2
#' @import knitr
#' @import rmarkdown
#'
#' @examples
#'
#' data(expr_data)
#'
#' ica_object <- run_ica(expr_data)
#' ica_object <- covar_association_check(ica_object, sample_info)
#' ica_object <- detect_clusters(ica_object)
#'
#' reportgen(ica_object, geneinfo_df = probe_info, output_path = "./", prefix ="example_report")
#'
#'
#'
#' @export
reportgen <- function(input_list = NULL,
                      n_comps = NULL,
                      prefix = NULL,
                      geneinfo_df = NULL,
                      output_path = NULL, file.ext = "png"){

    if(is.null(input_list)){
        stop("Please specify the input to generate a report. \n")
    }

    if(is.null(prefix)){
        stop("Please specify the prefix of the output file \n")
    }

    if(is.null(output_path)){
        message("<output_path> is not specified, using current working directory \n")
        output_path <- getwd()

    }

    if(is.null(n_comps)){
        message("<n_comps> not specified, plotting all components.")
        n_comps <- length(input_list$percent_var)

    }


    phenotype <- NULL # avoid R CMD check NOTES
    method <- attr(input_list, 'method')

    markdown.file <- system.file("templates/PICA_Component_Visualization_Report.Rmd", package="picaplot")

    outFile = paste(output_path,"/",prefix,"_",method,"_summary.html",sep="")

    suppressMessages(rmarkdown::render(markdown.file,output_file = outFile,output_format = "html_document"))


}
