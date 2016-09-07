#' @import ggplot2
#' @export
plot_comp_chr <- function(input_list = NULL,
                          geneinfo_df = NULL,
                          comp_idx,
                          x.axis = NULL){


    if(is.null(x.axis)){
        x.axis <- chr_axis_creator(geneinfo_df)
    }

    if(is.null(geneinfo_df)){
        stop("Gene position information is missing.")
    }

    if(class(input_list) == "ICAobject"){
        geneinfo_df$idx <- 1:nrow(geneinfo_df)
        geneinfo_df$loading <- input_list$S[as.character(geneinfo_df$phenotype),comp_idx]
        geneinfo_df$peaks <- 1 * ( as.character(geneinfo_df$phenotype) %in% names(input_list$peaks[[comp_idx]]) )

        plot.title <- paste("IC",comp_idx,"_",sum(geneinfo_df$peaks),"peaks", sep = " ")


    } else if (class(input_list) == "PCAobject"){
        geneinfo_df$idx <- 1:nrow(geneinfo_df)
        geneinfo_df$loading <- input_list$rotation[as.character(geneinfo_df$phenotype),comp_idx]
        geneinfo_df$peaks <- 1 * ( as.character(geneinfo_df$phenotype) %in% names(input_list$peaks[[comp_idx]]) )

        plot.title <- paste("PC",comp_idx,"_",sum(geneinfo_df$peaks),"peaks", sep = " ")

    }

    rects1 <- data.frame(xstart = x.axis[1,],
                         xend = x.axis[2,],
                         idx = rep(0),
                         loading = rep(0),
                         fill.col = factor(rep(c(0,1),
                                               length.out = length(x.axis[1,]))))


    p <- ggplot(geneinfo_df, aes_string(x = "idx", y = "loading")) +
        geom_rect(data = rects1, aes_string(xmin = "xstart",
                                            xmax = "xend",
                                            ymin = "-Inf",
                                            ymax = "Inf",
                                            fill = "fill.col"), alpha = 0.35) +
        geom_linerange(aes(ymin = 0, ymax = loading, colour = factor(peaks))) +
        theme_bw() +
        #geom_vline(xintercept=c(0,x.axis[2,]), linetype="dotted") +
        scale_x_continuous("Gene Chromosome Location",breaks=x.axis[3,],labels=unique(geneinfo_df$pheno_chr))+
        scale_y_continuous(expand = c(0, 0)) + labs(title = plot.title) +
        scale_color_manual(values=c("grey", "red")) +
        scale_fill_manual(values=c("skyblue", NA)) +
        ylab("Gene Weights") +
        theme(legend.position="none")

    return(p)
}

#' @import ggplot2
#' @export
plot_comp_no_geneinfo <- function(input_list = NULL,
                                  comp_idx){


    if(class(input_list) == "ICAobject"){

        comp_df <- data.frame("idx" = 1:nrow(input_list$S),
                              "loading" = input_list$S[,comp_idx],
                              "peaks" = 1 * ( as.character(rownames(input_list$S)) %in% names(input_list$peaks[[comp_idx]]) ))

        plot.title <- paste("IC",comp_idx,"_",sum(comp_df$peaks),"peaks", sep = " ")

    } else if (class(input_list) == "PCAobject"){

        comp_df <- data.frame("idx" = 1:nrow(input_list$rotation),
                              "loading" = input_list$rotation[,comp_idx],
                              "peaks" = 1 * ( as.character(rownames(input_list$rotation)) %in% names(input_list$peaks[[comp_idx]]) ))

        plot.title <- paste("PC",comp_idx,"_",sum(comp_df$peaks),"peaks", sep = " ")



    }


    p <- ggplot(comp_df, aes_string(x = "idx", y = "loading")) +
        geom_linerange(aes(ymin = 0, ymax = loading, colour = factor(peaks))) +
        theme_bw() + labs(title = plot.title) +
        scale_color_manual(values=c("grey", "red")) +
        ylab("Gene Weights") +
        theme(legend.position="none")

    return(p)
}

#' @export
comp_coeff_df <- function(input_list = NULL,
                          comp_idx){


    if(class(input_list) == "ICAobject"){

        comp_df <- data.frame("coeff" = input_list$A[comp_idx,],
                              "sample_idx" = 1:ncol(input_list$A))


    } else if (class(input_list) == "PCAobject"){

        comp_df <- data.frame("coeff" = input_list$x[,comp_idx],
                              "sample_idx" = 1:nrow(input_list$x))

    }

    if(attr(input_list, "covar_cor") == "yes"){
        if(length(input_list$comp_cov[[comp_idx]]) > 0){
            associated_covars <- input_list$comp_cov[[comp_idx]]
            max_assoc <- colnames(associated_covars)[which.min(associated_covars)]
            comp_df$covar <- input_list$covars[,max_assoc]
        }

    }

    return(comp_df)

}

#' @import ggplot2
#' @export
plot_single_component <- function(input_list = NULL,
                                  comp_idx,
                                  geneinfo_df = NULL,
                                  plot_here = TRUE){


    component_plots <- list()


    if(!is.null(geneinfo_df)){
        component_plots[[1]] <- plot_comp_chr(input_list, comp_idx = comp_idx, geneinfo_df = geneinfo_df)
    } else {
        component_plots[[1]] <- plot_comp_no_geneinfo(input_list, comp_idx = comp_idx)
    }

    coeff_plot_df <- comp_coeff_df(input_list, comp_idx)

    if(is.null(coeff_plot_df$covar)){
        component_plots[[2]] <- ggplot(coeff_plot_df, aes(x = sample_idx, y = coeff)) +
            geom_point() + theme_bw()+ theme(legend.position = "none")

        component_plots[[3]] <- ggplot(coeff_plot_df, aes(x = coeff)) +
            geom_histogram() + coord_flip() + theme_bw()

    } else if (!is.null(coeff_plot_df$covar)){

        component_plots[[2]] <- ggplot(coeff_plot_df, aes(x = sample_idx, y = coeff, col = covar)) +
            geom_point() + theme_bw()+ theme(legend.position = "none")

        component_plots[[3]] <- ggplot(coeff_plot_df, aes(x = coeff, fill = covar)) +
            geom_histogram() + coord_flip() + theme_bw()

    }

    if(plot_here){

        suppressMessages(multiplot(plotlist = component_plots,
                                   layout = matrix(c(1,1,1,1,2,2,3,3),
                                                   ncol =4 , byrow = TRUE)))

    }
    return(component_plots)


}

#' Function that generates the x axis ticks for chromosome positions
#'
#' The function takes in a gene position info data frame with the following columns:
#' "phenotypes" holding the names of the genes used in the analysis (rownames of phenotype.mx)
#' "pheno_chr" holding the chromosome each gene belongs to
#' "pheno_start" the starting position of a given gene
#' "pheno_end" the end position of a given gene
#' "idx" a non-overlapping index that shows orders the genes based on their position, which goes from
#' 1 to number of genes.
#'
#' @param geneinfo.df dataframe that holds the position of the genes.
#'
#' @export
chr_axis_creator <- function(geneinfo.df){
    pheno_chr <- NULL # to get rid of R CMD check NOTES
    chromosomes <- unique(geneinfo.df$pheno_chr)
    x.axis <- matrix(0, nrow = 3, ncol = length(chromosomes))
    geneinfo.df$idx <- 1:nrow(geneinfo.df)
    for(k in 1:length(chromosomes)){
        j <- as.character(chromosomes[k])

        if(is.na(j)){
            x.axis[1,k] <- max(x.axis[1,]) + 1
            x.axis[2,k] <- max(geneinfo.df[,"idx"])
        } else {
            x.axis[1,k] <- min(subset(geneinfo.df, pheno_chr == j)[,"idx"])
            x.axis[2,k] <- max(subset(geneinfo.df, pheno_chr == j)[,"idx"])
        }

    }

    x.axis[3,] <- (x.axis[1,] + x.axis[2,]) / 2
    return (x.axis)
}

#' Plotting multiple ggplots
#'
#' Equivalent function for ggplots to setting par(mfrow=c(x,y)) for
#' basic R plots.
#' Adapted from : http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'
#' @param plotlist list of plots
#' @param cols number of columns
#' @param layout layout matrix for custom layout
#' @return single plot object
#'
#' @keywords keywords
#'
#'
#' @export
multiplot <- function(plotlist=NULL, layout=NULL, cols=1) {

    # Make a list from the ... arguments and plotlist
    #plots <- c(list(...), plotlist)
    plots <- plotlist
    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
        }
    }
}
