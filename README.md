The `picaplot` package provides functions that perform principal or independent component analysis and easy visualization of the results. 

## Quick Start

### Important Note

The report generating feature of this package requires a program called ```pandoc``` version 1.12.3 or higher to be installed. This is not a problem when you are using the most recent version of Rstudio (0.98.1102 at the moment - 2 March 2015), since it comes with the required ```pandoc``` functionality. However, if you are running ```picaplot``` from good-old-fashioned-R (whether you are running it on a cluster or on your local machine) you might run into an error message that pandoc is not installed. In such a case, please review the following link to install the correct version of pandoc on your machine. (https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md)

### 1.Loading and Preparing an Example Dataset

If you already know what you want to use the package for, follow this simple example to get started!

#### 1) Installing the package through the function ```install_github``` from the package ```devtools```.

- `picaplot` has a few dependencies and you need to install them manually before you install `picaplot`. Simply run the code below to install the dependencies.

```r

picaplot_dependencies <- c("ggplot2", "knitr", "rmarkdown", "devtools")

install.packages(picaplot_dependencies)

``` 

- To install `picaplot` from github:

```r

devtools::install_github("jinhyunju/picaplot") #installing picaplot

library(picaplot)

```

#### 2) Loading an example dataset

- Here we are going to use a public dataset that is available on the Gene Expression Omnibus (GEO).

- You can use your own dataset, it just needs to be a matrix which has the dimension of (gene x samples). A dataframe with covariate information is optional with dimensions (samples x covariates). 

- To generate the example dataset, all you have to do is source the script included in the package. 

- If you are interested in the details of the script you can check the path to the script by printing out the value that is saved in the ```example.data.script``` object and open it up in any text editor. 

- Please be aware that the script will install two packages ```GEOquery``` and ```biomaRt``` if you don't already have it on your machine. 

```r
example.data.script <- system.file("templates/create_example_data.R", package="picaplot")

source(example.data.script)
```

- This will take a few minutes depending on your internet connection, since it is downloading data from GEO and biomaRt. 

- If everything ran correctly, it will generate 3 objects, ```expr.data```, ```sample.info```, and ```probe.info```.

- ```expr.data``` = gene expression measurements for 26391 probes and 47 samples. 

- ```sample.info``` = 5 covariates for each sample 

- ```probe.info``` = positional information for 26391 probes
  
  
- One thing that you have to watch out for is that the rownames of ```sample.info``` have to match the column names of the ```expr.data```. They don't necessarily have to be in the same order but they should have names that overlap. 


### 2. Core Functionality of the package

#### 1) Running PCA / ICA 

The functions for running PCA / ICA on an expression matrix are `run_pca()` and `run_ica()` respectively.  

```r 
# run PCA 
pca_object <- run_pca(expr.data)

# run ICA
ica_object <- run_ica(expr.data)

```

This generates a PCAobject / ICAobject with the outputs saved in the format of a list.

1-1) `run_pca()` outputs

The following entries will be generated in the output list `pca_object` after running the example above. 

  * `rotation` : Matrix of principal component gene weights where each column represents a single component. (standard `prcomp()` output)
  
  * `x` : Matrix of the projections of the original data onto principal componets. Each column holds a projection. (standard `prcomp()` output)
  
  * `sdev` : The standard deviation (square root of the eigen values) of each principal components. (standard `prcomp()` output)
  
  * `percent_var` : The percent variance each principal component is explaining. Calculated based on `sdev`
  
  * `peaks` : Indicating which gene has a gene weight larger than 2 standard deviations of its component gene weights. 
  
  * `center` : The mean values for each gene used to center the data. (standard `prcomp()` output)
  
  * `scale` : TRUE or FALSE value indicating whether the data was scaled. (standard `prcomp()` output)
  
  * Three attributes are set within the list object. "PCAobject" for `class`, "pca" for `method` and "no" for `covar_cor`.

1-2) `run_ica()` outputs

The following entries will be generated in the output list `ica_object` after running the example above. 

  * `A` : The IC coefficient matrix, with each row representing coefficients for the corresponding independent component. (standard `fastICA()` output)
  
  * `S` : Matrix of gene weights for each independent component. Each column holds a single component. (standard `fastICA()` output)

  * `percent_var` : The percent variance each independent component is explaining.
  
  * `peaks` : Indicating which gene has a gene weight larger than 2 standard deviations of its component gene weights. 
  
  * `order` : The order of independent components based on the variance that they explain. 
  
  * `X`, `K`, `W` : Standard outputs of `fastICA()`. `X` is the pre-processed data matrix, `K` is the pre-whitening matrix projecting the data onto the first n principal components, and `W` is the estimated unmixing matrix. 
  
  * Three attributes are set within the list object. "ICAobject" for `class`, "ica" for `method` and "no" for `covar_cor`.

2) Testing Associations Between Covariates and Components

- As an optional step, you can check whether any covariates are associated with any of the components. 
`covar_association_check()` can be applied to both PCA and ICAobjects.

- In this example, we are going to test the associations between the covariats in `sample.info` and each PC and IC. 

```r 

pca_object <- covar_association_check(pca_object, 
                                      covars = sample.info)

ica_object <- covar_association_check(ica_object, 
                                      covars = sample.info)

```

- This will add the following entries to the list.

  * `covar_pvals` : A matrix of p-values with the dimension of number of components x tested covariates. 
  
  * `comp_cov` : A list with length equal to the number of components that shows in each entry which covariates have a p-value lower than the set threshold. 
  
  * `covars` : A copy of the supplied `sample.info` for plotting. 
  
  * `covar_threshold` : The threshold for calling a covariate association significant. The default is set to 0.05 divided by the number of tests (= `length(covar_pvals)`).

3) Plotting Individual Components

- Individual components can be inspected by using the `plot_component()` function. 
This will generate 3 plots showing the gene loading on the component of interest and the component coefficients.

- You can specify which component to inspect by setting the component index in the option `comp_idx`.

```r

pca1_plot = plot_component(pca_object, 
                           comp_idx = 1)

ica1_plot = plot_component(ica_object,
                           comp_idx = 1)

```

- If you have information regarding the chromosome and position of each gene you can supply it to the function to color the gene loading plot by chromosome. 

- This `geneinfo_df` dataframe needs the following columns: `phenotype` which has an entry for all rownames of the `expr.data`, 
`pheno_chr` showing which chromosome the corresponding gene is on, `pheno_start` for the starting base position of the given phenotype, 
`pheno_end` for the end base position of the phenotype. 

- Note that for ideal results you want to have the `geneinfo_df` sorted by chromosome and position. 
`plot_component` uses the order of `phenotype` given in the dataframe to generate the plot, so the user has control over how to plot it by changing the order of the `geneinfo_df` dataframe.  

```r

pca1_color = plot_component(pca_object, 
                            comp_idx = 1, 
                            geneinfo_df = probe.info)

ica1_color = plot_component(ica_object, 
                            comp_idx = 1, 
                            geneinfo_df = probe.info)

```

4) Generate a Report for All Components

- To create an HTML report showing all components with more detail use the `reportgen()` function.

- You can control the number of components to be plotted by setting `n_comps`, if `n_comps` is not specified it will plot every component. The default order is set by the amount of variance each component explains

- If option `output_path` is not set, it will generate the report and plots in the current working directory.

```r

reportgen(pca_object, n_comps = 10, prefix = "PCAreport")

reportgen(ica_object, n_comps = 10, prefix = "ICAreport")


```
