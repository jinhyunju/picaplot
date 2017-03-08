# this is an example script to load expression data published on GEO

if(!require(biomaRt)){
  message("Package biomaRt is not installed")
  message("Installing biomaRt from bioconductor")
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomaRt")  # install biomaRt package
  library("biomaRt")
}

if(!require(GEOquery)){
  message("Package GEOquery is not installed")
  message("Installing GEOquery from bioconductor")
  source("http://bioconductor.org/biocLite.R")
  biocLite("GEOquery") # install the GEOquery package
  library("GEOquery")
}

if(!require(Biobase)){
    message("Package Biobase is not installed")
    message("Installing Biobase from bioconductor")
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase") # install the GEOquery package
    library("Biobase")
}
options('download.file.method'='curl')

gse60028 <- GEOquery::getGEO(GEO = "GSE60028")  # get the GEO data
# (this may take a few minutes depending on your network connection)

geo.eset <- gse60028$GSE60028_series_matrix.txt.gz # extract the expression set

expr_all_probes <- exprs(geo.eset)       # extract the expression matrix
covariate.data <- pData(geo.eset)   # extract sample information
feature.df <- fData(geo.eset)       # extract probe information
feature.df <- feature.df[,c("ID", "ENTREZ_GENE_ID")]

covariate.data <- covariate.data[,c("characteristics_ch1",
                                    "characteristics_ch1.1",
                                    "characteristics_ch1.2",
                                    "characteristics_ch1.3",
                                    "characteristics_ch1.4")]


# load biomart database
# mart.ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
mart.ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="hsapiens_gene_ensembl",
                        host = "www.ensembl.org")

biomart.output <- getBM(attributes=c('entrezgene',
                                     'hgnc_symbol',
                                     'gene_biotype',
                                     'chromosome_name',
                                     'start_position',
                                     'end_position'),
                        filters = 'entrezgene',
                        values = feature.df$ENTREZ_GENE_ID,
                        mart = mart.ensembl)

biomart.output <- subset(biomart.output,
                         chromosome_name %in% c(1:22, "X","Y"))

# match probe information with position information
probe_all_info <- merge(feature.df,
                        biomart.output,
                        by.x = "ENTREZ_GENE_ID",
                        by.y = "entrezgene")


# set the column names
colnames(probe_all_info)[5:7] <- c("pheno_chr",
                                   "pheno_start",
                                   "pheno_end")

colnames(probe_all_info)[2] <- c("phenotype")

probe_all_info <- probe_all_info[,c("phenotype",
                                    "pheno_chr",
                                    "pheno_start",
                                    "pheno_end")]

probe_all_info$pheno_chr <- factor(probe_all_info$pheno_chr,
                                   levels = gtools::mixedsort(unique(probe_all_info$pheno_chr)))

probe_all_info <- probe_all_info[order(probe_all_info$pheno_chr,
                                       decreasing = FALSE),]
# filter expression data down to probes with information
expr_all_probes <- expr_all_probes[as.character(probe_all_info$phenotype),]

# parsing the original covariate.data dataframe into a more natural form
sample_info <- data.frame(apply(covariate.data, 2,
                                 function(x) sapply(strsplit(x, ": "),
                                                    function(a) a[2]))) #extract the data

colnames(sample_info) <- apply(covariate.data, 2,
                               function(x) strsplit(x, ":")[[1]][1])

sample_info$level_of_reaction <- as.character(sample_info$level_of_reaction)

# correct a column which had missing data
sample_info$level_of_reaction[which(is.na(sample_info$level_of_reaction))] <- "0"

# correct a column which had missing data
sample_info$level_of_reaction <- factor(sample_info$level_of_reaction)

# subsetting data to only use 10% of the probes to make it small
set.seed(1984)

probe_info <- probe_all_info[sort(sample(x = c(1:nrow(probe_all_info)),
                                  size = round(nrow(probe_all_info) * 0.1),
                                  replace = FALSE),decreasing = FALSE),]

expr_data <- expr_all_probes[as.character(probe_info$phenotype),]

# remove redundant and temporary objects
rm(gse60028, biomart.output, covariate.data, mart.ensembl, geo.eset, feature.df)
