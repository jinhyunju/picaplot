# Pre-processing of the data
#
# The function takes in a g x N matrix, where g is the number of features and N is the number of samples,
# and removes any 0 variance entries followed by mean centering.
# If scaling option is true it will scale all rows to have unit standard deviation.
#
# @param input_matrix Input matrix with a dimension of g X N  \code{input_matrix}
# @param scale_pheno If TRUE features are scaled to have unit standard deviation.
# @return Pre-processed data matrix
# @keywords keywords
#
# @examples
# raw.data <- rbind(rnorm(100,2,3),rnorm(100,10,2))
# centered <- pre_process_data(raw.data, scale_pheno = FALSE)
# apply(centered, 1, mean)
# apply(centered, 1, sd)
# centered.scaled <- pre_process_data(raw.data, scale_pheno = TRUE)
# apply(centered.scaled, 1, mean)
# apply(centered.scaled, 1, sd)
# @export
pre_process_data <- function(input_matrix, scale_pheno){
    var.vec <- apply(input_matrix, 1, stats::var)
    var0.count <- sum(var.vec == 0)

    # Remove 0 variance  if existing
    if(var0.count != 0){
        message("- Removing ",var0.count," 0 variance features \n")
        input_matrix <- input_matrix[-c(which(var.vec == 0)),]
    } else {
        message("- No 0 variance features in dataset \n")
    }


    if(scale_pheno == TRUE){
        message("- Centering and Scaling Input Matrix \n")
        # scale Input matrix to have 0 mean and 1 sd
        input_matrix <- t(apply(input_matrix, 1, function(x) (x - mean(x))/stats::sd(x)))
    } else {
        message("- Centering Input Matrix without Scaling \n")
        input_matrix <- t(apply(input_matrix, 1, function(x) (x - mean(x))))
    }

    message("- Dimension of Input Matrix = [",
            nrow(input_matrix),",",ncol(input_matrix),"]\n")
    return(input_matrix)
}

