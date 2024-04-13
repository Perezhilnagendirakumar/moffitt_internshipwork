# script for performing correlation analysis of results from benchmark_fire_hnscc.R and benchmark_gapclust_hnscc.R 
# Results of the scripts are present in correlation_analysis1to2000_fire.xlsx and correlation_analysis1to2000_gapclust.xlsx 

#load("pca.RData")
# test has the information of PCA iteration result 
test<-test[1:20]
for (i in 1:(length(test)-1)) { # 2:(min(3, length(test))) to iterate only with particular range 
  # Extract the matrices from the list
  matrix1 <- test[[i]]$u
  v1<- test[[i]]$variance
  prop_variance_1<-test[[i]]$prop_variance
  matrix2 <- test[[i + 1]]$u
 v2<-test[[i + 1]]$variance
  prop_variance_2<-test[[i+1]]$prop_variance
  result_matrix_u <- matrix(NA, nrow = nrow(matrix1), ncol = nrow(matrix2))
  for (j in 1:nrow(matrix1)) {
    for (k in 1:nrow(matrix2)) {
      # Extract the rows from the two matrices
      row1 <- matrix1[j, ]
      row2 <- matrix2[k, ]
      
      # Perform the correlation test
      correlation_result <- cor.test(row1, row2, method = "pearson")
      
      # Extract and store the correlation coefficient
      correlation_coefficient <- correlation_result$estimate
      result_matrix_u[j, k] <- correlation_coefficient
    }
  }
  
  correlation_list <- list(v1 = v1, v2 = v2,prop_variance_1=prop_variance_1,prop_variance_2=prop_variance_2, matrix = result_matrix_u)
  # Update the list with named entries for u and variance
  correlation_matrices[[paste(i, i+1, sep = "_vs_")]] <-correlation_list
}
correlation_matrices<-list()
# to change negative correlations to positive 
correlation_matrices <- lapply(correlation_matrices, function(element) {
  element$matrix <- abs(element$matrix)
  return(element)
})
# Print or use the resulting list of correlation matrices (correlation_matrices)
str(correlation_matrices)
# correlation matrices --> further updated to dataframe and colnames and rownames added in the name of test 
save(test,file=paste("results/overlap_fire","correlation_matrix_1to 2000.RData"))

# change the matrix to dataframe 
test<-lapply(correlation_matrices, function(matrix) data.frame(matrix))

new_colnames <- paste("pc", 1:50, sep = "")
new_rownames <- paste("pc", 1:50, sep = "")

# Use lapply to update each dataframe in the list
test <- lapply(test, function(df) {
  # Set new column names
  #colnames(df) <- new_colnames
  colnames(df)[5:length(df)] <- new_colnames
  # Set new row names
  rownames(df) <- new_rownames
  return(df)
})

# Create a named list with dataframes and sheet names
dataset_names <- setNames(test, names(test))

# Export each dataframe to separate sheets in the same Excel file
openxlsx::write.xlsx(dataset_names, file = paste("results/correlation_analysis1to2000_fire.xlsx"))
openxlsx::write.xlsx(dataset_names, file ='correlation_analysis1to2000_fire.xlsx')
### compare one vs another 
result_matrix <- matrix(NA, nrow = nrow(test), ncol = nrow(test_1))
# for top 100 vs top 200 
for (i in 1:nrow(test)) {
  for (j in 1:nrow(test_1)) {
    # Extract the rows from the two matrices
    row1 <- test[i, ]
    row2 <- test_1[j, ]
    
    # Perform the correlation test
    correlation_result <- cor.test(row1, row2, method = "pearson")
    
    # Extract and store the correlation coefficient
    correlation_coefficient <- correlation_result$estimate
    result_matrix[i, j] <- correlation_coefficient
  }
}
result_matrix<-data.frame(result_matrix)

write.xlsx(result_matrix,"correlation_analysis.xlsx")
