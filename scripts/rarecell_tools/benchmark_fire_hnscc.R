# script include the benchmarking steps for FiRE
# Results of the scripts are included in benchmark_tools.pdf
load("seurat_integrated.rds")
seurat_obj<-seurat_integrated
require('Matrix')
require('plyr')
project ="hnscc"
data<-seurat_obj[["RNA"]]@counts
data<-as.matrix(data)
data<-t(data) #sample * fetures
genes <- c(1:dim(data)[2])
data_mat <- list(mat=data, gene_symbols=genes)
gene_keep <-seq(100,5000 ,by=100)


.normalize_by_umi_2 <-function(x, dataSave, minLibSize, verbose){
  mat  = x$mat
  gene_symbols = x$gene_symbols

  #Filter Cells
  if (minLibSize > 0){

    keepCellIndex <- c()
    for (i in c(1:dim(mat)[1])){
      count = sum(mat[i,] > 0)
      if (count > minLibSize){
        keepCellIndex <- c(keepCellIndex, i)
      }
    }

    mat <- mat[keepCellIndex,]
    if (verbose){
      cat(paste("Dimensions of matrix after cell filtering : ",dim(mat),"\n"))
    }
    write.csv(keepCellIndex, file = paste(dataSave,"keepCellIndexes.csv",sep=""), quote = F,row.names = F)
  }

  #Filter Genes
  cs <- colSums(mat>2)
  x_use_genes <- which(cs > 3)

  x_filt<-mat[,x_use_genes]
  gene_symbols = gene_symbols[x_use_genes]
  if (verbose){
    cat("Dimensions os filtered Matrix:")
    cat(paste(dim(x_filt),"\n"))
  }

  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=gene_symbols)
}


# --------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
.get_variable_gene<-function(m){

  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}

ranger_preprocess<-function(data_mat, ngenes_keep=ngenes_keep, dataSave='./', optionToSave=F, minLibSize=0, verbose=F){

  ngenes_keep = ngenes_keep
  #write.csv(x = data_mat$gene_symbols, file =//;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; "gene_symbols.csv",quote = F,row.names =F)
  #l<-.normalize_by_umi(data_mat)
  l<-.normalize_by_umi_2(data_mat, dataSave, minLibSize, verbose)
  m_n<-l$m

  if (verbose){
    cat("Select variable Genes...\n")
  }

  df<- .get_variable_gene(m_n)
  gc()

  if (verbose){
    cat("Sort Top Genes...\n")
  }

  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]

  if (verbose){
    cat("Cutoff Genes...\n")
  }

  df$used<-df$dispersion_norm >= disp_cut_off

  features = head(order(-df$dispersion_norm),ngenes_keep)
  #system("rm genes",ignore.stderr = T)

  if (optionToSave){
    write.csv(features, file = paste(dataSave,"genes",sep=""), quote = F,row.names = F)
    write.csv(l$use_genes[features], file = paste(dataSave,"genes_used_all",sep=""), quote = F,row.names = F)
  }

  #genes = read.csv(file = "genes")
  #features = genes$x

  #Final Data
  m_n_68K<-m_n[,features]
  m_filt<-Matrix(log2(m_n_68K+1),sparse = T)

  if (verbose){
    cat(paste("Writing Log Normalized whole_matrix, DIM:",dim(m_filt)[1], dim(m_filt)[2]))
  }
  #system("rm whole_matrix",ignore.stderr = T)
  if (optionToSave){
    writeMM(m_filt,file="whole_matrix")
  }

  list(preprocessedData=m_filt, selGenes=features)
}




result<-list()
for (i in gene_keep) {
  result[[length(result) + 1]]<- ranger_preprocess(data_mat, ngenes_keep = i, dataSave = './', optionToSave = FALSE, minLibSize = 0, verbose = TRUE)
  names(result)[length(result)] <- paste("top", i, sep = "_")
}
exprimentID="hpv+hnscc"
save(result,file=paste("results/", exprimentID,"_topgenesmatrix.RData",sep=""))

preprocessedData<-list()
for (i in gene_keep){
  element_name <- paste("top_", i, sep = "")
  if (element_name %in% names(result)) {
    # Extract preprocessedData from the specified element
    preprocessedData [[length(preprocessedData) + 1]]<- as.matrix(result[[element_name]]$preprocessedData)
    names(preprocessedData)[length(preprocessedData)] <- paste("top", i, sep = "_")
  }

}
save(preprocessedData,file=paste("results/", exprimentID,"_preprocessedData.RData",sep=""))
# to look at the dimension of the multiple matrix in the list
lapply(preprocessedData, dim)
rare_cells<-list()
for (i in names(preprocessedData)) {
  fit <- model$fit(preprocessedData[[i]])
  rareness_score <- model$score(preprocessedData[[i]])
  q3 <- quantile(rareness_score, 0.75)
  iqr <- IQR(rareness_score)
  th <- q3 + (1.5 * iqr)
  indIqr <- which(rareness_score >= th)
  #this will be keep tract of index if not initiated for loop will treat as a value and create a list
  #length 5000 which has null information
  rare_cells[[length(rare_cells)+1]] <- rownames(data)[indIqr]
  names(rare_cells)[length(rare_cells)] <- i
}
save(rare_cells,file=paste("results/", exprimentID,"_overlaprarecells.RData",sep=""))
# plot these results in line plot
library(ggplot2)
library(tidyr)
long_data <- pivot_longer(overlapfire_hnscc, cols = starts_with("top."), names_to = "TopGenes", values_to = "OverlapPercentage")
ggplot(long_data, aes(x = OverlapPercentage, y = as.numeric(gsub("top\\.", "", TopGenes)), color = TopGenes)) +
  geom_line() +
  geom_point() +
  labs(title = "Percentage of Overlapping Rare Cells vs. Top Genes",
       x = "Overlap Percentage",
       y = "Top Genes") +
  theme_minimal()


## pca iteration

pca_result<-list()
for (i in seq_along(preprocessedData)) {
  # Call the GapClust function for the current element in 'test'
  tmp<-preprocessedData[[i]]
  cat(dim(tmp))
  pca <- irlba(t(tmp), nv=min(c(50, dim(tmp)-1)))
  pca_result[[length(pca_result) + 1]]<-pca
  cat("Iteration:", i, "\n")
}
save(pca_result,file=paste("pca.RData"))
test<-list()
for (i in seq_along(pca_result)){
  u<-pca_result[[i]]$v
  u<-t(u)
  d<-pca_result[[i]]$d
  variance<-d^2
  prop_variance<-variance/sum(variance)
  iteration_info <- list(u = u, variance = variance,prop_variance=prop_variance)
  # Append the list to the 'test' list
  test[[i]] <- iteration_info

}

gene_keep <-seq(100,5000 ,by=100)
names(test) <- paste("top", gene_keep, sep = "_")
save(test,file="pca.RData")
