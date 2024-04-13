####preprocessing steps to retrieve the top 100 to 5000 gene * cells matrix
#load("seurat_integrated.rds")
#hnscc<- CreateSeuratObject(counts = seurat_integrated)
#hnscc <- FindVariableFeatures(object = hnscc, selection.method='vst', nfeatures=dim(data)[1], verbose = F)
#vst<- as.matrix(HVFInfo(pbmc))

# Results of the scripts are present in "Benchmarktools.pdf"

#vst<-data.frame(vst[, "variance.standardized"])
#colnames(vst)<-"variance.standardized"

#vst <- vst[order(-vst$variance.standardized), , drop = FALSE]

#save(vst,file=paste("results/", exprimentID="hnscc","_vst.RData",sep=""))

gene_keep <-seq(100,5000 ,by=100)
result<-list()
for (i in gene_keep) {
  genes<-head(rownames(vst),i)
  tmp<-GetAssayData(object = hnscc, slot = "data")[rownames(hnscc)%in% genes,]
  print(dim(tmp))
  result[[length(result) + 1]]<-tmp
  names(result)[length(result)] <- paste("top", i, sep = "_")
}
  exprimentID="hpv+hnscc"
  #obtaining the matrix
  save(result,file=paste("results/", exprimentID,"_topgenesmatrix.RData",sep=""))


  library(Seurat)
  library(irlba)
  library(rflann)
  library(e1071)
  GapClust <- function(data, k=200){

    rare_cells<-list()

    tmp<-data
    cat(dim(tmp))
    pca <- irlba(t(tmp), nv=min(c(50, dim(tmp)-1))) # More robust no error, contrast to calcul.pca
    pca$pca <-t(pca$d*t(pca$u))
    knn.res <- Neighbour(pca$pca, pca$pca, k=k)

    distance.diff <- (knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE])
    diff.left <- distance.diff[, -1, drop = FALSE] - distance.diff[, -ncol(distance.diff), drop = FALSE]
    diff.both <- diff.left[, -ncol(diff.left), drop=FALSE] - diff.left[, -1, drop=FALSE]
    diff.both[,1] <- diff.both[,1] + distance.diff[,1]  # Very important due to distance variation to the first neighbor.

    v1.k <- matrix(NA, dim(data)[2], k-3)
    skew <- c()
    top.values.ave <- c()
    for(j in 1:dim(diff.both)[2]){
      v <- diff.both[,j]
      v1 <- v
      for(m in 1:length(v)){
        v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
      }
      v1.k[, j] <- (v1)
      v2 <- v1[order(v1, decreasing = T)[(j+2):length(v1)]]
      v2[is.na(v2)] <- 0
      top.values <- v1[knn.res$indices[which.max(v1),1:(j+1)]]
      v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
      skew <- c(skew, skewness(v2))
      top.values.ave <- c(top.values.ave, mean(top.values))
    }

    ids <- which(skew > 2)
    col.mat <- matrix(0, length(ids), dim(tmp)[2])
    for(i in 1:length(ids)){
      top.cell <- which.max(v1.k[,(ids[i])])
      col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]] * top.values.ave[ids[i]]
    }

    id.max <- apply(col.mat, 2, which.max)
    max.val <- apply(col.mat, 2, max)
    id.max[max.val==0] <- 0
    cnt <- table(id.max)
    cnt <- cnt[names(cnt)!='0']
    id.max.match <- cnt[which(cnt == (ids[as.integer(names(cnt))] + 1))] - 1

    cls <- rep(0, dim(tmp)[2])
    for(id.match in id.max.match){
      cls[id.max==(id.match)] <- which(id.max.match %in% id.match)
    }

    rare.cells <- list()
    for(id.match in id.max.match){
      rare.cells[[as.character(id.match)]] <- knn.res$indices[which.max(v1.k[,id.match]), 1:(id.match+1)]
    }
    results <- list(skewness=skew, rare_cell_indices=rare.cells, rare_score=v1.k,pca=pca$u)
    cal<-unlist(results$rare_cell_indices)
    #cal<-colnames(hnscc)[cal]
    #rare_cells<-cal
    return(results)
  }

  #### iterating over the gapclust function to identify the rare cells
  gapclust_rarecells<-list()
  for (i in seq_along(result)) {
    # Call the GapClust function for the current element in 'test'
    gapclust_rarecells[[length(gapclust_rarecells)+ 1]] <- GapClust(result[[i]])
    cat("Iteration:", i, "\n")
  }


  names(gapclust_rarecells) <- paste("top", gene_keep, sep = "_")

  save(gapclust_rarecells,file=paste("results/", exprimentID,"_overlaprarecells.RData",sep=""))

  library(ggplot2)
  library(tidyr)
  long_data <- pivot_longer(overlapgapclust_hnscc, cols = starts_with("top."), names_to = "TopGenes", values_to = "OverlapPercentage")
  ggplot(long_data, aes(x = OverlapPercentage, y = as.numeric(gsub("top\\.", "", TopGenes)), color = TopGenes)) +
    geom_line() +
    geom_point() +
    labs(title = "Percentage of Overlapping Rare Cells gapclust vs. Top Genes",
         x = "Overlap Percentage",
         y = "Top Genes") +
    theme_minimal()

  # retrive the pca values

  pca<-list()

  for (i in seq_along(gapclust_rarecells)){
    element<-gapclust_rarecells[[i]]$pca
    element<-t(element)
    pca[[length(pca)+ 1]] <- element
  }

  names(pca)<-paste("top",gene_keep,sep = "_")





