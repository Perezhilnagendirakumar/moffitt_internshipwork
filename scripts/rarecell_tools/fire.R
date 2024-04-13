# script for running FIRE rare cell tool

Rfundir              = "/share/lab_teng/trainee/PerezhilNagendirakumar/HNSCC/fire/"
workdir              ="/share/lab_teng/trainee/PerezhilNagendirakumar/HNSCC/fire/"
setwd(workdir)
dir.create(file.path(workdir, "results"), showWarnings = FALSE)
source("preprocess.R")
#source(paste(Rfundir,"preprocess.R"))
project ="B_cells"
data<-Bcells@assays[["RNA"]]@counts
data<-t(data) #sample * fetures 
genes <- c(1:dim(data)[2])
data_mat <- list(mat=data, gene_symbols=genes)
preprocessedList <- ranger_preprocess(data_mat)
preprocessedData <- as.matrix(preprocessedList$preprocessedData)
model <- new(FiRE::FiRE, 100, 50, 1017881, 5489, 0)
model$fit(preprocessedData)
score <- model$score(preprocessedData)
write.csv(score,file=paste("results/",project,"_score.csv",sep=""),quote=F)

q3 <- quantile(score, 0.75)
iqr <- IQR(score)
th <- q3 + (1.5*iqr)

#Select indexes that satisfy IQR-based thresholding criteria.
indIqr <- which(score >= th)
rare_cells <- rownames(data)[indIqr]
p1<-DimPlot(object = Bcells, cells.highlight =rare_cells , group.by ='cell_type',cols.highlight = "yellow", cols = "gray", order = FALSE,label = TRUE,raster=FALSE)

save(rare_cells,file=paste("results/",project ,"_rarecell.RData",sep=""))