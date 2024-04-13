source(paste(Rfundir,"GiniClust2_preprocess.R",sep=""))
source(paste(Rfundir,"GiniClust2_filtering_RawCounts.R",sep=""))
source(paste(Rfundir,"GiniClust2_fitting.R",sep=""))
source(paste(Rfundir,"GiniClust2_Gini_clustering.R",sep=""))
table(P_G) #P_G is the Gini-based clustering result
source(paste(Rfundir,"GiniClust2_Gini_tSNE.R",sep="")) #visualization of GiniClust results using tSNE