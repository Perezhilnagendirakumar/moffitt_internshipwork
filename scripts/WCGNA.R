library(hdWGCNA)
library(WGCNA)

# Run Co expression network analysis 

Bcells <- SetupForWGCNA(
  Bcells,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells 
Bcells <- MetacellsByGroups(
  Bcells,
  group.by = c("seurat_clusters","cell_type"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)

# normalize the counts
Bcells<- NormalizeMetacells(Bcells)
# setup expression matrix 
Bcells<- SetDatExpr(Bcells,assay = "SCT",slot =  "data")

Bcells<- TestSoftPowers(Bcells, networkType = 'signed')
library(tidyverse)
library(cowplot)
library(patchwork)
# plot the results 

plot_list<-PlotSoftPowers(Bcells)

power_table <- GetPowerTable(Bcells)
head(power_table)
# construct wgcna network 
Bcells <- ConstructNetwork(Bcells,setDatExpr=F)
#plot the dendrogram
pdf("dendogram.pdf",height = 15, width = 15)
PlotDendrogram(Bcells)
dev.off()
Bcells <- ScaleData(Bcells,features=VariableFeatures(Bcells))


# compute ME and module connectivity 
Bcells <- ModuleEigengenes(
  Bcells,
  group.by.vars="cell_type"
)


Bcells<- ModuleConnectivity(
  Bcells
)

# Reset the module names 
Bcells <- ResetModuleNames(
  Bcells,
  new_name = "M")

pdf("KMEs.pdf",height = 15, width = 15)
PlotKMEs(Bcells, ncol=5)
dev.off()

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
Bcells <- ModuleExprScore(
  Bcells,
  n_genes = 25,
  method='UCell'
)
# ME feature plot 
pdf("ME_featureplot.pdf",height = 15, width = 15)
plot_list <- ModuleFeaturePlot(
  Bcells,
  features='scores',
  order='shuffle', # order so cells are shuffled
  ucell = TRUE ,
  label=TRUE
)
wrap_plots(plot_list, ncol=6)
dev.off()

ModuleCorrelogram(Bcells)


