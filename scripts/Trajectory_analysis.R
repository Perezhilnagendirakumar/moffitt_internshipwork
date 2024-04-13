# script to perform trajectory analysis


set.seed(1234)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)


# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3


# ...1 Convert to cell_data_set object ------------------------

cds <- as.cell_data_set(Bcells)
cds

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info

list_cluster <- Bcells@meta.data$SCT_snn_res.0.1
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- Bcells@reductions$umap@cell.embeddings



# plot

cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
           color_cells_by = "ident",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)
pdf("trajectory.pdf", height = 13, width = 7.5)
plot_cells(cds,
           color_cells_by = 'ident',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)
dev.off()


# ...4. Order the cells in pseudotime -------------------

root_group = colnames(cds)[pData(cds)$ident== "Germinal Center B cells"]

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells =root_group)

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) +
  geom_boxplot()




# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

FeaturePlot(b.seu, features = c('E2F2', 'STMN1', 'CD52'))


# visualizing pseudotime in seurat

b.seu$pseudotime <- pseudotime(cds)
Idents(b.seu) <- b.seu$redefined_cluster
FeaturePlot(b.seu, features = "pseudotime", label = T)



