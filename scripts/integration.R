# Script to perform integration

library(Seurat)
library(ggplot2)
library(Rcpp)
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
library(markdown)
library(XML)
library(RCurl)
library(AnnotationHub) # BiocManager::install("AnnotationHub")
library(gridExtra)
library(cowplot)
set.seed(1234)


# read in data -------------------
split_seurat <- readRDS('data/split_seurat.rds')

# Select the most variable features to use for integration ---------------
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000)

# run PCA on each object in the list, which is required for running the alternative reciprocal PCA workflow --------------------
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- ScaleData(x, features = integ_features, verbose = FALSE)
  x <- RunPCA(x, features = integ_features, verbose = FALSE)
})

# Prepare the SCTransform object for integration -------------------------
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)



# Perform CCA, find the best buddies or anchors and filter incorrect anchors --------------------
# CCA is prohibitively computationally expensive, hence using rpca reduction for better runtimes
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                        normalization.method = "SCT",
                                        reduction = 'rpca',
                                        anchor.features = integ_features)

# Integrate across samples ------------------------------------
seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")

# Save integrated seurat object ---------------------------
#saveRDS(seurat_integrated, "results/integrated_seurat.rds")
#seurat_integrated <- readRDS("results/integrated_seurat.rds")

# UMAP visualization - to visualize integrated data -------------------
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated,dims = 1:50)
pdf(filename = "PC_heatmap.pdf", width = 16, height = 8.135, units = "in", res = 300)
DimHeatmap(seurat_integrated,
           dims = 1:11,
           cells = 500,
           balanced = TRUE)
dev.off()
# Visualize the PC and select the appropriate PC for downstream analysis
pdf(filename = "PC_select.pdf", width = 16, height = 8.135, units = "in", res = 300)
ElbowPlot(object = seurat_integrated,
          ndims = 50)
dev.off()
pdf(filename = "pc_top_gene.pdf", width = 16, height = 8.135, units = "in", res = 300)
VizDimLoadings(seurat_integrated, dims = 1:2, reduction = "pca")
dev.off()

seurat_integrated<- FindNeighbors(seurat_integrated, dims = 1:11)
seurat_integrated <- FindClusters(merged_seurat,resolution = c(0.3, 0.5, 0.7, 0.9, 1.1))
Idents(object = seurat_integrated) <- "integrated_snn_res.0.5"

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated,
                             dims = 1:11,
                             reduction = "pca")

custom_palette <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
                    "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd",
                    "#5e4fa2", "#9e0142", "#d53e4f", "#f46d43", "#fdae61",
                    "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5")

samples <- seurat_integrated@meta.data$orig.ident

num_samples <- length(unique(samples))
# Generate a vector of random colors
set.seed(123)
random_colors <- rainbow(num_samples)

group_colors <- c("red", "blue")

# Plot UMAP
d1 <- DimPlot(seurat_integrated,reduction = "umap",label = TRUE,cols = custom_palette)
ggsave(d1, filename = 'figures/UMAP.pdf', width = 10, height = 10)
d2 <- DimPlot(seurat_integrated,reduction = "umap",group.by= "orig.ident" ,cols = random_colors)
ggsave(d2, filename = 'figures/UMAP_by_sample.pdf', width = 10, height = 10)
d3<-DimPlot(seurat_integrated, reduction = "umap",group.by= "group",cols = group_colors)
ggsave(d3, filename = 'figures/UMAP_by_group.pdf', width = 10, height = 10)

#plot TSNE
seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:11,check_duplicates = FALSE)
d4<-DimPlot(seurat_integrated, reduction = "tsne",label = TRUE,cols = custom_palette)
ggsave(d4, filename = 'figures/TSNE.pdf', width = 10, height = 10)
d5<-DimPlot(seurat_integrated, reduction = "tsne",group.by = "group",cols = group_colors)
ggsave(d5, filename = 'figures/TSNE_by_group.pdf', width = 10, height = 10)
d6<-DimPlot(seurat_integrated, reduction = "tsne",group.by= "orig.ident",cols = random_colors)
ggsave(d6, filename = 'figures/TSNE_by_sample.pdf', width = 10, height = 10)

# Save integrated Seurat Object
saveRDS(seurat_integrated, "results/integrated_seurat.rds")

#Plot the PCA colored by group
pdf(filename = "PCA_by_group.pdf", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "group",
        cols =group_colors)
dev.off()
pdf(filename = "PCA_by_sample.pdf", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated,
        reduction = "pca",
        group.by= "orig.ident",
        cols =random_colors)
dev.off()


# function to get the cell proportion based on group for each cluster
get_group_proportions <- function (seuratobj, group.by = "active.ident") {
  if (group.by == "active.ident") {
    seuratobj[["active.ident"]] <- seuratobj@metadata$seurat_clusters
  }
  # gets the total number of cells within each group
  total_populations <- seuratobj@meta.data %>% group_by(group) %>% summarize (total= n())
  # gets the proportion of cells for each clusters within a group by dividing by the total
  count_populations <- seuratobj@meta.data %>% group_by_at(vars(group.by, "seurat_clusters")) %>% summarize (count = n())
  count_populations <- left_join(count_populations, total_populations, by = "group")
  count_populations <- count_populations %>% mutate (ratio = n/total.pop)
  count_populations
}

group_prop<-get_group_proportions(seurat_integrated)
write.table(group_prop,"cell_ratio_by_clusters.pdf",sep="\t")


# function to plot the group proportion 

plot_group_proportions <- function (seuratobj, graph.type = "dodge") {
  # get the proportions using get_group_proportions()
  count_populations <- get_group_proportions(seuratobj)
  # plot the proportions
  if (graph.type == "dodge"){
    ggplot(count_populations, aes (x = active.ident, y = proportion))+
      geom_bar (aes(fill = orig.ident), stat = "identity", position = "dodge")+
      theme(axis.text.x = element_text(angle=45, hjust = 1))
  } else if (graph.type == "stacked") {
    ggplot(count_populations, aes (x = orig.ident, y = proportion))+
      geom_bar (aes(fill = active.ident), stat = "identity", position = "fill")
  }
  else
    print("invalid graph type")
}
pdf(filename = "cell_ratio_by_clusters.pdf")
plot_group_proportions(seurat_integrated)
dev.off()
