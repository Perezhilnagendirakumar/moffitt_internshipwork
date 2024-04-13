# Marker Identification

library(Seurat)
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
library(XML)
library(RCurl)
library(AnnotationHub) # BiocManager::install("AnnotationHub")
library(multtest) # BiocManager::install('multtest')
library(metap) # install.packages('metap')
library(reactome.db)
library(org.Hs.eg.db)
library(fgsea)
set.seed(1234)


# Identification of conserved markers in all conditions -------------------------

DefaultAssay(seurat_integrated) <- "RNA"

# testing out on one cluster
markers <- FindConservedMarkers(seurat_integrated, ident.1 = 0, grouping.var = "group", verbose = FALSE)
head(markers)

# running on all clusters

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "group")
}

# function call
# Iterate function across all clusters
conserved_markers <- map_dfr(c(0:19), get_conserved)

write.table(conserved_markers,"marker_gene.txt",sep = "\t")

# Extract top 5 markers per cluster
top5 <- conserved_markers %>%
  mutate(avg_fc = (HPV_pos_avg_log2FC +HPV_neg_avg_log2FC) /2) %>%
  group_by(cluster_id) %>%
  top_n(n = 5,
        wt = avg_fc)

# Visualize top 5 markers per cluster
View(top5)

conserved_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# heatmap of top 10 markers in each cluster
pdf("top10_marker_gene_heatmap.pdf")
DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend()
dev.off()


# visualization of top marker genes

pdf("top_marker_gene_exp.pdf")
FeaturePlot(seurat_integrated, features = c("IL7R", "GZMK", "CD79A","STMN1","GNLY","IL32","TIMP1",
                                            "IGHG1","ISG15","PTGDS","MSAA1","CXCL8","FSCN1","CST3","TPSAB1","GZMB",
                                            "TNFRSF4","C1QB","CXCL13","S100A9"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()






# Plot interesting marker gene expression
pdf("B_cell_marker_gene_exp.pdf")
FeaturePlot(seurat_integrated, features = c("CD79A", "IGKC", "IGLC3", "CD19",
                                            "MS4A1"),
                         ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()
pdf("CD8+_effector_T_cells_marker_gene_exp.pdf")
FeaturePlot(seurat_integrated, features = c("CD2", "CD3D", "CD3E", "CD3G",
                                       "CD4","GZMA"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
 dev.off()
pdf("CD4+_effector_T_cells_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("CD2", "CD3D", "CD3E", "CD3G",
                                            "CD8A","GZMk"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()

pdf( "Macrophages_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("CD14", "CD163", "CD68", "FCGR2A",
                                            "CD86","CXCL2"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()

pdf("Mast_cells_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("MS4A2", "TPSAB1", "TPSB2"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()

pdf( "Myeloid_Dendritic_Cell_1_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("LAMP3", "CD40", "CD83","CCR7"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()

pdf( "Naive_T_Cells_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("CD2", "CD3D", "CD3E","CD3G","CCR7","SELL"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()

pdf("Plasma_cells_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("CD79A", "IGKC", "IGLC3","JCHAIN","SDC1","XBP1"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()


pdf( "Plasmacytoid_Dendritic_cells_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("IL3RA", "LILRA4", "CLEC4C"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()


pdf( "Tregs_marker_gene_exp.pdf")

FeaturePlot(seurat_integrated, features = c("CD2", "CD3D", "CD3E","CD3G","CD4","FOXP3"),
            ncol = 4,label = TRUE) & theme(text = element_text(face = "bold"),
                                           axis.text.x=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.text.y=element_text(angle=0, hjust=1, size=12, face="plain"),
                                           axis.title = element_text(size=15,face="plain"),
                                           axis.title.y.right = element_text(size = 12),
                                           legend.text=element_text(size=12, face="plain"),
                                           legend.title=element_text(size=10, face="plain"),
                                           axis.line = element_line(size=0.5))
dev.off()



# Cluster annotation

seurat_integrated <- RenameIdents(object = seurat_integrated,
                                  "0" = "Naive T cells",
                                  "1" = "CD8+ effector T cells",
                                  "2" = "CD8+ effector T cells",
                                  "3" = "B cells",
                                  "4" = "Tregs",
                                  "5" = "Macrophages",
                                  "6" = "CD4+ effector T cells",
                                  "7" = "CD8+ effector T cells",
                                  "8" = "CD8+ effector T cells",
                                  "9" = "Tregs",
                                  "10" = "Macrophages",
                                  "11" = "Macrophages",
                                  "12" = "Plasma cells",
                                  "13" = "CD8+ effector T cells",
                                  "14" = "Plasmacytoid Dendritic Cells",
                                  "15" = "B cells",
                                  "16" = "CD8+ effector T cells",
                                  "17" = "Myeloid Dendritic Cell 1",
                                  "18" = "Macrophages",
                                  "19" = "CD8+ effector T cells"
                                        )

colors <- c("#0000FF", "#FF0000", "#808000", "#FFFF00", "#008080", "#A52A2A", "#ADD8E6", "#800080", "#90EE90", "#FFC0CB")
#UMAP
pdf( "UMAP_cell_type.pdf", width = 16, height = 8.135, units = "in", res = 600)
DimPlot(object = seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = TRUE,cols = colors)
dev.off()

#TSNE
pdf( "tsne_cell_type.pdf", width = 16, height = 8.135, units = "in", res = 600)

DimPlot(object = seurat_integrated,
        reduction = "tsne",
        label = TRUE,
        label.size = 3,
        repel = TRUE)
dev.off()

# Dot plot
markers.to.plot <- c("CD79A","IGKC","IGLC3","CD19","MS4A1","CD2","CD3D","CD3E","CD3G","CD4","GZMA","CD8A",
                     "GZMK","GZMH","NKG7","CD14","CD163","CD68","FCGR2A","CD86","CXCL2","MS4A2","TPSAB1",

                                          "TPSB2","LAMP3","CD40","CD83","CCR7","SELL","IL7R","JCHAIN","SDC1","XBP1","IL3RA","LILRA4","CLEC4C","FOXP3")
pdf("used_marker_gene_dotplot.pdf", width = 16, height = 8.135, units = "in", res = 600)

DotPlot(seurat_integrated, features = markers.to.plot, cols = c("blue", "grey"),
        dot.scale = 8) +
         RotatedAxis()
dev.off()



