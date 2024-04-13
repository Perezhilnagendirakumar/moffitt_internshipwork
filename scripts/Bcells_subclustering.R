# This script is used to cluster cancer cells and run trajectory analysis of
# individual patient using Monocle2 and CytoTRACE, followed by GeneSwitches analysis,
# in order to identify pre-nodal metastatic cells and modulators.

library(Seurat)
library(sctransform)
library(clustree)
library(ggplot2)
library(RColorBrewer)
library(dplyr)





##### ReClustering of Bcells using Seurat #####
load("seurat_integrated.RData")
Bcells <- subset(seurat_integrated, idents = "B cells")
# UMAP of B cell subset 
DimPlot(Bcells, reduction = "umap", label = TRUE, label.size = 4,cols = "#808000")

Bcells <- SCTransform(Bcells, vars.to.regress = c("nCount_RNA","mitoRatio"), verbose = TRUE)

#Perform linear dimensional reduction
npc <- 20
Bcells <- RunPCA(Bcells, features = VariableFeatures(object = Bcells), npcs = npc)
pdf("heatmap_PC.pdf", height = 13, width = 7.5)
DimHeatmap(Bcells, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()
upc <- 15

Bcells <- FindNeighbors(Bcells, dims = 1:upc)
for (i in seq(0,1,by = 0.1)) {
  Bcells <- FindClusters(Bcells, resolution = i, n.start = 10)
}
pdf("clustree.pdf", height = 10)
clustree(Bcells, prefix = "SCT_snn_res.")
dev.off()

Bcells<-FindNeighbors(Bcells)
Bcells <- FindClusters(Bcells, resolution = 0.1)
Bcells <- RunUMAP(Bcells, dims = 1:15, min_dist = 0.3)

colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", 
            "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", 
            "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", 
            "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", 
            "#5254a3")
##plot umap
p1 <- DimPlot(Bcells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4,cols = c("Blue","Red","Green","Purple","Orange"))
p2<- DimPlot(Bcells, reduction = "umap", group.by = "orig.ident",cols = colors)
p3 <- DimPlot(Bcells, reduction = "umap", group.by = "group", label = FALSE,cols = c("Blue","red"))

pdf("umap_Bcells_sct0.1.pdf", height = 10, width = 13)
CombinePlots(plots = list(p1, p2, p3))
dev.off()


# find markers for every cluster compared to all remaining cells,
# report only the positive ones
Bcell.markers <- FindAllMarkers(Bcells, assay="SCT",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
Bcell.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(Bcell.markers, "marker_genes_cluster.csv")
top10 <- Bcell.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("heatmap_degenes.pdf",height = 15, width = 15)
DoHeatmap(Bcells, features = top10$gene) + NoLegend()
dev.off()

saveRDS(Bcells,"Bcells_seurat.rds")

# perform DE test using scp 
library(SCP)

Bcells <- RunDEtest(srt = Bcells, group_by = "seurat_clusters", fc.threshold = 1.2, only.pos = FALSE)
VolcanoPlot(srt = Bcells, group_by = "seurat_clusters")

# GSEA 
library(fgsea)

Bcell_markers <- Bcell.markers[order(Bcell.markers$avg_log2FC, decreasing = T),]
Bcell_markers <- Bcell_markers[,c("gene", "avg_log2FC")]
rownames(Bcell_markers) <- NULL

prepare_ranked_list <- function(ranked_list) {
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

Bcell_markers <- prepare_ranked_list(Bcell_markers)

hallmark_pathway <- gmtPathways("h.all.v7.0.symbols.gmt.txt")

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results <- fgsea(pathways = hallmark_pathway,
                       stats = Bcell_markers,
                       minSize = 15,
                       maxSize = 500,
                       nperm= 1000)
library(dplyr)
fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

# clusterprofiler 

library(clusterProfiler)

library(org.Hs.eg.db)
Bcell_markers<-gene_diff_score
map = mapIds(org.Hs.eg.db, keys = names(gene_diff_score), keytype = "SYMBOL", column = "ENTREZID")
library(GSEAtraining)
gene_diff_score = convert_to_entrez_id(gene_diff_score)
#Gene ontology = BP 
gse = gseGO(geneList = gene_diff_score, ont = "BP", OrgDb = org.Hs.eg.db)
head(gse)

# plotting the results Gene set enrichment 
require(DOSE)
#dotplot 
dotplot(gse)
## remove redundent GO terms
gse2<-simplify(gse)


# heatmap plot 

A<-heatplot(gse2)
B<-heatplot(gse2, foldChange=gene_diff_score)

library(gridExtra)
library(grid)

grid.arrange(
  arrangeGrob(A, top = textGrob("A", x = 0, hjust = -0.5, vjust = 1, gp = gpar(fontsize = 12, fontface = "bold"))),
  arrangeGrob(B, top = textGrob("B", x = 0, hjust = -0.5, vjust = 1, gp = gpar(fontsize = 12, fontface = "bold"))),
  ncol = 2
)

# cnet plot 
C<-cnetplot(gse2,foldChange =gene_diff_score)
D<-cnetplot(gse2, circular=TRUE, foldChange =gene_diff_score,colorEdge=TRUE)

grid.arrange(rasterGrob(c), rasterGrob(d), nrow = 2)


# KEGG pathway analysis 

gse = gseKEGG(geneList = gene_diff_score, ont = "BP", OrgDb = org.Hs.eg.db)

# plotting the results Gene set enrichment 
require(DOSE)
#dotplot 
dotplot(gse)

# heatmap plot 
## remove redundent GO terms
gse2<-simplify(gse)

A<-heatplot(gse2)
B<-heatplot(gse2, foldChange=gene_diff_score)
grid.arrange(rasterGrob(c), rasterGrob(d), nrow = 2)



#ACtivated Bcells /ABC genes

pdf("ActivatedBcell_marker_gene_exp.pdf")
FeaturePlot(Bcells, 
            features = c("CD19", "MS4A1", "ID3", "GPR183"),
            ncol = 4, 
            label = TRUE, 
            cols = c("grey", "red")) +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "plain"),
        axis.title = element_text(size = 15, face = "plain"),
        axis.title.y.right = element_text(size = 12),
        legend.text = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 10, face = "plain"),
        axis.line = element_line(size = 0.5))
dev.off()
# Antibody Secreting cells /ASC

pdf("AntibodysecretingBcell_marker_gene_exp.pdf")
FeaturePlot(Bcells, 
            features = c("CD38", "MZB1", "XBP1", "PRDM1"),
            ncol = 4, 
            label = TRUE, 
            cols = c("grey", "red"),split.by = "group") +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "plain"),
        axis.title = element_text(size = 15, face = "plain"),
        axis.title.y.right = element_text(size = 12),
        legend.text = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 10, face = "plain"),
        axis.line = element_line(size = 0.5))
dev.off()

# Germinal B cells /GCB
pdf("GerminalBcell_marker_gene_exp.pdf")
FeaturePlot(Bcells, 
            features = c("TCL1A", "CXCR5", "AICDA", "MME"),
            ncol = 4, 
            label = TRUE, 
            cols = c("grey", "red"),split.by = "group") +
  theme(text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face = "plain"),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "plain"),
        axis.title = element_text(size = 15, face = "plain"),
        axis.title.y.right = element_text(size = 12),
        legend.text = element_text(size = 12, face = "plain"),
        legend.title = element_text(size = 10, face = "plain"),
        axis.line = element_line(size = 0.5))

dev.off()

# rename Idents 
Bcells <- RenameIdents(object = Bcells,
                       "0" = "ABC",
                       "1" = "GCB",
                       "2" = "ASC",
                       "3" = "ASC",
                       "4" = "ABC")

new_labels <- c("Activated B cells", "Germinal Center B cells", "Antibody Secreting Cells")
Idents(Bcells) <- plyr::mapvalues(Idents(Bcells), 
                                              from = c("ABC", "GCB", "ASC"),
                                              to = new_labels)

pdf( "UMAP_Bcell_suptypes.pdf")

DimPlot(object = Bcells,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = TRUE,cols = colors)
dev.off()

saveRDS(Bcells,"Bcells_seu.rds")

