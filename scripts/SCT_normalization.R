# Script to normalize and perform PCA

library(Seurat)
library(ggplot2)
library(harmony) # 1.0
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
library(cowplot)
set.seed(1234)

#load('data/seurat_filtered.RData')
str(filtered_seurat)

# Exploring Sources of unwanted variation
# 1. Normalize the counts -------------------------
seurat_phase <- NormalizeData(filtered_seurat)
str(seurat_phase)


# 2. Evaluating effects of cell cycle --------------

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv")
cell_cycle_genes <- read.csv(text = cc_file)
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <-annotations %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Score cells for cell cycle
# Assign each cell a score based on its expression og G2/M and S phase markers. The function calculates cell cycle phase scores based on canonical markers
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)


# Next we would like to determine whether cell cycle is a major source of variation in our dataset using PCA
# To perform PCA: Find most variable features > scale data > runPCA
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)
top10 <- head(VariableFeatures(seurat_phase), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_phase)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)


# Perform PCA to evaulate similarities/differences between cell cycle phase
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
cc1 <- DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")+ ggtitle("Grouped all cells together")

cc2 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "Phase",
               split.by = "Phase")+ ggtitle("Split by Cell Cycle Phase")
cc3 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "Phase",
               split.by = "group")+ ggtitle("Split by Source")
plot <- plot_grid(cc1, cc2, cc3,  ncol = 2, nrow = 2)


ggsave(plot, filename = 'figures/dimplot_cellCycle_before_regressing.pdf', width = 20, height = 10)

# Since cells separate entirely by phase, we need to regress out cell cycle scores

seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase))


# Now, a PCA on the variable genes no longer returns components associated with cell cycle
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), nfeatures.print = 10)




# 3. Evaluating effects of mitochondrial gene expression --------------
# Next, we want to evaluate expression of mitochondrial genes' expression
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks=c(-Inf, 0.003867, 0.011655, 0.013933, Inf),
                                     labels=c("Low","Medium","Medium high", "High"))



# Plot the PCA colored by mitoFr

m1 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "mitoFr")+ggtitle("Grouped all cells together")

m2 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "mitoFr",
               split.by = "mitoFr")+ggtitle("Split by Mito Fractions")
m3 <- DimPlot(seurat_phase,
               reduction = "pca",
               group.by= "mitoFr",
               split.by = "group")+ggtitle("Split by Source")
plot <- plot_grid(m1, m2, m3, ncol = 2, nrow = 2)


ggsave(plot, filename = 'figures/dimplot_mitoFr_before_regressing.pdf', width = 20, height = 10)

# since they show different pattern of scatter, mitoRatio will be required to be regressed out

# 4. SCTransform --------------
# SCTransform automatically accounts for cellular sequencing depth by regressing out sequencing depth (nUMIs)
# Split seurat object by condition to regress cell cycle scoring, mitoRatio and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "seq_folder")
options(future.globals.maxSize = 4000 * 1024^2)


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("S.Score","G2M.Score","mitoRatio"))
}


# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")




