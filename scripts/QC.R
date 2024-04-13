# Script to preprocess and QC snRNASeq data


library(Seurat)
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
library(DoubletFinder) # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
library(markdown)
set.seed(1234)


files = list.files()
# read in count data and create Seurat objects --------------------------------------------------------------------------------------

# Create a Seurat object for each sample
for (file in files) {
  print(paste0("Reading file: ", file))

  # Read 10X data into Seurat object
  seurat_data <- Read10X(data.dir = paste0("~/Desktop/assessment/seurat/data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 500,
                                   min.cells = 10,
                                   project = file)

  # Assign Seurat object to a variable with the file name
  assign(file, seurat_obj)
}




# Create a merged Seurat object
# This will make it easier to run the QC steps for both sample groups together and enable us to easily compare the data quality for all the samples.
merged_seurat <- merge(x =GSM4138111 ,
                       y = c(GSM4138113,
                             GSM4138115,
                             GSM4138117,
                             GSM4138119,
                             GSM4138121,
                             GSM4138123,
                             GSM4138125,
                             GSM4138127,
                             GSM4138129,
                             GSM4138131,
                             GSM4138133,
                             GSM4138135,
                             GSM4138137,
                             GSM4138139,
                             GSM4138141,
                             GSM4138143,
                             GSM4138145,
                             GSM4138147,
                             GSM4138149,
                             GSM4138151,
                             GSM4138153,
                             GSM4138155,
                             GSM4138157,
                             GSM4138159,
                             GSM4138161
                       ),
                       #Because the same cell IDs can be used for different samples, we add a sample-specific prefix
                       # to each of our cell IDs using the add.cell.id argument.
                       add.cell.id = files)
# keep only merged seurat object in workspace 
rm(list = setdiff(ls(), "merged_seurat"))
metadata<-merged_seurat@meta.data

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# create sample column
metadata<-merged_seurat@meta.data
# create group column  column
metadata$group <- ifelse(metadata$orig.ident %in% c("GSM4138111", "GSM4138113", "GSM4138115", "GSM4138117", "GSM4138119",
                                                    "GSM4138121", "GSM4138123", "GSM4138125", "GSM4138127", "GSM4138129",
                                                    "GSM4138131", "GSM4138133", "GSM4138135", "GSM4138137", "GSM4138139",
                                                    "GSM4138141", "GSM4138143", "GSM4138145"), "HPV_neg", "HPV_pos")
merged_seurat@meta.data$group <- metadata$group

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")


# Assessing Quality metrics --------------------------------------------------------------------------------------



# 1. cell counts
# Visualize the number of cell counts per group
library(ggplot2)
p1 <- metadata %>%
  ggplot(aes(x=group, fill=group)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

p2 <- metadata %>%
  ggplot(aes(x=group, fill=orig.ident)) +
  geom_bar(position = 'stack') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

p1 + p2


# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500, that is the low end of what we expect.
# If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.

# Visualize the number UMIs/transcripts per cell
metadata %>%
  ggplot(aes (x=nCount_RNA, fill= group)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  facet_wrap(~seq_folder) +
  ylab("Cell density")




# genes detected per cell
# Visualize the distribution of genes detected per cell
metadata %>%
  ggplot(aes(x=nFeature_RNA, fill= group)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  facet_wrap(~seq_folder) +
  scale_x_log10()


# Visualize the distribution of genes detected per cell via boxplot
metadata %>%
  ggplot(aes(x=group, y=log10(nFeature_RNA), fill=group)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NGenes: HPV + and HPV -")

# UMIs vs. genes detected
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata <- merged_seurat@meta.data
metadata %>%
  filter(orig.ident %in% samples_to_plot) %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) +
  geom_point(alpha = 0.4) +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~seq_folder)

# Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot.
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs.


# mitochondrial counts ratio
# This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells.
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>%
  ggplot(aes( x=mitoRatio, fill=group)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 0.2)
# High levels of expression from mitochondria indicate dying or dead cells. Basically, quality samples are those that surpass 0.2 mitochondria ratio mark.


# complexity
# We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score.
# The novelty score is computed by taking the ratio of nGenes over nUMI.
# we expect the novelty score to be above 0.80 for good quality cells.
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, fill=group)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
# Almost all our cells seem to have a good complexity.

preQC_meta <- metadata

# Filtering --------------------------------------
# ...Cell-level Filtering ---------

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, subset = nFeature_RNA > 500
                          & nFeature_RNA <=3000
                          & nCountRNA>=1000
                          & log10GenesPerUMI > 0.80
                          & mitoRatio < 0.10
                          )


# ...Gene-level Filtering ---------
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene (to keep genes expressed in 10 or more cells)
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")
