# scripts for performing enrichment analysis for particular modules M1 and M2 that depicts the 
# minor cell population from the rare cell tools 
library(hdWGCNA)
library(WGCNA)
library(enrichR)

dbs <-c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021','LINCS_L1000_Chem_Pert_up','LINCS_L1000_Chem_Pert_down', 'WikiPathway_2021_Human', 'KEGG_2021_Human')
Bcells <- RunEnrichr(Bcells, dbs=dbs)
enrichr_df <- GetEnrichrTable(Bcells) %>% subset(P.value < 0.05)
pdf("GO_BP.pdf",height = 15, width = 15)
EnrichrDotPlot(
  Bcells,
  mods ="all", 
  database = "GO_Biological_Process_2021", 
  n_terms=1 
)
dev.off()
