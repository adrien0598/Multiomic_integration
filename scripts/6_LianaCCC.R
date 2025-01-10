
############################ LianaCCC ########################################
# Preparing data for Liana+ run (python)
# Analysing output from Liana+
################################################################################

# libraries
library(reticulate)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(tidyverse)
library(ggsci)
library(ggprism)

# Import conversion table for Mouse -> Human genes
traduction_table = read.csv("../data/mouse_human_marmoset_macaque_orthologs_20231113.csv") %>%
  select(c(human_Symbol, mouse_Symbol)) %>% 
  na.omit()

# import Seurat object
NAME = "GL261"
# NAME = "CT2A"
SEURAT_OBJ = "20250106_GL261_final_sharedAnnot.rds"
# SEURAT_OBJ = "20240628_CT2A_deep_annotation_myeloid_homo.rds"
sc <- readRDS(paste0("../data/", SEURAT_OBJ))
sc <- subset(sc, cond.day != "RES_D7")
sc$cell_types = as.character(sc$cell_types)

# Add tranlated count table to object
counts = sc@assays$SCT@counts
counts = counts[(row.names(counts) %in% traduction_table$mouse_Symbol),]
row.names(counts) = row.names(counts) %>% sapply(function(x){
  traduction_table$human_Symbol[traduction_table$mouse_Symbol == x]
}) %>% unname()
Hum_gene_assay <- CreateAssayObject(counts)
sc[['hum_gene_assay']] <- Hum_gene_assay

# Export for Liana+
SaveH5Seurat(sc, filename = paste0("../data/", NAME, "_anndata.h5Seurat"))
Convert(paste0("../data/", NAME, "_anndata.h5Seurat"), dest = "h5ad", assay = "hum_gene_assay")
writeMM(sc@assays$hum_gene_assay@counts, paste0("../data/SCT_", NAME, "_counts.mtx"))

################################################################################
########################### Running Liana.py ###################################

# Loading Liana results
lr_res = read.csv(paste0("../data/lr_res_", NAME, ".csv")) %>% 
  filter(ligand_padj < 0.05) %>%
  filter(receptor_padj < 0.05) %>%
  filter(interaction_stat > 0)

# Plot general results
lr_res %>% 
  select(source, target) %>% 
  pivot_longer(cols = c(source, target)) %>% 
  ggplot()+
  aes(x = name, fill = value)+
  geom_bar()+
  scale_fill_npg()+
  theme_prism()
ggsave(paste0("../figures/", NAME, "_CCC.png"), 
       height = 2*400, width = 2*423, dpi = 2*100, 
       unit = "px")

# Other plot
s = lr_res %>% 
  select(ligand) %>% 
  count(ligand) %>% 
  arrange(-n)
s = s$ligand
lr_res %>% 
  select(ligand, source) %>% 
  count(ligand, source) %>% 
  arrange(-n) %>% 
  mutate(ligand = factor(ligand, levels = s)) %>%
  ggplot()+
  aes(x = ligand, y = n, fill = source)+
  geom_col()+
  theme_prism()+
  scale_fill_npg()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

s = lr_res %>% 
  select(receptor) %>% 
  count(receptor) %>% 
  arrange(-n)
s = s$receptor
lr_res %>% 
  select(receptor, target) %>% 
  count(receptor, target) %>% 
  arrange(-n) %>% 
  mutate(receptor = factor(receptor, levels = s)) %>%
  ggplot()+
  aes(x = receptor, y = n, fill = target)+
  geom_col()+
  theme_prism()+
  scale_fill_npg()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Export as RDS
lr_res %>% 
  saveRDS(paste0("../data/", NAME, "_CCC.rds"))
