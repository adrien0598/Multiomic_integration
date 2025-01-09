
############################ PseudoBulk ########################################
# Build pseudobulk from single cell seurat object
################################################################################

# libraries
library(Seurat)
library(dplyr)
library(edgeR)
library(tibble)

# import Seurat object
NAME = "GL261"
# NAME = "CT2A"
SEURAT_OBJ = "20250106_GL261_final_sharedAnnot.rds"
# SEURAT_OBJ = "20240628_CT2A_deep_annotation_myeloid_homo.rds"
CONDITION = "cond.day"
CELL_TYPE = "cell_types"
IDENT = "orig.ident"
sc <- readRDS(paste0("../data/", SEURAT_OBJ))

# computing pb
counts = sc@assays[["RNA"]]@counts
mice = sc[[IDENT]][,1]
ct = sc[[CELL_TYPE]][,1]
cond = sc[[CONDITION]][,1]
pb = data.frame(geneID = row.names(counts))
mm = paste0(cond, "-", ct)
condition = c()
for (m in (paste0(mice,ct) %>% unique())){
  tmp = counts[,m == paste0(mice,ct)]
  tmp_pb = tmp %>% apply(1, sum)
  pb[[m]] = tmp_pb
  condition = c(condition, unique(mm[m == paste0(mice,ct)]))
}

# TMM normalisation with edgeR
pb = pb %>% 
  column_to_rownames("geneID") 
y = pb %>%
  DGEList(group = condition)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)

for (i in c(1:ncol(pb))){
  pb[,i] = (pb[,i]/ (y$samples$lib.size[i]/y$samples$norm.factors[i]))*1e6
}

# Export
write.table(pb, paste0("../data/PseudoBulk_ct_", NAME, ".tsv"), 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
