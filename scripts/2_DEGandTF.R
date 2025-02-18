
############################ DEGandTF ########################################
# Perform differential analysis of gene expression 
# Compute TF activity prediction
################################################################################

# libraries
library(decoupleR)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(limma)
library(vsn)
library(ggprism)
source("CustomFunctions.R")

# PKN
net <- decoupleR::get_collectri(organism = 'mouse', 
                                split_complexes = FALSE)

# PseudoBulk
#NAME = "GL261"
NAME = "CT2A"
PB = paste0("../data/PseudoBulk_ct_", NAME, ".tsv")
counts <- read.table(PB, header = TRUE) 
counts = counts[,!grepl("RES_D7", colnames(counts))]

# Design
if (NAME == "GL261"){
  design = data.frame(sample = colnames(counts),
                      condition = substring(colnames(counts), 10, 12))
  ct = colnames(counts) %>% sapply(function(x){
    substring(x, 16, nchar(x))
  }) %>% unname()
}else {
  design = data.frame(sample = colnames(counts),
                      condition = substring(colnames(counts), 9, 11))
  ct = colnames(counts) %>% sapply(function(x){
    substring(x, 15, nchar(x))
  }) %>% unname()
}

## DEG and TF function for a given cell type
Run_deg = function(cell_type, cell_line, a){
  counts_tmp = counts[,ct == cell_type]
  gene.size = counts_tmp %>% apply(1, sum)
  counts_tmp = counts_tmp[gene.size > a,]
  design_tmp = design[ct == cell_type,]
  print(design_tmp)
  fit <- vsnMatrix(as.matrix(counts_tmp))
  meanSdPlot(fit) # quality control
  counts_vsn <- as.data.frame(vsn::predict(fit,as.matrix(counts_tmp)))
  
  limmaRes <- runLimma(counts_vsn, design_tmp, comparisons = list(c(2,-1)))

  deg <- topTable(limmaRes[[1]], coef = 1, number = 20000, adjust.method = "fdr") %>%
    select(logFC, t, P.Value) %>% 
    rownames_to_column("ID") %>%
    filter(!is.na(t))
  
  deg %>% saveRDS(paste0("../data/DEGct/DEG_", cell_line, "_", cell_type, ".rds"))
  
  n_tfs = 12 
  row.names(deg) = deg$ID
  contrast_acts <- decoupleR::run_ulm(mat = deg[, 't', drop = FALSE], 
                                      net = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor='mor', 
                                      minsize = 10)
  
  f_contrast_acts <- contrast_acts %>%
    mutate(rnk = NA)
  
  msk <- f_contrast_acts$score > 0
  
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  
  tfs <- f_contrast_acts %>%
    arrange(rnk) %>%
    head(n_tfs) %>%
    pull(source)
  
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs) %>%
    filter(p_value < 0.01)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p <- f_contrast_acts %>% ggplot()+
    aes(x = stats::reorder(source, score), y = score) + 
    geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    theme_minimal() +
    theme(axis.title = element_text(face = "bold", 
                                    size = 12),
                   axis.text.x = element_text(angle = 45, 
                                              hjust = 1, 
                                              size = 10,
                                              face = "bold"),
                   axis.text.y = element_text(size = 10, 
                                              face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    labs(x = "Transcription Factors", y = "Predicted activity score", title = cell_type)
  
  print(p)
  paste0("../figures/", cell_line, "TF_activity_", cell_type, ".png") %>%
    ggsave(plot = p,
           width = 10,
           height = 7,
           unit = "cm")
  
  contrast_acts %>% 
    as.data.frame() %>% 
    filter(condition == "t") %>% 
    saveRDS(paste0("../data/TFct/", cell_line, "TF_activity_", cell_type, ".rds"))
}

# run function for all cell type in selected cell line
unique(ct)
Run_deg(cell_type = "Myeloid", cell_line = NAME, a = 4)
Run_deg(cell_type = "T_cells", cell_line = NAME, a = 5)
Run_deg(cell_type = "Tumor", cell_line = NAME, a = 6)
Run_deg(cell_type = "Vascular", cell_line = NAME, a = 5)
Run_deg(cell_type = "NK_cells", cell_line = NAME, a = 12)
Run_deg(cell_type = "Neuroglia", cell_line = NAME, a = 10)
Run_deg(cell_type = "Neutrophils", cell_line = NAME, a = 30)
Run_deg(cell_type = "B_cells", cell_line = NAME, a = 14)




