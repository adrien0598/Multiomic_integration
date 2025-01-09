
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
source("scripts/CustomFunctions.R")

# PKN
net <- decoupleR::get_collectri(organism = 'mouse', 
                                split_complexes = FALSE)

# PseudoBulk
NAME = "GL261"
# NAME = "CT2A"
PB = paste0("../data/PseudoBulk_ct_", NAME, ".tsv")
counts <- read.table(PB, header = TRUE) %>%
  select(-grep("RES_D7", colnames(counts)))

# Design
if (NAME == "GL261"){
  design = data.frame(sample = colnames(counts),
                      condition = substring(colnames(counts), 10, 12))
  ct = colnames(counts) %>% sapply(function(x){
    substring(x, 15, nchar(x))
  }) %>% unname()
}else {
  design = data.frame(sample = colnames(counts),
                      condition = substring(colnames(counts), 9, 11))
  ct = colnames(counts) %>% sapply(function(x){
    substring(x, 16, nchar(x))
  }) %>% unname()
}

## DEG and TF function for a given cell type
Run_deg = function(cell_type, cell_line){
  counts_tmp = counts[,ct == cell_type]
  gene.size = counts_tmp %>% apply(1, sum)
  counts_tmp = counts_tmp[gene.size > 25,]
  design_tmp = design[ct == cell_type,]
  fit <- vsnMatrix(as.matrix(counts_tmp))
  meanSdPlot(fit)
  counts_vsn <- as.data.frame(vsn::predict(fit,as.matrix(counts_tmp)))
  limmaRes <- runLimma(counts_vsn, design_tmp, comparisons = list(c(2,-1)))
  table <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 20000, adjust.method = "fdr"))
  
  deg <- table %>%
    dplyr::select(logFC, t, P.Value) %>% 
    dplyr::filter(!is.na(t)) %>% 
    as.matrix()
  
  deg %>% as.data.frame() %>% saveRDS(paste0("data/DEGct/DEG_", cell_line, "_", cell_type, ".rds"))
  
  n_tfs = 12 
  
  contrast_acts <- decoupleR::run_ulm(mat = deg[, 't', drop = FALSE], 
                                      net = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor='mor', 
                                      minsize = 10)
  
  f_contrast_acts <- contrast_acts %>%
    dplyr::mutate(rnk = NA)
  
  msk <- f_contrast_acts$score > 0
  
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  
  tfs <- f_contrast_acts %>%
    dplyr::arrange(rnk) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs) %>%
    filter(p_value < 0.01)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p <- ggplot2::ggplot(data = f_contrast_acts, 
                       mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                              y = score)) + 
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    ggplot2::scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1, 
                                                       size = 10, 
                                                       face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 10, 
                                                       face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    ggplot2::xlab("TFs")+
    labs(x = "Transcription Factors", y = "Predicted activity score", title = cell_type)
  
  print(p)
  paste0("figures/", cell_line, "TF_activity_", cell_type, ".png") %>%
    ggsave(plot = p,
           width = 10,
           height = 7,
           unit = "cm")
  
  contrast_acts %>% 
    as.data.frame() %>% 
    filter(condition == "t") %>% 
    saveRDS(paste0("data/TFS_ct/", cell_line, "TF_activity_", cell_type, ".rds"))
}

# GL261
for (cell_t in unique(ct)){
  Run_deg(cell_t, cell_line = NAME)
}



