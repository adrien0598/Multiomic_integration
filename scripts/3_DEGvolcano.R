
############################ DEGvolcano ########################################
# Plot a volcano/jitter plot combining all deg per cell type
################################################################################

# library
library(tidyverse)
library(ggprism)
library(ggsci)

# Recover DEG
cell_line = "GL261"
# cell_line = "CT2A"
p = 0.05
files = list.files(path = "../data/DEGct/", pattern = cell_line)
gg = NULL
for (f in files) {
  ct = f %>% substring((6+nchar(cell_line)), (nchar(f)-4))
  tmp = readRDS(paste0("../data/DEGct/", f)) %>% 
    filter(P.Value < p)
  gg = gg %>% rbind(data.frame(Gene = tmp$ID, 
                 Cell_type = rep(ct, nrow(tmp)),
                 logFC = tmp$logFC))
}

# select labels to be plotted in the graph
n = 10 # number of labels per cluster (2*n)
gg_selection = NULL
for (f in files) {
  ct = f %>% substring((6+nchar(cell_line)), (nchar(f)-4))
  tmp = readRDS(paste0("../data/DEGct/", f)) %>% 
    filter(P.Value < p) 
  up = tmp %>%
    arrange(-logFC) %>%
    slice(c(1:n))
  down = tmp %>%
    arrange(logFC) %>%
    slice(c(1:n))
  gg_selection = gg_selection %>% 
    rbind(data.frame(Gene = c(up$ID, down$ID), 
                     Cell_type = rep(ct, (2*n)),
                     logFC = c(up$logFC, down$logFC)))
}

# Order by number of DEG
cluster_label = gg %>% 
  count(Cell_type) %>% 
  arrange(-n) %>% 
  select(Cell_type)
gg = gg %>% 
  mutate(Cell_type = factor(gg$Cell_type, levels = cluster_label[,1]))


ggplot(gg) +
  aes(x = Cell_type, y = logFC, color = Cell_type)+
  geom_jitter(size = .5, alpha = .8)+
  geom_hline(yintercept = 0.75)+
  geom_hline(yintercept = -0.75)+
  geom_label_repel( 
    data= gg_selection, # Filter data first
    aes(label=Gene))+
  theme_prism()+
  scale_color_npg()

ggsave(paste0("../figures/DEGvolcano_", cell_line, ".png"),
       height= 3*730, width = 3*1000, unit = "px", dpi = 3*100)