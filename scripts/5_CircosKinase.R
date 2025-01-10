
############################ CircosKinase ########################################
# Draw a nice circos plot summarizing kinase activity
################################################################################

# libraries
library(circlize)
library(tidyverse)
library(dichromat)
library(viridis)

# Import data
## Kinase activity
Activity <- read.csv("../data/phospho_kinase_activities.csv") %>%
  select("ID", "score")

## Network
Net <- read.csv("../data/kinase_substrate_network.csv") %>%
  merge(Activity, by.x = "tf", by.y = "ID", all = TRUE)
tmp = c(c(1:4), c(1:15), c(1:20), c(1:4), c(1:6), c(1:4), c(1:4), c(1:11), 
        c(1:8), c(1:10), c(1:7), c(1:5), c(1:6), c(1:7), c(1:6), c(1:4), 
        c(1:9), c(1:5), c(1:7))
Net$Order = tmp

# Color set up
## Blue red
colfunc1 <- colorRampPalette(c("blue", "white"))
colfunc2 <- colorRampPalette(c("white", "red"))
Blue_Red_color = data.frame(n = c(1:1000), colo = c(colfunc1(500), colfunc2(500)))
tmp = (((Net$logFC-min(Net$logFC))/(max(Net$logFC)-min(Net$logFC)))*999) %>% as.integer()
Net$logFC_color = tmp %>% sapply(function(x){
  Blue_Red_color$colo[Blue_Red_color$n == (x+1)]
}) %>% unlist() %>% unname()

# Viridis
vir = rev(viridis::viridis(1000))
viridi = data.frame(n = c(1:1000), colo = vir)
tmp = (((Net$score-min(Net$score))/(max(Net$score)-min(Net$score)))*999) %>% as.integer()
Net$NES_color = tmp %>% sapply(function(x){
  viridi$colo[Blue_Red_color$n == (1000-x)]
}) %>% unlist() %>% unname()


png("../figures/Circos_kinase.png", width = 1500, height = 1500, pointsize = 40) # writing to png

# Draw circos plot
circos.par("track.height" = 0.2, canvas.xlim = c(-1.2,1.2), canvas.ylim = c(-1.2,1.2))
circos.initialize(as.factor(Net$tf), x = Net$Order)

color_nes = Net$NES_color[!duplicated(Net$tf)]
names(color_nes) = Net$tf[!duplicated(Net$tf)]

circos.track(as.factor(Net$tf), y = Net$logFC, panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  n = c(1:CELL_META$xlim[2])
  i = get.cell.meta.data("sector.numeric.index")
  if (i %in% c(1,2,3, 18,19,17,16, 15, 14)){
    circos.text(x = n, 
                y = CELL_META$cell.ylim[2] + 0.4*mm_y(5), 
                labels = Net$target[Net$tf == CELL_META$sector.index],
                facing = "clockwise", cex = 0.5, adj = 0)
  }else{
    circos.text(x = n, 
                y = CELL_META$cell.ylim[2] + 0.4*mm_y(5), 
                labels = Net$target[Net$tf == CELL_META$sector.index],
                facing = "reverse.clockwise", cex = 0.5, adj = 1)
  }
})
circos.trackPoints(Net$tf, Net$Order, Net$logFC, col = Net$logFC_color, pch = 16, cex = 0.5)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = c(1,xlim[2])
  i = get.cell.meta.data("sector.numeric.index")
  circos.rect(c(1), c(0), (breaks[2]), c(1),
              col = color_nes[[i]], border = NA)
  if (i %in% c(1,2,3, 18,19,17,16, 15, 14)){
    circos.text(CELL_META$xcenter, 
                CELL_META$cell.ylim[1] - 0.5*mm_y(5), 
                CELL_META$sector.index,
                facing = "clockwise", cex = 0.8, adj = 1)
  }else{
    circos.text(CELL_META$xcenter, 
                CELL_META$cell.ylim[1] - 0.5*mm_y(5), 
                CELL_META$sector.index,
                facing = "reverse.clockwise", cex = 0.8, adj = 0)
  }
})

dev.off()

