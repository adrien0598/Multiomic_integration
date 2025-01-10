
############################ KinaseActivity ########################################
# Compute Kinase / phosphatase activity prediction from phosphoproteomic
################################################################################

# libraries
library(decoupleR)
library(tidyverse)
library(ggrepel)
library(limma)
library(vsn)
library(ggprism)
library(OmnipathR)
source("CustomFunctions.R")

# Import phosphoproteomic data
phospho <- read.table("../data/phospho.csv", sep = ";", dec = ",", header = TRUE)
gg <- phospho %>%
  filter(adj.P.Val < 0.01) %>%
  select(ID, logFC) %>%
  as.data.frame() %>%
  arrange(desc(logFC))
plot_top(gg, 10, "logFC")

# Import PKN
## Kinases
uniprot_kinases <- OmnipathR::annotations(resources = "UniProt_keyword") %>%
  filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
  distinct() %>%
  pull(genesymbol) %>%
  unique() # 636 kinases

## Kinase - phosphosite links
omnipath_ptm <- OmnipathR::signed_ptms() %>%
  filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
  filter(!(stringr::str_detect(sources, "ProtMapper") & n_resources == 1)) %>%
  mutate(p_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset),
                mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
  transmute(p_site, enzyme_genesymbol, mor) %>%
  filter(enzyme_genesymbol %in% uniprot_kinases)

omnipath_ptm$likelihood <- 1

# removing ambiguous modes of regulations
omnipath_ptm$id <- paste(omnipath_ptm$p_site,omnipath_ptm$enzyme_genesymbol, sep ="")
omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id),]
omnipath_ptm <- omnipath_ptm[,-5]

# renaming KSN to fit decoupler format
names(omnipath_ptm)[c(1,2)] <- c("target","tf")

# Compute kinase act.
p = 0.05
phospho_filtered = phospho %>% 
  filter(adj.P.Val < p) %>% 
  select(ID, logFC) %>% 
  mutate(logFC = as.numeric(logFC)) %>%
  tibble::column_to_rownames("ID") %>%
  as.matrix()
kin_activity <- run_wmean(
  minsize = 4,
  mat = phospho_filtered, 
  network = omnipath_ptm, 
  .source = "tf",
  .target = "target",
  times = 1000
)
kin_activity <- kin_activity[kin_activity$statistic == "norm_wmean",c(2,4)] 
colnames(kin_activity)[1] = "ID"
kin_activity = kin_activity %>% arrange(-score)

gg <- kin_activity %>% # ploting top kinases
  select(ID, score) %>%
  as.data.frame() %>%
  arrange(desc(score))
plot_top(gg, 5, "score")

# Export
## Kinase activities
write.csv(kin_activity, "../data/phospho_kinase_activities.csv")

substrat = omnipath_ptm %>% 
  filter(tf %in% kin_activity$ID) %>% 
  filter(target %in% row.names(phospho_filtered))
tmp = phospho[phospho$ID %in% substrat$target,] %>% select(ID, logFC)
substrat = substrat %>% merge(tmp, by.x = "target", by.y = "ID", all = TRUE)

## Kinase phosphosite network
write.csv(substrat, "../data/kinase_substrate_network.csv", row.names = FALSE)


