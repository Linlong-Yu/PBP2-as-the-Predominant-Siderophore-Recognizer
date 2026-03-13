# 加载包
library(phytools)
library(dplyr)
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
#library(RColorBrewer)
library(Biostrings)
library(phangorn)
library(ggtree)
library(readxl)
alignment <- readAAStringSet("../input/siderophore_rec_MSA_clustalo.fasta")
phy <- as.phyDat(alignment, type = "AA")
dm <- dist.ml(phy, model = "WAG")  # use WAG model
tree <- NJ(dm)
tree <- midpoint.root(tree)
boot <- bootstrap.phyDat(phy, FUN = function(x) NJ(dist.ml(x)), bs=1000)
tree <- plotBS(tree, boot, p=50)
metadata <- read_excel("../input/PBP2_info_experiment.xlsx") %>%
  dplyr::select(
    protein_id = `Protein_name`,  
    species = `Organism name`,
    Gram = `Gram-postive/Gram-negative`,
    ligand = `Siderophore`,
    Genus = `Genus`
  )
metadata$pair = paste(metadata$protein_id,'-',metadata$ligand)
tree_tips <- data.frame(protein_id = tree$tip.label)
tree_metadata <- left_join(tree_tips, metadata, by = "protein_id")
tree_metadata$rec_id = tree_metadata$pair
Gram_colors <- c("#4b8cff", "#ff6b6b")

p <- ggtree(tree) %<+% tree_metadata +
  geom_tiplab(
    aes(label = Genus, color = Gram),  
    size = 3,
    align = TRUE,
    offset = 0.33 
  ) +
  geom_tiplab(
    aes(label = rec_id),
    size = 3,
    align = FALSE,
    color = "black",
    offset = -0.33,
    vjust = -0.5
  ) +
 
  scale_color_manual(
    name = "Gram postive/negative",
    values = Gram_colors
  ) +

  theme_tree2() +
  labs(title = "Phylogenetic Tree of SBP Receptors") +
  xlim(0, 2)  
print(p)
#ggsave("../Figure/phylogenetic_tree.pdf", width = 10, height = 8) 
