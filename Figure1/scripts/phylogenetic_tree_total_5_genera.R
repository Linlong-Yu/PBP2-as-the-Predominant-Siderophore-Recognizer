library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(dplyr)

library(ggnewscale)
library(ggnewscale)
tree <- read.tree("../input/RAxML_result.GCF_data_refined.tre")
genus_info <- read.csv("../input/genus_info.csv")
sider_info <- read.csv('../input/summ_genus.csv')
#  get strainname 
sider_info$strain_name <- sub("^(([^_]+_[^_]+)).*$", "\\1", sider_info$dataname)
genus_info$number_NIS <- 0
genus_info$number_NRPS <- 0
for (i in seq_len(nrow(genus_info))) {
  strain <- genus_info$strain_name[i]
  
  
  matched_rows <- sider_info[sider_info$strain_name == strain, ]
  
  
  genus_info$number_NIS[i] <- sum(grepl("NI-siderophore", matched_rows$syn_type, ignore.case = TRUE), na.rm = TRUE)
  
  genus_info$number_NRPS[i] <- sum(grepl("NRP-metallophore", matched_rows$syn_type, ignore.case = TRUE), na.rm = TRUE)
}
genus_info$total_sider_number <- genus_info$number_NIS+genus_info$number_NRPS
genus_info$genus <- factor(genus_info$Genus)
genus_colors <- c('Bacillus' = '#66cc8f', 
                  'Corynebacterium' = '#4b545e', 
                  'Rhodococcus' = '#f25641', 
                  'Staphylococcus' = '#ffb300', 
                  'Streptomyces' = '#66b3d9')



library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(treeio)
library(dplyr)

p_base <- ggtree(tree, layout = 'fan', open.angle = 15, branch.length = "none") %<+% genus_info

p_base
p_genus <- p_base +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = strain_name, fill = genus),
    width = 5,
    offset = 0.1,
    color = NA
  ) +
  scale_fill_manual(values = genus_colors, name = "Genus")
p_genus
# add number_NIS ring
p_NIS <- p_genus +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = strain_name, x = number_NIS),
    stat = "identity",
    orientation = "y",
    width = 0.3,
    offset = 0.07,
    fill = "#b71e36"
  )
p_NIS
# add number_NRPS ring
p_NRPS <- p_NIS +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = strain_name, x = number_NRPS),
    stat = "identity",
    orientation = "y",
    width = 0.3,
    offset = 0.02,
    fill = "#084b4a"
  )
p_NRPS


