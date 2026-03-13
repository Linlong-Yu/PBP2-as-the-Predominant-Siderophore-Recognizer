library(dplyr)
library(tidyr)      
library(ggplot2)
library(ggrepel)
library(ggExtra)   
library(patchwork)  
library(grid)      

data <- read.csv('../input/bgc_pfam_matrix.csv')

gram_pos_strain <- read.csv("../input/genome_bama_screening.csv", stringsAsFactors = FALSE)
gram_pos_strain <- gram_pos_strain[gram_pos_strain$Has_BamA == 'False', ]

all_gram_pos_strains <- data.frame(Strain_Name = unique(gram_pos_strain$Strain_Name))

data_pos <- data[data$Subfolder %in% gram_pos_strain$Strain_Name, ]

PBP2_data <- read.csv("//192.168.1.150/share/NAS_923/yll/PBP2_analysis/pbp2_analysis_results/genome_pbp2_summary.csv", stringsAsFactors = FALSE)

species_info <- read.csv('../../../gtdb_meta_info/NCBI_all_gtdb_meta_info.csv')

species_info <- species_info[!duplicated(species_info$Strain_name),]
colnames(species_info)[2] = 'Strain_Name' 

rm(data)

bgc_stats <- data_pos %>%
  group_by(Subfolder) %>%
  summarise(
    Total_Siderophores = n(), 
    NRPS_Count = sum(Category == "NRP-metallophore", na.rm = TRUE),
    NIS_Count  = sum(Category != "NRP-metallophore", na.rm = TRUE),
    Siderophores_with_PBP2 = sum(Peripla_BP_2 > 0, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  rename(Strain_Name = Subfolder)

pbp2_gene_stats <- PBP2_data %>%
  group_by(Strain_Name) %>%
  summarise(Genome_Total_PBP2 = n()) %>%
  ungroup()

strain_overview <- all_gram_pos_strains %>%

  left_join(select(species_info, Strain_Name, Genus), by = "Strain_Name") %>%
  
  left_join(bgc_stats, by = "Strain_Name") %>%
  
  left_join(pbp2_gene_stats, by = "Strain_Name") %>%
 
  mutate(
    Total_Siderophores = replace_na(Total_Siderophores, 0),
    NRPS_Count = replace_na(NRPS_Count, 0),
    NIS_Count = replace_na(NIS_Count, 0),
    Siderophores_with_PBP2 = replace_na(Siderophores_with_PBP2, 0),
    Genome_Total_PBP2 = replace_na(Genome_Total_PBP2, 0)
  )
saveRDS(strain_overview,'strain_sider_PBP2_overiw.rds')
strain_overview <-readRDS('strain_sider_PBP2_overiw.rds')

MIN_GENUS_SIZE <- 20 

genus_overview <- strain_overview %>%
  filter(!is.na(Genus) & Genus != "") %>%
  group_by(Genus) %>%
  summarise(
    Genus_Size = n(), 
    
    Avg_Total_Siderophores = mean(Total_Siderophores, na.rm = TRUE),
    Avg_Genome_PBP2        = mean(Genome_Total_PBP2, na.rm = TRUE),
    Avg_Sidero_with_PBP2   = mean(Siderophores_with_PBP2, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  filter(Genus_Size >= MIN_GENUS_SIZE) %>% 
  mutate(
    
    
    Denominator = ifelse(Avg_Total_Siderophores == 0, 1, Avg_Total_Siderophores), # 防止除以0
    Ratio_PBP2_In_Sidero = ifelse(Avg_Total_Siderophores == 0, 0, Avg_Sidero_with_PBP2 / Denominator),
    Ratio_GenomePBP2_vs_Sidero = ifelse(Avg_Total_Siderophores == 0, 0, Avg_Genome_PBP2 / Denominator),
    Ratio_cheat_PBP2 = Avg_Genome_PBP2-Avg_Total_Siderophores,
    
    Genus_Clean = gsub("^g__", "", Genus)
  )


my_clean_theme <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

genus_overview_expanded <- genus_overview %>% 
  uncount(Genus_Size, .remove = FALSE) 

bw_x <- stats::bw.nrd0(na.omit(genus_overview_expanded$Avg_Genome_PBP2))
bw_y <- stats::bw.nrd0(na.omit(genus_overview_expanded$Avg_Total_Siderophores))

top_sidero <- genus_overview %>% arrange(desc(Avg_Total_Siderophores)) %>% slice(1:5) %>% pull(Genus)
top_pbp2   <- genus_overview %>% arrange(desc(Avg_Genome_PBP2)) %>% slice(1:5) %>% pull(Genus)
top_ratio  <- genus_overview %>% arrange(desc(Ratio_PBP2_In_Sidero)) %>% slice(1:5) %>% pull(Genus) 
target_labels <- unique(c(top_sidero, top_pbp2, top_ratio))
genus_overview$Label_Show <- ifelse(genus_overview$Genus %in% target_labels, 
                                    genus_overview$Genus_Clean, NA)

p_scatter_main <- ggplot(genus_overview_expanded, aes(x = Avg_Genome_PBP2, y = Avg_Total_Siderophores)) +
  
 
  geom_point(data = genus_overview, 
             aes(size = Genus_Size, color = Ratio_PBP2_In_Sidero), 
             alpha = 0.8) +
  
  geom_text_repel(data = genus_overview,
                  aes(label = Label_Show), 
                  size = 3.5, fontface = "italic", 
                  max.overlaps = Inf, min.segment.length = 0) +
  
  scale_color_gradient(low = "grey80", high = "#b30000", name = "PBP2% inside Sidero") +
  
  scale_size(
    range = c(2, 10), 
    name = "Genus Count",
    breaks = c(20,100,500,2000) 
  ) +

scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  
  labs(x = "Average Genome PBP2 Genes", 
       y = "Average Total Siderophores") +
  
  my_clean_theme +
  theme(legend.position = "bottom", legend.box = "horizontal")


p_scatter_final <- ggMarginal(
  p_scatter_main, 
  type = "density", 
  margins = "both", 
  color = "black", 
  alpha = 0.6,
  size = 4,
  
  xparams = list(fill = "#377EB8", bw = bw_x), 
  yparams = list(fill = "#4DAF4A", bw = bw_y)  
)
p_scatter_final

ggsave("../Figure/Genus_Scatter_Marginal_Weighted.pdf", p_scatter_final, width = 8, height = 7)

global_limits <- range(genus_overview$Genus_Size) 


my_scale_size <- scale_size(
  range = c(2, 10),                    
  name = "Genus Count",
  limits = global_limits,               
  breaks = c(20, 100, 500, 2000)       
)


df_top10_pbp2 <- genus_overview %>% 
  arrange(desc(Avg_Genome_PBP2)) %>% 
  slice(1:10)

p_bubble_pbp2 <- ggplot(df_top10_pbp2, aes(x = Avg_Genome_PBP2, y = reorder(Genus_Clean, Avg_Genome_PBP2))) +
  geom_segment(aes(yend = Genus_Clean, xend = 0), color = "grey90", linetype = "dashed") +
  geom_point(aes(size = Genus_Size), color = "#377EB8", alpha = 0.8) +
  geom_text(aes(label = round(Avg_Genome_PBP2, 2)), vjust = -0.8, size = 3.5) +
  
  my_scale_size + 
  
  labs(x = "Average Genome PBP2 Count", y = "") +
  my_clean_theme +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right")

df_top10_sidero <- genus_overview %>% 
  arrange(desc(Avg_Total_Siderophores)) %>% 
  slice(1:10)

p_bubble_sidero <- ggplot(df_top10_sidero, aes(x = Avg_Total_Siderophores, y = reorder(Genus_Clean, Avg_Total_Siderophores))) +
  geom_segment(aes(yend = Genus_Clean, xend = 0), color = "grey90", linetype = "dashed") +
  geom_point(aes(size = Genus_Size), color = "#4DAF4A", alpha = 0.8) +
  geom_text(aes(label = round(Avg_Total_Siderophores, 2)), vjust = -0.8, size = 3.5) +
  
  my_scale_size + 
  
  labs(x = "Average Total Siderophores", y = "") +
  my_clean_theme +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right")

genus_ratio_filtered <- genus_overview %>% filter(Avg_Total_Siderophores > 0.5)

# --- Top 10 ---
df_ratio_top <- genus_ratio_filtered %>% arrange(desc(Ratio_GenomePBP2_vs_Sidero)) %>% slice(1:10)

p_bubble_ratio_top <- ggplot(df_ratio_top, aes(x = Ratio_GenomePBP2_vs_Sidero, y = reorder(Genus_Clean, Ratio_GenomePBP2_vs_Sidero))) +
  geom_segment(aes(yend = Genus_Clean, xend = 0), color = "grey90", linetype = "dashed") +
  geom_point(aes(size = Genus_Size), color = "#E41A1C", alpha = 0.8) +
  geom_text(aes(label = round(Ratio_GenomePBP2_vs_Sidero, 2)), vjust = -0.8, size = 3.5) +
  
  my_scale_size + 
  
  labs(x = "Ratio (Genome PBP2 / Siderophores)", y = "", subtitle = "Top 10 High Ratio (> 0.5 Siderophores)") +
  my_clean_theme +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right")

# --- Bottom 10 ---
df_ratio_bot <- genus_ratio_filtered %>% arrange(Ratio_GenomePBP2_vs_Sidero) %>% slice(1:10)

p_bubble_ratio_bot <- ggplot(df_ratio_bot, aes(x = Ratio_GenomePBP2_vs_Sidero, y = reorder(Genus_Clean, -Ratio_GenomePBP2_vs_Sidero))) +
  geom_segment(aes(yend = Genus_Clean, xend = 0), color = "grey90", linetype = "dashed") +
  geom_point(aes(size = Genus_Size), color = "grey50", alpha = 0.8) +
  geom_text(aes(label = round(Ratio_GenomePBP2_vs_Sidero, 2)), vjust = -0.8, size = 3.5) +
  
  my_scale_size + 
  
  labs(x = "Ratio (Genome PBP2 / Siderophores)", y = "", subtitle = "Bottom 10 Low Ratio (> 0.5 Siderophores)") +
  my_clean_theme +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right")

genus_cheat_analysis <- genus_overview %>%
  filter(Avg_Genome_PBP2 > 1) %>% 
  mutate(Difference_PBP2_minus_Sidero = Avg_Genome_PBP2 - Avg_Total_Siderophores)

# --- Top 10 ---
df_diff_top <- genus_cheat_analysis %>% arrange(desc(Difference_PBP2_minus_Sidero)) %>% slice(1:10)

p_bubble_diff_top <- ggplot(df_diff_top, aes(x = Difference_PBP2_minus_Sidero, y = reorder(Genus_Clean, Difference_PBP2_minus_Sidero))) +
  geom_segment(aes(yend = Genus_Clean, xend = 0), color = "grey90", linetype = "dashed") +
  geom_point(aes(size = Genus_Size), color = "#984EA3", alpha = 0.8) + 
  geom_text(aes(label = round(Difference_PBP2_minus_Sidero, 2)), vjust = -0.8, size = 3.5) +
  
  my_scale_size + 
  
  labs(x = "Difference", y = "", subtitle = "Top 10: Excess Receptors (Avg PBP2 > 1)") +
  my_clean_theme +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right")

# --- Bottom 10 ---
df_diff_bot <- genus_cheat_analysis %>% arrange(Difference_PBP2_minus_Sidero) %>% slice(1:10)

p_bubble_diff_bot <- ggplot(df_diff_bot, aes(x = Difference_PBP2_minus_Sidero, y = reorder(Genus_Clean, -Difference_PBP2_minus_Sidero))) +
  geom_segment(aes(yend = Genus_Clean, xend = 0), color = "grey90", linetype = "dashed") +
  geom_point(aes(size = Genus_Size), color = "#FF7F00", alpha = 0.8) +
  geom_text(aes(label = round(Difference_PBP2_minus_Sidero, 2)), vjust = -0.8, size = 3.5) +
  
  my_scale_size + 
  
  labs(x = "Difference", y = "", subtitle = "Bottom 10: Excess Production (Avg PBP2 > 1)") +
  my_clean_theme +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right")

if(!dir.exists("../Figure")) dir.create("../Figure", recursive = TRUE)

ggsave("../Figure/Bubble_Rank_Top10_GenomePBP2.pdf", p_bubble_pbp2, width = 7, height = 5)
ggsave("../Figure/Bubble_Rank_Top10_Siderophores.pdf", p_bubble_sidero, width = 7, height = 5)
ggsave("../Figure/Bubble_Rank_Ratio_Top10.pdf", p_bubble_ratio_top, width = 7, height = 5)
ggsave("../Figure/Bubble_Rank_Ratio_Bottom10.pdf", p_bubble_ratio_bot, width = 7, height = 5)
ggsave("../Figure/Bubble_Diff_Top10_Cheaters.pdf", p_bubble_diff_top, width = 7, height = 5)
ggsave("../Figure/Bubble_Diff_Bottom10_Producers.pdf", p_bubble_diff_bot, width = 7, height = 5)

cat("\n\n>>> 1. Top 10 Genera by Avg Genome PBP2 Count <<<\n")
print(
  df_top10_pbp2 %>% 
    select(Genus_Clean, Avg_Genome_PBP2, Genus_Size) %>% 
    as.data.frame() # 
)

cat("\n\n>>> 2. Top 10 Genera by Avg Total Siderophores <<<\n")
print(
  df_top10_sidero %>% 
    select(Genus_Clean, Avg_Total_Siderophores, Genus_Size) %>% 
    as.data.frame()
)

cat("\n\n>>> 3. Top 10 High Ratio (Genome PBP2 / Siderophores) <<<\n")
print(
  df_ratio_top %>% 
    select(Genus_Clean, Ratio_GenomePBP2_vs_Sidero, Genus_Size) %>% 
    as.data.frame()
)

cat("\n\n>>> 4. Bottom 10 Low Ratio (Genome PBP2 / Siderophores) <<<\n")
print(
  df_ratio_bot %>% 
    select(Genus_Clean, Ratio_GenomePBP2_vs_Sidero, Genus_Size) %>% 
    as.data.frame()
)

cat("\n\n>>> 5. Top 10 Difference (Excess Receptors / Potential Cheaters) <<<\n")
print(
  df_diff_top %>% 
    select(Genus_Clean, Difference_PBP2_minus_Sidero, Genus_Size) %>% 
    as.data.frame()
)

cat("\n\n>>> 6. Bottom 10 Difference (Excess Production / Strong Producers) <<<\n")
print(
  df_diff_bot %>% 
    select(Genus_Clean, Difference_PBP2_minus_Sidero, Genus_Size) %>% 
    as.data.frame()
)