data <- read.csv('../input/bgc_pfam_matrix.csv')
species_info <- read.csv('../input/NCBI_all_gtdb_meta_info.csv')
gram_pos_strain <- read.csv("../input/genome_bama_screening.csv", stringsAsFactors = FALSE)
gram_pos_strain <- gram_pos_strain[gram_pos_strain$Has_BamA =='False',]
data_pos <- data[data$Subfolder %in% gram_pos_strain$Strain_Name,]
A = unique(data_pos$Subfolder)
SBP_domain <- c("SBP_bac_5", "SBP_bac_3", "Peripla_BP_2", "ZnuA", "SBP_bac_1", "Bmp", "Lipoprotein_9", "DctP", "OpuAC", "PBP_like_2", "Phosphonate_bd", "SBP_bac_6", "Peripla_BP_4", "SBP_bac_8", "Peripla_BP_6", "SBP_bac_11", "NMT1_3", "NMT1", "NosD", "ABC_sub_bind", "TctC", "PBP_like", "Peripla_BP_1", "Cypl", "Peripla_BP_7", "Peripla_BP_5", "LppC", "DUF3798", "DUF3834", "SBP_bac_10", "ANF_receptor", "ABC_transp_aux", "Peripla_BP_3")

valid_cols <- intersect(colnames(data_pos), c(SBP_domain,''))

domain_counts <- sapply(data_pos[, valid_cols], function(x) sum(x > 0, na.rm = TRUE))

df_domain_plot <- data.frame(
  Domain = names(domain_counts),
  Frequency = as.numeric(domain_counts)
)

df_domain_plot <- df_domain_plot[df_domain_plot$Frequency > 0, ]
df_domain_plot$occur <- df_domain_plot$Frequency/nrow(data_pos)

library(ggplot2)
library(dplyr)

df_domain_plot <- df_domain_plot %>% arrange(desc(Frequency))

p_domain <- ggplot(df_domain_plot, aes(x = reorder(Domain, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "#4E84C4", width = 0.7, color = "black") +
  geom_text(aes(label = Frequency), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Frequency of SBP Domains in BGCs",
       x = "SBP Domain",
       y = "Frequency (Count)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

print(p_domain)

ggsave("../Figure/SBP_Domain_Frequency.pdf", p_domain, width = 8, height = max(6, 0.3 * nrow(df_domain_plot)))