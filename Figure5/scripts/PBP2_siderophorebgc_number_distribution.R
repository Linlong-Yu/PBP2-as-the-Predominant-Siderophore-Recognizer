library(dplyr)
library(ggplot2)
library(patchwork)
library(ggsci) 
library(scales)
library(tidyr)

strain_overview <- readRDS("strain_sider_PBP2_overiw.rds")
bgc_data <- read.csv("../input/bgc_pbp2_summary_typed.csv", stringsAsFactors = FALSE)
species_info <- read.csv('../input/NCBI_all_gtdb_meta_info.csv', stringsAsFactors = FALSE)

species_info_clean <- species_info %>%
  distinct(Strain_name, .keep_all = TRUE) %>%
  select(Strain_name, Phylum, Class, Order)

strain_full <- strain_overview %>%
  left_join(species_info_clean, by = c("Strain_Name" = "Strain_name")) %>%
  mutate(Order = gsub("^o__", "", Order))

order_stats <- strain_full %>%
  group_by(Order) %>%
  summarise(Count = n()) %>%
  mutate(Percent = Count / sum(Count)) %>%
  arrange(desc(Percent))

major_orders <- order_stats %>%
  filter(Percent >= 0.005) %>% 
  pull(Order)

strain_plot_data <- strain_full %>%
  mutate(
    Taxonomy_Show = ifelse(Order %in% major_orders & !is.na(Order) & Order != "", Order, "Other")
  )
taxa_levels <- c(major_orders, "Other")
strain_plot_data$Taxonomy_Show <- factor(strain_plot_data$Taxonomy_Show, levels = taxa_levels)

pbp2_total  <- sum(strain_overview$Genome_Total_PBP2, na.rm = TRUE)
pbp2_inside <- sum(bgc_data$PBP2_Count, na.rm = TRUE)
pbp2_outside <- pbp2_total - pbp2_inside
pbp2_pie_data <- data.frame(Category = c("Inside BGC", "Scattered"), Count = c(pbp2_inside, pbp2_outside)) %>% 
  mutate(Fraction = Count/sum(Count), Label = paste0(round(Fraction*100,1),"%"))

total_nrps <- sum(strain_overview$NRPS_Count, na.rm = TRUE)
total_nis  <- sum(strain_overview$NIS_Count, na.rm = TRUE)
sidero_pie_data <- data.frame(Category = c("NRPS", "NIS"), Count = c(total_nrps, total_nis)) %>% 
  mutate(Fraction = Count/sum(Count), Label = paste0(round(Fraction*100,1),"%"))

theme_ai_optimized <- theme_classic(base_size = 14) +
  theme(
   
    text = element_text(family = "sans", color = "black", face = "plain"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "plain"),
    plot.title = element_blank(), 
  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    
    legend.position = "top",
    legend.justification = "center"
  )

p_pbp2_main <- ggplot(strain_plot_data, aes(x = Genome_Total_PBP2, fill = Taxonomy_Show)) +
  geom_histogram(binwidth = 1, color = "black", linewidth = 0.1, alpha = 0.9) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_d3("category20", name = "Order") +
  labs(x = "Number of PBP2 Genes per Genome", y = "Number of Strains") +
  theme_ai_optimized +
  
  guides(fill = guide_legend(nrow = 3, title.position = "top", title.hjust = 0.5))


p_pie_1 <- ggplot(pbp2_pie_data, aes(x = "", y = Fraction, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Inside BGC" = "#D62728", "Scattered" = "#7F7F7F")) + 
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  theme_void() + labs(title = "Location") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "bottom", legend.title = element_blank())

final_pbp2 <- p_pbp2_main + 
  inset_element(p_pie_1, left = 0.5, bottom = 0.4, right = 0.98, top = 0.95)

max_sidero <- max(strain_plot_data$Total_Siderophores, na.rm=TRUE)

p_sidero_main <- ggplot(strain_plot_data, aes(x = Total_Siderophores, fill = Taxonomy_Show)) +
  geom_bar(width = 0.7, color = "black", linewidth = 0.1, alpha = 0.9) +
  scale_x_continuous(breaks = seq(0, max_sidero, 1), limits = c(-0.5, max_sidero + 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_d3("category20", name = "Order") +
  labs(x = "Number of Siderophore BGCs per Genome", y = "Number of Strains") +
  theme_ai_optimized +
  
  guides(fill = guide_legend(nrow = 3, title.position = "top", title.hjust = 0.5))

p_pie_2 <- ggplot(sidero_pie_data, aes(x = "", y = Fraction, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("NRPS" = "#2CA02C", "NIS" = "#98DF8A")) + 
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  theme_void() + labs(title = "Type") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "bottom", legend.title = element_blank())

final_sidero <- p_sidero_main + 
  inset_element(p_pie_2, left = 0.5, bottom = 0.4, right = 0.98, top = 0.95)

combined_plot <- final_pbp2 + final_sidero + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top")    
combined_plot

if(!dir.exists("../Figure")) dir.create("../Figure", recursive = TRUE)
#ggsave("../Figure/Fig_Distribution_TopLegend.pdf", combined_plot, width = 11, height = 5)
