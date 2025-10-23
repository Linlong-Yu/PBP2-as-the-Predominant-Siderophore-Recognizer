library(dplyr)
data = read.csv('../input/summ_genus.csv',row.names = 1)
data <- data %>%
  mutate(
    NIS = ifelse(grepl("NI-siderophore", syn_type), 1, 0),
    NRPS = ifelse(grepl("NRP-metallophore", syn_type), 1, 0)
  )
data_summary <- data %>%
  group_by(Genus) %>%
  summarise(
    NIS_count = sum(NIS),
    NRPS_count = sum(NRPS),
    Total = n()
  )
SM_info = read.csv("../input/genus_gbk_summary.csv")
SM_info$Genus[3] = 'Bacillus'
SM_info_ordered <- SM_info[match(data_summary$Genus, SM_info$Genus), ]
SecondMeta_sider_info <- cbind(data_summary, SM_info_ordered)
SecondMeta_sider_info$Total_secondmeta = SecondMeta_sider_info$Total_gbk_Files - SecondMeta_sider_info$Subfolder_Count
SecondMeta_sider_info$average_sm = SecondMeta_sider_info$Total_secondmeta/SecondMeta_sider_info$Subfolder_Count
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
SecondMeta_sider_info = SecondMeta_sider_info[,-5]
data_long <- SecondMeta_sider_info %>%
  
  mutate(Other_BGC = Total_secondmeta - NIS_count - NRPS_count) %>%
  
  select(Genus, NIS_count, NRPS_count, Other_BGC, Subfolder_Count) %>%
 
  pivot_longer(cols = c(NIS_count, NRPS_count, Other_BGC),
               names_to = "BGC_Type",
               values_to = "Count") %>%

  group_by(Genus) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(Genus = fct_reorder(Genus, Total))

# add labels
data_long$BGC_Type <- factor(data_long$BGC_Type,
                             levels = c("Other_BGC", "NIS_count", "NRPS_count"),
                             labels = c("Other BGCs", "NIS-siderophore BGC", "NRPS-siderophore BGC"))
data_long$average_counts = data_long$Count/data_long$Subfolder_Count
# plot
p <- ggplot(data_long, aes(x = average_counts, y = Genus, fill = BGC_Type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "Other BGCs" = "grey70",
    "NIS-siderophore BGC" = "#d52a4b",
    "NRPS-siderophore BGC" = "#248c9a"
  )) +
  scale_x_continuous(expand = c(0, 0)) +  
  theme_minimal() +
  labs(
    title = "BGC Distribution by Genus",
    x = "Count",
    y = "",
    fill = "BGC Type"
  ) +
  theme(
    panel.grid = element_blank(),                        
    axis.line = element_line(color = "black"),          
    axis.ticks = element_line(color = "black"),          
    panel.border = element_rect(color = "black", fill = NA, size = 0.8), 
    legend.position = "bottom"
  )

p

