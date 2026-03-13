library(dplyr)
library(ggplot2)
library(ggsci)    
library(rstatix)  
library(ggpubr)    
library(scales)    

bgc_data <- read.csv("../input/bgc_pbp2_summary_typed.csv", stringsAsFactors = FALSE)

plot_data <- bgc_data %>%
  mutate(
    Category = case_when(
      grepl("NRP-metallophore", BGC_Type) ~ "NRPS-siderophore",
      grepl("NI-siderophore", BGC_Type)   ~ "NIS-siderophore",
      TRUE                                ~ "Other BGC"
    ),
    Has_PBP2 = ifelse(PBP2_Count > 0, "Present", "Absent")
  )

stats_summary <- plot_data %>%
  group_by(Category) %>%
  summarise(
    Total_Count = n(),
    PBP2_Count  = sum(Has_PBP2 == "Present"),
    Ratio       = PBP2_Count / Total_Count
  ) %>%
  mutate(Category = factor(Category, levels = c("NRPS-siderophore", "NIS-siderophore", "Other BGC")))

print(stats_summary)

contingency_table_2d <- table(plot_data$Category, plot_data$Has_PBP2)

pairwise_res <- pairwise_fisher_test(contingency_table_2d, p.adjust.method = "bonferroni") %>%
  filter(!((group1 == "NRPS-siderophore" & group2 == "NIS-siderophore") | 
             (group1 == "NIS-siderophore" & group2 == "NRPS-siderophore"))) %>%
  mutate(
    label_final = case_when(
      p.adj == 0 ~ "p < 2.2e-16",
      p.adj < 0.001 ~ paste0("p = ", formatC(p.adj, format = "e", digits = 2)),
      TRUE ~ paste0("p = ", sprintf("%.3f", p.adj))
    )
  )
pairwise_res_test <- pairwise_res %>%

  filter(!((group1 == "NRPS-siderophore" & group2 == "NIS-siderophore") | 
             (group1 == "NIS-siderophore" & group2 == "NRPS-siderophore"))) %>%
  mutate(
    
    label_final = paste0("p = ", formatC(p.adj, format = "e", digits = 2))
  )
max_ratio <- max(stats_summary$Ratio)
pairwise_res <- pairwise_res %>%
  mutate(y.position = seq(from = max_ratio * 1.05, by = 0.08, length.out = n()))

theme_ai_fig3 <- theme_classic(base_size = 14) +
  theme(
    
    text = element_text(family = "sans", color = "black", face = "plain"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "plain"),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none" 
  )

p_fig3 <- ggplot(stats_summary, aes(x = Category, y = Ratio, fill = Category)) +

  geom_bar(stat = "identity", color = "black", width = 0.6, linewidth = 0.5, alpha = 0.9) +
  
  geom_text(aes(label = scales::percent(Ratio, accuracy = 0.1)), 
            vjust = -0.5, size = 4, fontface = "plain") +
  
  stat_pvalue_manual(
    pairwise_res, 
    label = "label_final", 
    tip.length = 0.02,
    size = 4,
    bracket.size = 0.5 
  ) +
  
  scale_y_continuous(
    labels = scales::percent, 
    limits = c(0, max(pairwise_res$y.position) * 1.15),
    expand = c(0, 0)
  ) +
  
  scale_fill_npg() + 
  
  labs(x = "", y = "Percentage of BGCs containing PBP2") +
  
  theme_ai_fig3


if(!dir.exists("../Figure")) dir.create("../Figure", recursive = TRUE)
ggsave("../Figure/Figure3_BGC_PBP2_Ratio_AI_Ready.pdf", p_fig3, width = 6, height = 7)

print(p_fig3)
