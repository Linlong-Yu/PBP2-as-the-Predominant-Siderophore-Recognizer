data = read.csv('../output/128_resulte_of_correlation_and_counts_one2many.csv')
data$proportion = data$num_genes/128
data = data[-c(457,492,40,580,145,500),]#remove synthetase
library(ggplot2)
#highlight SBP
highlight_domains <- c("NMT1", "NMT1_3", "OpuAC", "PBP_like_2", "Peripla_BP_1", "Peripla_BP_2", 
                       "Peripla_BP_3", "Peripla_BP_4", "Peripla_BP_5", "Peripla_BP_6", 
                       "SBP_bac_1", "SBP_bac_11", "SBP_bac_3", "SBP_bac_5", "SBP_bac_6", 
                       "SBP_bac_8", "TctC", "ZnuA")

library(ggrepel)
data$highlight <- ifelse(data$aSDomain %in% highlight_domains, "SBP subtype gene", "Other gene")
library(grid)
data_transport <- subset(data, highlight == "SBP subtype gene")
data_nontransport <- subset(data, highlight != "SBP subtype gene")

p <- ggplot() +
 
  geom_point(data = data_nontransport,
             aes(x = correlation, y = proportion, color = highlight),
             alpha = 0.6, size = 5) +
  
  geom_point(data = data_transport,
             aes(x = correlation, y = proportion, color = highlight),
             alpha = 0.6, size = 5) +
  
  geom_text_repel(data = data_transport,
                  aes(x = correlation, y = proportion, label = aSDomain),
                  size = 5, color = "black",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = 10) +
  
  scale_color_manual(values = c("SBP subtype gene" = "#ff6666", 
                                "Other gene" = "#81d4fa")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_minimal() +
  labs(
    x = "Correlation",
    y = "Proportion",
    color = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(5, "pt")
  )

print(p)

