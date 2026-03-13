library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
final_plot_data <- readRDS('../input/PBP2_number_info.rds')
final_plot_data$Genus <- factor(final_plot_data$Genus, 
                                levels = c("Streptomyces", "Corynebacterium", 
                                           "Rhodococcus", "Bacillus", "Staphylococcus"))

p <- ggplot(final_plot_data, aes(x = pbp2_number)) +
 
  geom_histogram(binwidth = 1, fill = "skyblue", color = "white", alpha = 0.8) +
 
  facet_wrap(~ Genus, nrow = 1) + 
  
  labs(y = "Strain number", x = "Copy number of PBP2 genes") +

  scale_x_continuous(
    breaks = pretty_breaks(),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = pretty_breaks() 
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),                  
    axis.line = element_line(color = "black"),     
    axis.ticks = element_line(color = "black"),    
    axis.ticks.length = unit(-0.1, "cm"), 
    axis.text.x  = element_text(color = "black", margin = margin(t = 5)), # 调整文字间距防止重叠
    axis.text.y  = element_text(color = "black", margin = margin(r = 5)),
    strip.text = element_text(face = "bold", size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


print(p)

ggsave('../Figure/pbp2_distribution_standardized.pdf', plot = p, width = 15, height = 7)
plot_data_normalized <- final_plot_data %>%
  
  group_by(Genus, pbp2_number) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  
  mutate(proportion = count / sum(count)) %>%
  ungroup()

plot_data_normalized$Genus <- factor(plot_data_normalized$Genus, 
                                     levels = c("Streptomyces", "Corynebacterium", 
                                                "Rhodococcus", "Bacillus", "Staphylococcus")) 
 
library(scales)

p_norm <- ggplot(plot_data_normalized, aes(x = pbp2_number, y = proportion)) +
  
  geom_bar(stat = "identity", fill = "skyblue", color = "white", alpha = 0.8, width = 0.8) +
  
  geom_text(
    aes(label = scales::percent(proportion, accuracy = 1)), 
    vjust = -0.5,  
    size = 3.5,    
    color = "black"
  ) +
  
  facet_wrap(~ Genus, nrow = 1) + 
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
   
    expand = expansion(mult = c(0, 0.15)) 
  ) +
  
  scale_x_continuous(
    breaks = pretty_breaks(),
    expand = c(0, 0)
  ) +
  
  labs(
    y = "Percentage of Strains", 
    x = "PBP2 Gene Count",
    title = "Distribution of PBP2 Gene Counts by Genus (Normalized)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),      
    axis.ticks = element_line(color = "black"),     
    axis.ticks.length = unit(-0.1, "cm"), 
    axis.text.x  = element_text(color = "black", margin = margin(t = 5)),
    axis.text.y  = element_text(color = "black", margin = margin(r = 5)),
    strip.text = element_text(face = "bold", size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3) # 注意: size在新版ggplot2中建议用linewidth
  )


print(p_norm)
ggsave('../Figure/pbp2_distribution_standardized_percentage.pdf', plot = p_norm, width = 15, height = 5)
data <- read.csv('../input/NCBI_all_gtdb_meta_info.csv')
data_Rhodococcus <- data[data$Species == 's__Rhodococcus erythropolis' | data$Species =='s__Rhodococcus qingshengii',]
final_plot_high_PBP2 <- final_plot_data[final_plot_data$Strain_Name %in% data_Rhodococcus$Strain_name,]
mean(final_plot_high_PBP2$pbp2_number)
test <- final_plot_data[final_plot_data$pbp2_number >= 4 & final_plot_data$pbp2_number <= 13,]
