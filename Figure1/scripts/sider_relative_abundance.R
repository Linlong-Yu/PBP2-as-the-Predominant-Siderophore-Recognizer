data <- read.csv('../input/syn_info_5_genera.csv')
library(dplyr)
library(ggplot2)
data <- data %>%
  group_by(Genus) %>%
  mutate(total_size = sum(group_size),  
         relative_abundance = group_size / total_size)  
#calculate α-diversity
shannon_index <- data %>%
  group_by(Genus) %>%
  summarise(H = -sum(relative_abundance * log(relative_abundance), na.rm = TRUE))  # 计算香农指数
ggplot(shannon_index, aes(x = Genus, y = 1, fill = H)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient(low = "#81c3d7", high = "#f25e7a") +
  labs(x = "Index", y = "", fill = "Value") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())


library(ggrepel)
# set seed
set.seed(123)
# create fixed position
pos <- position_jitter(width = 0.2, height = 0, seed = 123)
data$Genus <- factor(data$Genus, levels = c("Bacillus", "Corynebacterium", "Rhodococcus", "Staphylococcus", "Streptomyces"))
p <- ggplot(data, aes(x = Genus, y = relative_abundance, color = Genus)) +
  geom_point(
    aes(group = Genus),
    position = pos, 
    alpha = 0.7, 
    size = 3
  ) +
  geom_text_repel(
    aes(label = label),
    position = pos,
    size = 3, 
    color = "black", 
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "black",
    segment.size = 0.3,
    force = 10,
    force_pull = 1,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      'Bacillus' = '#66cc8f', 
      'Corynebacterium' = '#4b545e', 
      'Rhodococcus' = '#f25641', 
      'Staphylococcus' = '#ffb300', 
      'Streptomyces' = '#66b3d9'
    ),
    
  ) +
  theme_minimal() +
  labs(x = "Genus", y = "Relative Abundance") +
  theme(axis.text.x = element_text(size = 11))

print(p)
