data = read.csv('../input/summ_genus.csv',row.names = 1)

library(ggplot2)
library(dplyr)
# create groups
data <- data %>%
  mutate(Known_Category = case_when(
    Similarity >= 0.8 ~ "High confidence (>80%)",
    Similarity >= 0.5 ~ "Medium confidence (50%-80%)",
    Similarity >= 0.3 ~ "Low confidence (30%-50%)",
    TRUE ~ "Unknown(<30%)"
  ))

# calculate proportion in different genera
data_summary <- data %>%
  group_by(Genus, Known_Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Genus) %>%
  mutate(Proportion = Count / sum(Count))

data_summary$Known_Category <- factor(data_summary$Known_Category, 
                                      levels = c("High confidence (>80%)", 
                                                 "Medium confidence (50%-80%)", 
                                                 "Low confidence (30%-50%)", 
                                                 "Unknown(<30%)"))

genus_labels <- c(
  
  Bacillus = expression(italic("Bacillus")),
  Corynebacterium = expression(italic("Corynebacterium")),
  Rhodococcus = expression(italic("Rhodococcus")),
  Staphylococcus = expression(italic("Staphylococcus")),
  Streptomyces = expression(italic("Streptomyces"))
)

p <- ggplot(data_summary, aes(x = Genus, y = Proportion, fill = Known_Category)) +
  geom_col(width = 0.7, position = "stack") +
  scale_fill_manual(values = c(
    "High confidence (>80%)" = "#007474",
    "Medium confidence (50%-80%)" = "#26a6a6", 
    "Low confidence (30%-50%)" = "#80c2c2", 
    "Unknown(<30%)" = "#b3b3b3")) +
  scale_x_discrete(labels = genus_labels) +
  scale_y_continuous(expand = c(0, 0)) +   
  labs(x = "Genus", 
       y = "Number of siderophore BGCs", 
       fill = "Siderophore BGC similarity to MIBiG") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, face = "plain"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )
print(p)

