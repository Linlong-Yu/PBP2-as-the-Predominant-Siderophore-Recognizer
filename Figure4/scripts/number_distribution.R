library(stringr)
data <- read.csv('../input/rec_in_5genera.csv')
sider_info <- read.csv('../input/sider_info.csv')

library(dplyr)
rec_info <- data %>%
  count(Strain_Name, name = "rec_number")
sider_count <- sider_info %>%
  count(strainName, name = "sider_number")
final_info <- rec_info %>%
  rename(strainName = Strain_Name) %>%
  left_join(sider_count, by = "strainName")

final_info$sider_number[is.na(final_info$sider_number)] <- 0

genus_info <- data %>%
  select(Strain_Name, Genus) %>%
  distinct() 

final_info <- final_info %>%
  left_join(genus_info, by = c("strainName" = "Strain_Name"))

library(ggplot2)
library(dplyr)
final_info <- final_info %>%
  mutate(Genus = factor(Genus, levels = c("Streptomyces", "Corynebacterium", "Rhodococcus", "Bacillus", "Staphylococcus")))
library(ggplot2)
library(ggplot2)
library(scales)
ggplot(final_info, aes(x = rec_number)) +
  geom_histogram(bins = 10, fill = "skyblue", alpha = 0.7) +
  facet_wrap(~ Genus, scales = "free", nrow = 1) +
  labs(y = "Strain number") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),                 
    axis.line = element_line(color = "black"),    
    axis.ticks = element_line(color = "black"),   
    axis.ticks.length = unit(-0.1, "cm"),  
    axis.text.x  = element_text(color = "black"),
    axis.text.y  = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )





