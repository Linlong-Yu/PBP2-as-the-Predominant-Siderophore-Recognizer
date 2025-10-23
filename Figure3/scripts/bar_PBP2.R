setwd('E:/Project/2024/Gram_pos_updata/Bacillus_new/coding/Feature_sequence')
library(ggplot2)

total_length <- 305

bar_data <- data.frame(
  pos = 1:total_length,
  fill = "black"  
)

bar_data$fill[1:18] <- "lightblue"
bar_data$fill[45:271] <- "pink"

bar_data$y <- 1  

ggplot(bar_data, aes(x = pos, y = y, fill = fill)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_identity() +
  theme_void() +
  coord_fixed(ratio = 1/10)  # 控制长宽比

