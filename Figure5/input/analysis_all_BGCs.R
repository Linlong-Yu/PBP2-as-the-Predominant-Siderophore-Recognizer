data <- read.csv('../input/all_syn_stat_df.csv')
data_with_TonB <- data[data$genome_tonb > 0,]
data_without_TonB <- data[data$genome_tonb == 0,]
Tonb <- as.data.frame(table(data_with_TonB$tonb_num)) 
PBP2_in_neg <-as.data.frame(table(data_with_TonB$pb_num))
data_without_TonB_NRPS <- data_without_TonB[data_without_TonB$sid_type == 'NRP',]
PBP2 <- as.data.frame(table(data_without_TonB$pb_num))
PBP2_NRPS_inPOS <- as.data.frame(table(data_without_TonB_NRPS$pb_num))

data <- read.csv('../input/bgc_pfam_matrix.csv')
phylogeny_info <- read.csv('../../../gtdb_meta_info/NCBI_all_gtdb_meta_info.csv')


numeric_data <- data[sapply(data, is.numeric)]
numeric_data <- as.data.frame(lapply(numeric_data, function(x) ifelse(x > 1, 1, x)))
col_sums <- colSums(numeric_data)
filtered_data <- numeric_data[, col_sums >= 1]
result <- rbind(Sum = colSums(filtered_data), filtered_data)
data$TonB_dep_Rec
data_transporter <- result[,c('SBP_bac_1','SBP_bac_11','SBP_bac_3','SBP_bac_5','SBP_bac_6','SBP_bac_8','Peripla_BP_1','Peripla_BP_2','Peripla_BP_3','Peripla_BP_4','Peripla_BP_5','Peripla_BP_6','ABC_membrane','ABC_tran','FecCD','TonB_dep_Rec','MFS_3','MFS_1')]
library(ggplot2)
library(tidyr)
df_plot <- data.frame(
  Domain = colnames(data_transporter),
  Value = as.numeric(data_transporter[1, ])
)

df_plot$Value <- df_plot$Value/(nrow(data_transporter)-1)
df_plot <- df_plot[order(df_plot$Value, decreasing = TRUE), ]
df_plot$Domain <- factor(df_plot$Domain, levels = df_plot$Domain)  # 保持排序

# 绘制柱状图
ggplot(df_plot, aes(x = Domain, y = Value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Domain name", y = "proportion")#, title = "First Row Barplot of data_transporter")
ggsave("../Figure/domain_proportion_in_all_BGC.pdf", width = 10, height = 6)
