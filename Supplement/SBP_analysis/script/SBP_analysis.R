data <- read.csv('../output/summ_genus.csv',row.names = 1)
data$filename <- paste0(sub("^[^_]+_[^_]+_", "", data$dataname), ".gbk")
data$path <- sub("index\\.html$", "", data$html_path)
data$gbk_address <- paste0(data$path, data$filename)
write.csv(data,'../output/all_siderophoreBGC.csv')
###随后转到GBK_extract_feature.ipynb里面提取domain信息
data <- read.csv('../output/all_siderophore_BGC_all_domain_info.csv')
all_domain <- as.data.frame(table(data$aSDomain))
SBP_domain <- c("SBP_bac_5", "SBP_bac_3", "Peripla_BP_2", "ZnuA", "SBP_bac_1", "Bmp", "Lipoprotein_9", "DctP", "OpuAC", "PBP_like_2", "Phosphonate_bd", "SBP_bac_6", "Peripla_BP_4", "SBP_bac_8", "Peripla_BP_6", "SBP_bac_11", "NMT1_3", "NMT1", "NosD", "ABC_sub_bind", "TctC", "PBP_like", "Peripla_BP_1", "Cypl", "Peripla_BP_7", "Peripla_BP_5", "LppC", "DUF3798", "DUF3834", "SBP_bac_10", "ANF_receptor", "ABC_transp_aux", "Peripla_BP_3"
)
library(dplyr)
SBP_data <- data[data$aSDomain %in% SBP_domain, ]
#计算共存信息
summary_list <- list()
coexist_list <- list()
##先提取了这些siderophore BGC中的所有含有SBPdomain的基因，然后查看有多少基因还含有其他domain
for (dom in unique(SBP_data$aSDomain)) {
  
  domain_rows <- SBP_data[SBP_data$aSDomain == dom, ]
  locus_tags <- unique(domain_rows$locus_tag)
  
  coexist_count <- 0
  coexist_domain_list <- c()
  
  for (ltag in locus_tags) {
    domains_in_locus <- data$aSDomain[data$locus_tag == ltag]
    others <- setdiff(unique(domains_in_locus), dom)
    
    if (length(others) > 0) {
      coexist_count <- coexist_count + 1
      coexist_domain_list <- c(coexist_domain_list, others)
    }
  }
  
  total_count <- length(unique(domain_rows$locus_tag))
  coexist_ratio <- round(coexist_count / total_count, 3)
  
  # 添加到 summary 列表
  summary_list[[length(summary_list) + 1]] <- data.frame(
    domain = dom,
    total = total_count,
    coexist = coexist_count,
    coexist_ratio = coexist_ratio
  )
  
  # 添加到 coexist 列表（按频率展开）
  if (coexist_count > 0) {
    coexist_freq_table <- table(coexist_domain_list)
    coexist_freq <- as.numeric(coexist_freq_table) / coexist_count
    coexist_names <- names(coexist_freq_table)
    
    coexist_list[[length(coexist_list) + 1]] <- data.frame(
      domain = dom,
      coexist_with = coexist_names,
      freq = round(coexist_freq, 3)
    )
  }
}
summary_df <- do.call(rbind, summary_list)
coexist_df <- do.call(rbind, coexist_list)
##随后在excel里面把这两个表格整理到了一起
write.csv(summary_df, "../output/SBP/all_SBP_summary.csv", row.names = FALSE)
write.csv(coexist_df, "../output/SBP/SBP_coexist_detail.csv", row.names = FALSE)

###再看每种SBP在多少BGC中有出现
SBP_data <- data[data$aSDomain %in% SBP_domain, ]
bgc_sbp_presence <- SBP_data %>%
  select(id, aSDomain) %>%
  distinct() %>%                      # 每个BGC-domain只保留一次
  mutate(present = 1) %>%              # 标记 presence
  tidyr::pivot_wider(
    names_from = aSDomain,
    values_from = present,
    values_fill = 0
  )
siderophore_BGC_meta_data <- read.csv('../output/summ_genus.csv',row.names = 1)

#high_confidence_BGC_info <- siderophore_BGC_meta_data[siderophore_BGC_meta_data$Similarity>=0.8,]
#high_confidence_BGC_name <- as.data.frame(table(high_confidence_BGC_info$Most.similar.known.cluster))
#write.csv(high_confidence_BGC_name,'../output/high_confidence_BGC_name.csv')
siderophore_name <- c('coelichelin', 'griseobactin', 'scabichelin', 'albomycin delta2','desferrioxamin B','desferrioxamine',
                      'desferrioxamin B/desferrioxamine E', 'dehydroxynocardamine','desferrioxamine E','fuscachelin A/fuscachelin B/fuscachelin C',
                      'heterobactin A/heterobactin S2', 'heterobactin B/heterobactin S2','legonoxamine A/desferrioxamine B/legonoxamine B',
                      'petrobactin', 'bacillibactin/bacillibactin E/bacillibactin F','bacillibactin','paenibactin','EDHA','qinichelins',
                      'staphyloferrin B', 'staphyloferrin A', 'peucechelin')
#siderophore_info <- read.csv('../input/SIDERITE_unique_structures_20250131.csv')
high_confidence_BGC <- siderophore_BGC_meta_data$dataname[
  siderophore_BGC_meta_data$Similarity >= 0.8 &
    siderophore_BGC_meta_data$Most.similar.known.cluster %in% siderophore_name
]
#看一下对应的菌株数目
high_confidence_BGC_strain_name <- sub("_N.*", "", high_confidence_BGC)
high_confidence_BGC_strain_name <- high_confidence_BGC_strain_name[!duplicated(high_confidence_BGC_strain_name)]
library(stringr)
#siderophore_info$Siderophore.name <- str_replace(siderophore_info$Siderophore.name, "^([A-Z])", ~ tolower(.x))
#A = as.data.frame(table(high_confidence_BGC_info$Most.similar.known.cluster))
#B = intersect(A$Var1,siderophore_info$Siderophore.name)
#siderophore_name = c(B,'albomycin delta2','bacillibactin/bacillibactin E/bacillibactin F','dehydroxynocardamine','desferrioxamin B','desferrioxamin B/desferrioxamine E','desferrioxamine','griseobactin','heterobactin A/heterobactin S2','heterobactin B/heterobactin S2')
#high_confidence_BGC_info_detail <- high_confidence_BGC_info[high_confidence_BGC_info$Most.similar.known.cluster %in% siderophore_name,]
#这俩没啥区别啊 就差了100+

# 1. 统计 total appearance
total_appearance_df <- SBP_data %>%
  select(id, aSDomain) %>%
  distinct() %>%
  group_by(aSDomain) %>%
  summarise(total_appearance = n(), .groups = "drop")

# 2. 统计 high_confidence_BGC 出现次数
high_conf_count_df <- SBP_data %>%
  select(id, aSDomain) %>%
  distinct() %>%
  filter(id %in% high_confidence_BGC) %>%
  group_by(aSDomain) %>%
  summarise(high_conf_count = n(), .groups = "drop")
library(tidyr)  # 加载tidyr包

# 3. 合并两个结果
sbp_domain_stats <- total_appearance_df %>%
  left_join(high_conf_count_df, by = "aSDomain") %>%
  mutate(high_conf_count = replace_na(high_conf_count, 0)) # 没有匹配的填0
write.csv(sbp_domain_stats,'../output/SBP/SBP_33_info.csv')
###101个的信息手动在excel里面添加了



##查看各种Most.similar.known.cluster SBP_domain出现的比例
library(dplyr)
# 把 high_confidence_BGC_info 里的 Most.similar.known.cluster 关联到 SBP_data
high_conf_with_domain <- SBP_data %>%
  select(id, aSDomain) %>%
  distinct() %>%
  inner_join(
    high_confidence_BGC_info %>% select(dataname, Most.similar.known.cluster),
    by = c("id" = "dataname")
  )

# 按 Most.similar.known.cluster + aSDomain 统计出现次数
domain_count_by_cluster <- high_conf_with_domain %>%
  group_by(Most.similar.known.cluster, aSDomain) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(Most.similar.known.cluster, desc(count))
domain_count_PBP2 <- domain_count_by_cluster[domain_count_by_cluster$aSDomain == 'Peripla_BP_2',]
