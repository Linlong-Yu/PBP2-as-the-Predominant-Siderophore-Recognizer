data <- read.csv('../input/summ_genus.csv',row.names = 1)
library(stringr)
data$strain_name <- str_extract(data$dataname, "[^_]+_[^_]+")
data$combined <- paste(data$strain_name, data$Most.similar.known.cluster, sep = "_")
library(dplyr)
# filter high confidence siderophore BGC
data = data[data$Similarity >= 0.8,]
siderophore_name <- c('coelichelin', 'griseobactin', 'scabichelin', 'albomycin delta2',
                      'desferrioxamin B/desferrioxamine E', 'dehydroxynocardamine',
                      'heterobactin A/heterobactin S2', 'qinichelins','paenibactin',
                      'petrobactin', 'bacillibactin/bacillibactin E/bacillibactin F',
                      'staphyloferrin B', 'staphyloferrin A', 'peucechelin',
                      'fuscachelin A/fuscachelin B/fuscachelin C','EDHA')
# filter product is siderophore
filtered_data <- data %>%
  filter(Most.similar.known.cluster %in% siderophore_name)
# there are some siderophore BGC copies in a genome,filter
filtered_data = filtered_data[-c(396,399,570,588),]
# filter 10 siderophore BGC for each siderophore type
final_result <- data.frame()
for (siderophore in siderophore_name) {
  subset_data <- filtered_data %>%
    filter(Most.similar.known.cluster == siderophore)
  row_count <- nrow(subset_data)
  
  if (row_count <= 10) {
    selected_data <- subset_data
  } else {
    if (siderophore == 'griseobactin') {
      # 优先选择 syn_group 为 41 的
      preferred <- subset_data %>% filter(syn_group == 41)
      if (nrow(preferred) >= 10) {
        selected_data <- preferred %>% arrange(desc(Similarity)) %>% head(10)
      } else {
        selected_data <- preferred
        remaining_needed <- 10 - nrow(preferred)
        remaining_data <- subset_data %>%
          filter(syn_group != 41) %>%
          arrange(desc(Similarity)) %>%
          head(remaining_needed)
        selected_data <- bind_rows(selected_data, remaining_data)
      }
    } else {
      # 其他 siderophore 的通用策略
      strict_match <- subset_data %>%
        filter(syn_type %in% c("NI-siderophore", 'NRP-metallophore,NRPS'))
      if (nrow(strict_match) >= 10) {
        selected_data <- strict_match %>% arrange(desc(Similarity)) %>% head(10)
      } else {
        selected_data <- strict_match
        remaining_needed <- 10 - nrow(strict_match)
        remaining_data <- subset_data %>%
          filter(!syn_type %in% c("NI-siderophore", "NRP-metallophore,NRPS")) %>%
          arrange(desc(Similarity)) %>%
          head(remaining_needed)
        selected_data <- bind_rows(selected_data, remaining_data)
      }
    }
  }
  final_result <- bind_rows(final_result, selected_data)
}
#save
write.csv(final_result, "../output/filtered_siderophores_BGC.csv", row.names = FALSE)

