data_FpuA <- read.table('../input/blast_result/FpuA_results.txt',sep = "\t")
colnames(data_FpuA) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data_FatB <- read.table('../input/blast_result/FatB_results.txt',sep = "\t")
colnames(data_FatB) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data_DesE <- read.table('../input/blast_result/DesE_results.txt',sep = "\t")
colnames(data_DesE) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data_CdtB <- read.table('../input/blast_result/CdtB_results.txt',sep = "\t")
colnames(data_CdtB) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data_FeuA <- read.table('../input/blast_result/FeuA_results.txt',sep = "\t")
colnames(data_FeuA) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")
data_HtbH <- read.table('../input/blast_result/HtbH_results.txt',sep = "\t")
colnames(data_HtbH) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore")
valid_FpuA <- data_FpuA[data_FpuA$pident > 40 & data_FpuA$length > 200,]
valid_FpuA<- valid_FpuA[!grepl("^C", valid_FpuA$sseqid), ]
valid_FpuA$strain_name <- sub('^([^_]+_[^_]+)_.+$','\\1',valid_FpuA$sseqid)
valid_FpuA <- valid_FpuA[!duplicated(valid_FpuA$strain_name),]
valid_FatB <- data_FatB[data_FatB$pident > 40 & data_FatB$length > 200,]
valid_FatB <- valid_FatB[!grepl("^C", valid_FatB$sseqid), ]
valid_FatB$strain_name <- sub("^([^_]+_[^_]+)_.+$", "\\1", valid_FatB$sseqid)
valid_FatB <- valid_FatB[!duplicated(valid_FatB$strain_name),]
valid_DesE <- data_DesE[data_DesE$pident > 40 & data_DesE$length > 200,]
valid_DesE <- valid_DesE[!grepl("^C", valid_DesE$sseqid), ]
valid_DesE$strain_name <- sub("^([^_]+_[^_]+)_.+$", "\\1", valid_DesE$sseqid)
valid_DesE <- valid_DesE[!duplicated(valid_DesE$strain_name),]
valid_CdtB <- data_CdtB[data_CdtB$pident > 40 & data_CdtB$length > 200,]
valid_CdtB <- valid_CdtB[!grepl("^C", valid_CdtB$sseqid), ]
valid_CdtB$strain_name <- sub("^([^_]+_[^_]+)_.+$", "\\1", valid_CdtB$sseqid)
valid_CdtB <- valid_CdtB[!duplicated(valid_CdtB$strain_name),]
valid_FeuA<- data_FeuA[data_FeuA$pident > 40 & data_FeuA$length > 200,]
valid_FeuA <- valid_FeuA[!grepl("^C", valid_FeuA$sseqid), ]
valid_FeuA$strain_name <- sub("^([^_]+_[^_]+)_.+$", "\\1", valid_FeuA$sseqid)
valid_FeuA <- valid_FeuA[!duplicated(valid_FeuA$strain_name),]
valid_HtbH<- data_HtbH[data_HtbH$pident > 40 & data_HtbH$length > 200,]
valid_HtbH <- valid_HtbH[!grepl("^C", valid_HtbH$sseqid), ]
valid_HtbH$strain_name <- sub("^([^_]+_[^_]+)_.+$", "\\1", valid_HtbH$sseqid)
valid_HtbH <- valid_HtbH[!duplicated(valid_HtbH$strain_name),]
sider_info <- read.csv('../input/sidero_info_accesion.csv')

sider_info$Most.similar.known.cluster[sider_info$Most.similar.known.cluster == 'desferrioxamin B/desferrioxamine E'| sider_info$Most.similar.known.cluster == 'desferrioxamin B'| sider_info$Most.similar.known.cluster == 'desferrioxamine E'| sider_info$Most.similar.known.cluster == 'legonoxamine A/desferrioxamine B/legonoxamine B' ] = 'desferrioxamine'
sider_info$Most.similar.known.cluster[sider_info$Most.similar.known.cluster == 'bacillibactin/bacillibactin E/bacillibactin F'] = 'bacillibactin'
sider_info$Most.similar.known.cluster[sider_info$Most.similar.known.cluster == 'heterobactin A/heterobactin S2'| sider_info$Most.similar.known.cluster == 'heterobactin B/heterobactin S2' ] = 'heterobactin'
desferrioxamine_info <- sider_info[sider_info$Most.similar.known.cluster == 'desferrioxamine',]

petrobactin_info <- sider_info[sider_info$Most.similar.known.cluster == 'petrobactin',]
petrobactin_info <- petrobactin_info[!duplicated(petrobactin_info$strainName),]
heterobactin_info <- sider_info[sider_info$Most.similar.known.cluster == 'heterobactin',]
#heterobactin_info <- heterobactin_info[!duplicated(heterobactin_info$strainName),]
bacillibactin_info <- sider_info[sider_info$Most.similar.known.cluster == 'bacillibactin',]
bacillibactin_info <- bacillibactin_info[!duplicated(bacillibactin_info$strainName),]
#bacillibactin_duplicated <- bacillibactin_info[duplicated(bacillibactin_info$strainName),]
rec_info <- read.csv('../input/rec_with_accession_version.csv')
rec_info_distance <- read.csv('../input/rec_in_5genera.csv')
rec_info_distance$VERSION <- rec_info$VERSION[match(rec_info_distance$Tag,rec_info$Tag)]
rec_info_distance$gbk_path <- rec_info$gbk_path[match(rec_info_distance$Tag,rec_info$Tag)]
rec_info$Border_Start <- rec_info_distance$Border_Start[match(rec_info$Tag, rec_info_distance$Tag)]
rec_info$Border_End <- rec_info_distance$Border_End[match(rec_info$Tag, rec_info_distance$Tag)]
rec_info$Strain_Name <- rec_info_distance$Strain_Name[match(rec_info$Tag, rec_info_distance$Tag)]
rec_info$locus_tag <- paste(rec_info$Strain_Name, rec_info$Tag, sep = "_")
HtbH_info <- rec_info[rec_info$locus_tag %in% valid_HtbH$sseqid,]
DesE_info <- rec_info[rec_info$locus_tag %in% valid_DesE$sseqid,]
CdtB_info <- rec_info[rec_info$locus_tag %in% valid_CdtB$sseqid,]
FeuA_info <- rec_info[rec_info$locus_tag %in% valid_FeuA$sseqid,]
FpuA_info <- rec_info[rec_info$locus_tag %in% valid_FpuA$sseqid,]
FatB_info <- rec_info[rec_info$locus_tag %in% valid_FatB$sseqid,]

calculate_strain_distances <- function(data1, data2) {
  
  if (!"strainName" %in% colnames(data1)) {
    stop("data1中缺少strainName列")
  }
  if (!"Strain_Name" %in% colnames(data2)) {
    stop("data2中缺少Strain_Name列")
  }
  if (!"locus_tag" %in% colnames(data2)) {
    stop("data2中缺少locus_tag列")
  }
  common_strains <- intersect(data1$strainName, data2$Strain_Name)
  
  if (length(common_strains) == 0) {
    warning("没有找到重合的strain名称")
    return(data.frame(
      sider_index = integer(0),
      strain_name = character(0),
      receptor_count = integer(0),
      locus_tag = character(0),
      absolute_distance = numeric(0),
      relative_distance = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  data1_filtered <- data1[data1$strainName %in% common_strains, ]
  data2_filtered <- data2[data2$Strain_Name %in% common_strains, ]
  
  results <- list()
  result_index <- 1
  
  for (i in 1:nrow(data1_filtered)) {
    sider_info <- data1_filtered[i, ]
    current_strain <- sider_info$strainName
    
    matching_receptors <- data2_filtered[data2_filtered$Strain_Name == current_strain, ]
    
    if (nrow(matching_receptors) == 0) {
      results[[result_index]] <- data.frame(
        sider_index = i,
        strain_name = current_strain,
        receptor_count = 0,
        locus_tag = NA,
        absolute_distance = NA,
        relative_distance = NA,
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1
      next
    }
    
    bgc_start <- as.numeric(sider_info$start) 
    bgc_end <- as.numeric(sider_info$end)     
    bgc_length <- as.numeric(sider_info$length) 
    
   
    if (is.na(bgc_start) || is.na(bgc_end) || is.na(bgc_length) || bgc_length <= 0) {
     
      for (j in 1:nrow(matching_receptors)) {
        receptor <- matching_receptors[j, ]
        results[[result_index]] <- data.frame(
          sider_index = i,
          strain_name = current_strain,
          receptor_count = nrow(matching_receptors),
          locus_tag = as.character(receptor$locus_tag),
          absolute_distance = NA,
          relative_distance = NA,
          stringsAsFactors = FALSE
        )
        result_index <- result_index + 1
      }
      next
    }
    
    for (j in 1:nrow(matching_receptors)) {
      receptor <- matching_receptors[j, ]
     
      receptor_start <- as.numeric(receptor$Border_Start)  
      receptor_end <- as.numeric(receptor$Border_End)     
     
      if (is.na(receptor_start) || is.na(receptor_end)) {
       
        results[[result_index]] <- data.frame(
          sider_index = i,
          strain_name = current_strain,
          receptor_count = nrow(matching_receptors),
          locus_tag = as.character(receptor$locus_tag),
          absolute_distance = NA,
          relative_distance = NA,
          stringsAsFactors = FALSE
        )
        result_index <- result_index + 1
        next
      }
      
      absolute_distance <- 0
     
      if (receptor_start <= bgc_end && receptor_end >= bgc_start) {
       
        absolute_distance <- 0
      } else {
        
        dist1 <- abs(receptor_start - bgc_start)
        dist2 <- abs(receptor_start - bgc_end)
        dist3 <- abs(receptor_end - bgc_start)
        dist4 <- abs(receptor_end - bgc_end)
        
        absolute_distance <- min(dist1, dist2, dist3, dist4)
      }
      
      relative_distance <- absolute_distance / bgc_length
      
      results[[result_index]] <- data.frame(
        sider_index = i,
        strain_name = current_strain,
        receptor_count = nrow(matching_receptors),
        locus_tag = as.character(receptor$locus_tag),
        absolute_distance = absolute_distance,
        relative_distance = relative_distance,
        stringsAsFactors = FALSE
      )
      result_index <- result_index + 1
    }
  }
 
  if (length(results) > 0) {
    final_result <- do.call(rbind, results)
    rownames(final_result) <- NULL
    return(final_result)
  } else {
    return(data.frame(
      sider_index = integer(0),
      strain_name = character(0),
      receptor_count = integer(0),
      locus_tag = character(0),
      absolute_distance = numeric(0),
      relative_distance = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
}

heterobactin <- calculate_strain_distances(heterobactin_info,HtbH_info)
desferrioxamine_DesE <-  calculate_strain_distances(desferrioxamine_info,DesE_info)
desferrioxamine_CdtB<-  calculate_strain_distances(desferrioxamine_info,CdtB_info)
bacillibactin <- calculate_strain_distances(bacillibactin_info,FeuA_info)
petrobactin_FpuA <- calculate_strain_distances(petrobactin_info,FpuA_info)
petrobactin_FatB <- calculate_strain_distances(petrobactin_info,FatB_info)

library(dplyr)
library(ggplot2)

heterobactin$group <- "heterobactin"
desferrioxamine_DesE$group <- "desferrioxamine_DesE"
desferrioxamine_CdtB$group <- "desferrioxamine_CdtB"
bacillibactin$group <- "bacillibactin"
petrobactin_FpuA$group <- "petrobactin_FpuA"
petrobactin_FatB$group <- "petrobactin_FatB"

all_data <- bind_rows(
  heterobactin,
  desferrioxamine_DesE,
  desferrioxamine_CdtB,
  bacillibactin,
  petrobactin_FpuA,
  petrobactin_FatB
)

all_data$group <- factor(all_data$group,
                         levels = c(
                                    "desferrioxamine_DesE", "desferrioxamine_CdtB",
                                    "bacillibactin","heterobactin",
                                    "petrobactin_FpuA", "petrobactin_FatB")
)

p <- ggplot(all_data, aes(x = log10(absolute_distance + 10))) +
  geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +  # 两列布局
  labs(
    #x = "Distance to nearest PBP2 gene (log10bp)",
    y = "Strain number"
  ) +
  theme_minimal()

p
#ggsave("../Figure/self_receptor_distance_by_cluster.pdf", width = 4, height = 6.5)
write.csv(all_data,'../output/self_rec_distance_info.csv')
sider_info_with_nearest <- read.csv('../output/sider_info_with_nearest_rec.csv')
# 计算overlap ratio
datasets <- list(
  heterobactin = heterobactin,
  desferrioxamine_DesE = desferrioxamine_DesE,
  desferrioxamine_CdtB = desferrioxamine_CdtB,
  bacillibactin = bacillibactin,
  petrobactin_FpuA = petrobactin_FpuA,
  petrobactin_FatB = petrobactin_FatB
)

# 计算overlap ratios
overlap_results <- sapply(datasets, function(data) {
  overlap_count <- length(intersect(data$locus_tag, sider_info_with_nearest$nearest_receptor_index))
  total_count <- nrow(data)
  ratio <- overlap_count / total_count
  return(c(overlap_count = overlap_count, total_count = total_count, ratio = ratio))
})

library(reshape2)
overlap_df <- data.frame(t(overlap_results))
overlap_df$dataset <- rownames(overlap_df)
stacked_data <- overlap_df %>%
  mutate(overlap = overlap_count,
         non_overlap = total_count - overlap_count) %>%
  select(dataset, overlap, non_overlap) %>%
  melt(id.vars = "dataset", variable.name = "category", value.name = "count")



library(dplyr)
library(ggplot2)

stacked_data <- stacked_data %>%
  mutate(category = factor(category, levels = c( "non_overlap","overlap")))


stacked_data_pct <- stacked_data %>%
  group_by(dataset) %>%
  mutate(count_pct = count / sum(count)) %>%
  ungroup()

dataset_order <- c("petrobactin_FpuA", "petrobactin_FatB","bacillibactin","desferrioxamine_CdtB","desferrioxamine_DesE", "heterobactin")

stacked_data_pct <- stacked_data_pct %>%
  mutate(dataset = factor(dataset, levels = dataset_order))

ggplot(stacked_data_pct, aes(x = dataset, y = count_pct, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("overlap" = "steelblue", "non_overlap" = "lightgray")) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Overlap vs Non-overlap Proportion",
       x = "Dataset", y = "Percentage") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())   
#ggsave("../Figure/overlap_self_nearest_receptor_distance_by_cluster.pdf", width = 7, height = 5)
