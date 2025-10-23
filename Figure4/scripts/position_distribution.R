sider_info <- read.csv('../input/high_confidence_sider_info.csv')
rec_info <- read.csv('../input/rec_with_accession_version.csv')
rec_info_distance <- read.csv('../input/rec_in_5genera.csv')
rec_info_distance$VERSION <- rec_info$VERSION[match(rec_info_distance$Tag,rec_info$Tag)]
rec_info_distance$gbk_path <- rec_info$gbk_path[match(rec_info_distance$Tag,rec_info$Tag)]
rec_info$Border_Start <- rec_info_distance$Border_Start[match(rec_info$Tag, rec_info_distance$Tag)]
rec_info$Border_End <- rec_info_distance$Border_End[match(rec_info$Tag, rec_info_distance$Tag)]
rec_info$Strain_Name <- rec_info_distance$Strain_Name[match(rec_info$Tag, rec_info_distance$Tag)]
# compute distance between siderophore BGC and PBP2

library(dplyr)
calculate_sider_receptor_stats <- function(sider_info, rec_info) {
 
  results <- list()
  
  for (i in 1:nrow(sider_info)) {
    bgc <- sider_info[i, ]
    
    matching_receptors <- rec_info[rec_info$gbk_path == bgc$gbk_path, ]
    
    result <- list(
      sider_index = i,
      dataname = bgc$dataname,  
      gbk_path = bgc$gbk_path,
      receptor_count = nrow(matching_receptors),
      distances = c(),
      relative_distances = c()
    )
    
    if (nrow(matching_receptors) > 0) {
   
      bgc_start <- as.numeric(bgc$start)
      bgc_end <- as.numeric(bgc$end)
      bgc_length <- as.numeric(bgc$length)
      
      if (!is.na(bgc_start) && !is.na(bgc_end) && !is.na(bgc_length) && bgc_length > 0) {
        
        for (j in 1:nrow(matching_receptors)) {
          receptor <- matching_receptors[j, ]
          
          if (as.character(bgc$accesion_version) == as.character(receptor$VERSION)) {
            
            receptor_start <- as.numeric(receptor$Border_Start)
            receptor_end <- as.numeric(receptor$Border_End)
            
            if (!is.na(receptor_start) && !is.na(receptor_end)) {
              
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
              
              if (!is.na(relative_distance)) {
               
                result$distances <- c(result$distances, absolute_distance)
                result$relative_distances <- c(result$relative_distances, relative_distance)
              }
            }
          }
        }
      }
    }
    
    results[[i]] <- result
  }
  
  return(results)
}



convert_results_to_dataframe <- function(results) {
  
  basic_df <- data.frame(
    sider_index = sapply(results, function(x) x$sider_index),
    dataname = sapply(results, function(x) x$dataname),
    gbk_path = sapply(results, function(x) x$gbk_path),
    receptor_count = sapply(results, function(x) x$receptor_count),
    stringsAsFactors = FALSE
  )
  
  detail_list <- list()
  for (i in 1:length(results)) {
    result <- results[[i]]
    if (length(result$distances) > 0) {
      for (j in 1:length(result$distances)) {
        detail_list[[length(detail_list) + 1]] <- data.frame(
          sider_index = result$sider_index,
          dataname = result$dataname,
          gbk_path = result$gbk_path,
          receptor_index = j,
          absolute_distance = result$distances[j],
          relative_distance = result$relative_distances[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(detail_list) > 0) {
    detail_df <- do.call(rbind, detail_list)
  } else {
    detail_df <- data.frame(
      sider_index = integer(0),
      dataname = character(0),
      gbk_path = character(0),
      receptor_index = integer(0),
      absolute_distance = numeric(0),
      relative_distance = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(summary = basic_df, details = detail_df))
}

stats_results <- calculate_sider_receptor_stats(sider_info, rec_info)
formatted_results <- convert_results_to_dataframe(stats_results)
siderophore_rec <- formatted_results$summary
siderophore_rec_detail <- formatted_results$details
siderophore_rec_detail$most_similar_siderophore <- sider_info$Most.similar.known.cluster[match(siderophore_rec_detail$dataname, sider_info$dataname)]
siderophore_rec_detail$most_similar_siderophore[siderophore_rec_detail$most_similar_siderophore == 'bacillibactin/bacillibactin E/bacillibactin F'] = 'bacillibactin'
siderophore_rec_detail$most_similar_siderophore[siderophore_rec_detail$most_similar_siderophore == 'fuscachelin A/fuscachelin B/fuscachelin C'] = 'fuscachelin'
siderophore_rec_detail$most_similar_siderophore[siderophore_rec_detail$most_similar_siderophore == 'heterobactin A/heterobactin S2'| siderophore_rec_detail$most_similar_siderophore == 'heterobactin B/heterobactin S2' ] = 'heterobactin'
siderophore_rec_detail$most_similar_siderophore[siderophore_rec_detail$most_similar_siderophore == 'desferrioxamin B/desferrioxamine E'| siderophore_rec_detail$most_similar_siderophore == 'desferrioxamin B'| siderophore_rec_detail$most_similar_siderophore == 'desferrioxamine E'| siderophore_rec_detail$most_similar_siderophore == 'legonoxamine A/desferrioxamine B/legonoxamine B' ] = 'desferrioxamine'
write.csv(siderophore_rec_detail,'../output/all_sider_rec_distance.csv')
siderophore_plot <- siderophore_rec_detail
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'peucechelin'] = 'peucechelin(1)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'scabichelin'] = 'scabichelin(107)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'dehydroxynocardamine'] = 'dehydroxynocardamine(5)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'heterobactin'] = 'heterobactin(38)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'EDHA'] = 'EDHA(1)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'desferrioxamine'] = 'desferrioxamine(956)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'paenibactin'] = 'paenibactin(38)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'qinichelins'] = 'qinichelins(9)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'albomycin delta2'] = 'albomycin delta2(6)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'coelichelin'] = 'coelichelin(348)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'fuscachelin'] = 'fuscachelin(8)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'griseobactin'] = 'griseobactin(105)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'bacillibactin'] = 'bacillibactin(1699)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'petrobactin'] = 'petrobactin(464)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'staphyloferrin A'] = 'staphyloferrin A(2087)'
siderophore_plot$most_similar_siderophore[siderophore_plot$most_similar_siderophore == 'staphyloferrin B'] = 'staphyloferrin B(1640)'

cluster_levels <- c("heterobactin(38)","peucechelin(1)", "scabichelin(107)",
                    'EDHA(1)','desferrioxamine(956)','paenibactin(38)','qinichelins(9)',                 
                    "albomycin delta2(6)","bacillibactin(1699)", "coelichelin(348)",
                    "fuscachelin(8)", "griseobactin(105)",
                    
                    "petrobactin(464)","dehydroxynocardamine(5)","staphyloferrin A(2087)", "staphyloferrin B(1640)" 
                    
)

siderophore_plot$most_similar_siderophore <- factor(
  siderophore_plot$most_similar_siderophore,
  levels = cluster_levels
)


library(ggplot2)

ggplot(siderophore_plot, aes(x = relative_distance)) + 
  geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
  facet_wrap(~ most_similar_siderophore, scales = "free_y") +
  labs(
    x = "Distance between synthetase and PBP2 genes (log10 bp)",
    y = "Strain number"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(color = "black"),
    axis.text.y  = element_text(color = "black"),
   
    panel.background = element_blank(),
  
    strip.text = element_text(face = "bold", color = "black")
  )


# compute distance of nearest PBP2
calculate_nearest_receptor <- function(sider_info, rec_info) {
  
  sider_info$nearest_receptor_index <- NA
  sider_info$nearest_receptor_sequence <- NA
  sider_info$nearest_receptor_distance <- NA
  sider_info$nearest_receptor_start <- NA
  sider_info$nearest_receptor_end <- NA

  for (i in 1:nrow(sider_info)) {
    bgc <- sider_info[i, ]
    
    matching_receptors <- rec_info[rec_info$gbk_path == bgc$gbk_path, ]
    
    matching_receptors <- matching_receptors[
      as.character(matching_receptors$VERSION) == as.character(bgc$accesion_version), 
    ]
    
    if (nrow(matching_receptors) == 0) {
      next  
    }
    
    bgc_start <- as.numeric(bgc$start)
    bgc_end <- as.numeric(bgc$end)
    
    if (is.na(bgc_start) || is.na(bgc_end)) {
      next
    }
   
    min_distance <- Inf
    nearest_receptor_idx <- NA
    
    for (j in 1:nrow(matching_receptors)) {
      receptor <- matching_receptors[j, ]
      
      receptor_start <- as.numeric(receptor$Border_Start)
      receptor_end <- as.numeric(receptor$Border_End)
      
      if (is.na(receptor_start) || is.na(receptor_end)) {
        next
      }
     
      actual_distance <- 0
      
      if (receptor_start <= bgc_end && receptor_end >= bgc_start) {
        
        actual_distance <- 0
      } else {
        
        dist1 <- abs(receptor_start - bgc_start)
        dist2 <- abs(receptor_start - bgc_end)
        dist3 <- abs(receptor_end - bgc_start)
        dist4 <- abs(receptor_end - bgc_end)
        
        actual_distance <- min(dist1, dist2, dist3, dist4)
      }
     
      if (actual_distance < min_distance) {
        min_distance <- actual_distance
        nearest_receptor_idx <- j
      }
    }
    
    if (!is.na(nearest_receptor_idx)) {
      nearest_receptor <- matching_receptors[nearest_receptor_idx, ]
      sider_info$nearest_receptor_index[i] <- nearest_receptor$Rec_Name
      sider_info$nearest_receptor_sequence[i] <- nearest_receptor$Sequence
      sider_info$nearest_receptor_distance[i] <- min_distance
      sider_info$nearest_receptor_start[i] <- nearest_receptor$Border_Start
      sider_info$nearest_receptor_end[i] <- nearest_receptor$Border_End
    }
  }
  
  return(sider_info)
}


sider_info_with_nearest <- calculate_nearest_receptor(sider_info, rec_info_distance)
write.csv(sider_info_with_nearest,'../output/sider_info_with_nearest_rec.csv')

library(ggplot2)
library(dplyr)

# remove NA
sider_info_plot <- sider_info_with_nearest %>%
  filter(!is.na(nearest_receptor_distance) & !is.na(Most.similar.known.cluster))
sider_info_plot$Most.similar.known.cluster[sider_info_plot$Most.similar.known.cluster == 'bacillibactin/bacillibactin E/bacillibactin F'] = 'bacillibactin'
sider_info_plot$Most.similar.known.cluster[sider_info_plot$Most.similar.known.cluster == 'heterobactin A/heterobactin S2'| sider_info_plot$Most.similar.known.cluster == 'heterobactin B/heterobactin S2' ] = 'heterobactin'
sider_info_plot$Most.similar.known.cluster[sider_info_plot$Most.similar.known.cluster == 'desferrioxamin B/desferrioxamine E'| sider_info_plot$Most.similar.known.cluster == 'desferrioxamin B'| sider_info_plot$Most.similar.known.cluster == 'desferrioxamine E'| sider_info_plot$Most.similar.known.cluster == 'legonoxamine A/desferrioxamine B/legonoxamine B' ] = 'desferrioxamine'
sider_info_plot$Most.similar.known.cluster[sider_info_plot$Most.similar.known.cluster == 'fuscachelin A/fuscachelin B/fuscachelin C'] = 'fuscachelin'
filtered_data <- sider_info_plot[!(sider_info_plot$Most.similar.known.cluster %in% c("bacillibactin","petrobactin")), ]
filtered_data_2w <- filtered_data[filtered_data$nearest_receptor_distance <20000,]
filtered_data_bacillibactin <- sider_info_plot[sider_info_plot$Most.similar.known.cluster=='bacillibactin',]
filtered_data_petrobactin <- sider_info_plot[sider_info_plot$Most.similar.known.cluster=='petrobactin',]

cluster_levels <- c("heterobactin","peucechelin", "scabichelin",
                    'EDHA','desferrioxamine','paenibactin','qinichelins',                 
                    "albomycin delta2", "bacillibactin","coelichelin",
                    "fuscachelin", "griseobactin",
                    "petrobactin","dehydroxynocardamine",
                    "staphyloferrin A", "staphyloferrin B" 
                    
)

sider_info_plot$Most.similar.known.cluster <- factor(
  sider_info_plot$Most.similar.known.cluster,
  levels = cluster_levels
)

p1 <- ggplot(sider_info_plot, aes(x =log10(nearest_receptor_distance + 10))) +  # +1 避免log(0)
  geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
  facet_wrap(~ Most.similar.known.cluster, scales = "free_y") +
  
  labs(
    #title = "Distribution of Nearest Receptor Distance by Siderophore Cluster Type (Log Scale)",
    x = "Distance to nearest PBP2 gene (log10bp)",
    y = "Strain numebr"
  ) +
  theme_minimal()

print(p1)


#ggsave("../Figure/receptor_distance_by_cluster.pdf", width = 12, height = 9)
#ggsave("../Figure/receptor_distance_by_cluster_plot.pdf", width = 12, height = 9)

  
  
  
  
  
  

