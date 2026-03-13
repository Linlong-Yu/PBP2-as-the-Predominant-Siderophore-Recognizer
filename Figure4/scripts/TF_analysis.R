data <- read.csv('../output/self_rec_distance_info.csv')
bacillibactin <- data[data$group == 'bacillibactin',]
petrobactin_FpuA <- data[data$group == 'petrobactin_FpuA',]
petrobactin_FatB <- data[data$group == 'petrobactin_FatB',]
desferrioxamine_DesE <- data[data$group == 'desferrioxamine_DesE',]
desferrioxamine_CdtB <- data[data$group == 'desferrioxamine_CdtB',]
distant_bacillibactin <- bacillibactin[bacillibactin$absolute_distance > 0,]
distant_petrobactin_FpuA <- petrobactin_FpuA[petrobactin_FpuA$absolute_distance > 0,]
distant_petrobactin_FatB <- petrobactin_FatB[petrobactin_FatB$absolute_distance > 0,]
distant_desferrioxamine_DesE <- desferrioxamine_DesE[desferrioxamine_DesE$absolute_distance > 6000,]
distant_desferrioxamine_CdtB <- desferrioxamine_CdtB[desferrioxamine_CdtB$absolute_distance > 0,]
Fur_fimo_result <- read.csv('../input/fimo_result_Fur.csv')
Fur_fimo_result <- Fur_fimo_result[!duplicated(Fur_fimo_result$sequence_name),]
Fur_fimo_result$Locus_tag_complete <- paste0(
  sub("\\.gbk$", "", Fur_fimo_result$strainname), 
  "_",
  Fur_fimo_result$Locus_id
)
DmdR1_fimo_result <- read.csv('../input/fimo_result_DmdR1.csv')
hist(DmdR1_fimo_result$score)
DmdR1_fimo_result <- DmdR1_fimo_result[!duplicated(DmdR1_fimo_result$sequence_name),]
DmdR1_fimo_result$locus <- sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", DmdR1_fimo_result$sequence_name)
#write.csv(DmdR1_fimo_result,'../output/DmdR1_fimo_result.csv')
DmdR1_fimo_result$strainname <- sub("\\.gbk.*$", "", DmdR1_fimo_result$sequence_name)
DmdR1_fimo_result$Locus_tag_complete <- paste0(DmdR1_fimo_result$strainname,"_",DmdR1_fimo_result$locus)

FeuA_Fur <- intersect(Fur_fimo_result$Locus_tag_complete,distant_bacillibactin$locus_tag)
FpuA_Fur <- intersect(Fur_fimo_result$Locus_tag_complete,distant_petrobactin_FpuA$locus_tag)
FatB_Fur <- intersect(Fur_fimo_result$Locus_tag_complete,distant_petrobactin_FatB$locus_tag)
strain_FpuA <- sub("^([^_]*_[^_]*)_.*$", "\\1", FpuA_Fur) 
strain_FatB <- sub("^([^_]*_[^_]*)_.*$", "\\1", FatB_Fur) 
intersect(strain_FatB,strain_FpuA)
Des_DmdR1 <- intersect(DmdR1_fimo_result$Locus_tag_complete,distant_desferrioxamine_DesE$locus_tag)
data_operon = read.csv('../input/operon_analysis_all_genome.csv')
data_operon$operon_Locus_tag_complete <- paste0(
  sub("\\.gbk$", "", data_operon$filename), 
  "_",
  data_operon$operon_first_locus
)
data_operon$Locus_tag_complete <- paste0(
  sub("\\.gbk$", "", data_operon$filename), 
  "_",
  data_operon$locus_tag
)
data_operon_clean <- data_operon[!is.na(data_operon$operon_first_ocus), ]
distant_desferrioxamine_CdtB$locus_operon <- data_operon$operon_Locus_tag_complete[match(distant_desferrioxamine_CdtB$locus_tag,data_operon$Locus_tag_complete)]
DesE_DmdR1 <- intersect(DmdR1_fimo_result$Locus_tag_complete,distant_desferrioxamine_DesE$locus_tag)
CdtB_DmdR1 <- intersect(DmdR1_fimo_result$Locus_tag_complete,distant_desferrioxamine_CdtB$locus_operon)
distant_desferrioxamine_CdtB_noDesE <- distant_desferrioxamine_CdtB[!desferrioxamine_CdtB$strain_name %in% desferrioxamine_DesE$strain_name,]
distant_desferrioxamine_CdtB_haveDesE <- distant_desferrioxamine_CdtB[desferrioxamine_CdtB$strain_name %in% desferrioxamine_DesE$strain_name,]
CdtB_DmdR1_noDesE <- intersect(DmdR1_fimo_result$Locus_tag_complete,distant_desferrioxamine_CdtB_noDesE$locus_operon)
CdtB_DmdR1_haveDesE <- intersect(DmdR1_fimo_result$Locus_tag_complete,distant_desferrioxamine_CdtB_haveDesE$locus_operon)

bacillus_BGC_fimo_results <- read.csv('../input/detect_all_BGC_all_genome_Fur.csv')
bacillus_BGC_fimo_results$dataname <- paste(bacillus_BGC_fimo_results$Subfolder,bacillus_BGC_fimo_results$File,sep = '_')
bacillus_BGC_fimo_results$dataname <- sub("\\.gbk$", "", bacillus_BGC_fimo_results$dataname)
bacillus_BGC_fimo_results_Fur <- bacillus_BGC_fimo_results[bacillus_BGC_fimo_results$Subset == 1,]
Str_BGC_fimo_results <- read.csv('../input/detect_all_BGC_all_genome_DmdR1.csv')
Str_BGC_fimo_results$dataname <- paste(Str_BGC_fimo_results$Subfolder,Str_BGC_fimo_results$File,sep = '_')
Str_BGC_fimo_results$dataname <- sub("\\.gbk$", "",Str_BGC_fimo_results$dataname)
Str_BGC_fimo_results_Fur <- Str_BGC_fimo_results[Str_BGC_fimo_results$Subset == 1,]
sider_info <- read.csv('../input/sidero_info_accesion.csv')
sider_info$Most.similar.known.cluster[sider_info$Most.similar.known.cluster == 'desferrioxamin B/desferrioxamine E'| sider_info$Most.similar.known.cluster == 'desferrioxamin B'| sider_info$Most.similar.known.cluster == 'desferrioxamine E'| sider_info$Most.similar.known.cluster == 'legonoxamine A/desferrioxamine B/legonoxamine B' ] = 'desferrioxamine'
sider_info$Most.similar.known.cluster[sider_info$Most.similar.known.cluster == 'bacillibactin/bacillibactin E/bacillibactin F'] = 'bacillibactin'
petrobactin_info <- sider_info[sider_info$Most.similar.known.cluster == 'petrobactin',]
#有一个错了，所以要去重，其他的检查过没问题
petrobactin_info <- petrobactin_info[!duplicated(petrobactin_info$strainName),]
bacillibactin_info <- sider_info[sider_info$Most.similar.known.cluster == 'bacillibactin',]
desferrioxamine_info <- sider_info[sider_info$Most.similar.known.cluster == 'desferrioxamine',]
library(dplyr)
distant_bacillibactin_syn <- bacillibactin_info[bacillibactin_info$strainName %in% distant_bacillibactin$strain_name,]
distant_petrobactin_FpuA_syn <- petrobactin_info[petrobactin_info$strainName %in% distant_petrobactin_FpuA$strain_name,]
distant_petrobactin_FatB_syn <- petrobactin_info[petrobactin_info$strainName %in% distant_petrobactin_FatB$strain_name,]
distant_desferrioxamine_DesE_syn <- desferrioxamine_info[desferrioxamine_info$strainName %in% distant_desferrioxamine_DesE$strain_name,]
distant_desferrioxamine_CdtB_noDesE_syn <- desferrioxamine_info[desferrioxamine_info$strainName %in% distant_desferrioxamine_CdtB_noDesE$strain_name,]
distant_desferrioxamine_CdtB_haveDesE_syn <- desferrioxamine_info[desferrioxamine_info$strainName %in% distant_desferrioxamine_CdtB_haveDesE$strain_name,]

bacillibactin_Fur <- intersect(distant_bacillibactin_syn$dataname,bacillus_BGC_fimo_results_Fur$dataname)
Petrobactin_FpuA_Fur <- intersect(distant_petrobactin_FpuA_syn$dataname,bacillus_BGC_fimo_results_Fur$dataname)
Petrobactin_FatB_Fur <- intersect(distant_petrobactin_FatB_syn$dataname,bacillus_BGC_fimo_results_Fur$dataname)
desferrioxamine_DesE_DmdR1 <- intersect(distant_desferrioxamine_DesE_syn$dataname,Str_BGC_fimo_results_Fur$dataname)
desferrioxamine_CdtB_noDesE_DmdR1 <- intersect(distant_desferrioxamine_CdtB_noDesE_syn$dataname,Str_BGC_fimo_results_Fur$dataname)
desferrioxamine_CdtB_haveDesE_DmdR1 <- intersect(distant_desferrioxamine_CdtB_haveDesE_syn$dataname,Str_BGC_fimo_results_Fur$dataname)


strain_CdtB <- sub("^([^_]*_[^_]*)_.*$", "\\1", CdtB_DmdR1)
near_Des_strain <- desferrioxamine_DesE$strain_name[desferrioxamine_DesE$absolute_distance < 6000]
no_near_Des<- desferrioxamine_info[!desferrioxamine_info$strainName %in% near_Des_strain, ]
no_near_DesE_have_CdtB <- no_near_Des[no_near_Des$strainName %in% desferrioxamine_CdtB$strain_name,]
strain_DmdR1_CdtB <- intersect(strain_CdtB,no_near_DesE_have_CdtB$strainName)

vec_list <- list(bacillibactin_Fur,FeuA_Fur,distant_bacillibactin$strain_name,
                 Petrobactin_FpuA_Fur,FpuA_Fur,distant_petrobactin_FpuA$strain_name,
                 Petrobactin_FatB_Fur,FatB_Fur,distant_petrobactin_FatB$strain_name,
                 desferrioxamine_DesE_DmdR1,DesE_DmdR1,distant_desferrioxamine_DesE$strain_name,
                 desferrioxamine_CdtB_haveDesE_DmdR1,CdtB_DmdR1_haveDesE,distant_desferrioxamine_CdtB_haveDesE$strain_name,
                 desferrioxamine_CdtB_noDesE_DmdR1,CdtB_DmdR1_noDesE,distant_desferrioxamine_CdtB_noDesE$strain_name)

lengths_list <- sapply(vec_list, length)

df <- data.frame(
  lengths_list[seq(1,18,by=3)],
  lengths_list[seq(2,18,by=3)],
  lengths_list[seq(3,18,by=3)]
)
colnames(df) <- c('synthetase','self_receptor','Total strain number')  

rownames(df) <- c('bacillibactin-FeuA','petrobactin-FpuA','petrobactin-FatB','desferrioxamine-DesE','desferrioxamine-CdtB(with DesE)','desferrioxamine-CdtB(without DesE)')
df$synthetase_ratio <- df$synthetase/df$`Total strain number`
df$szelf_receptor_ratio <- df$self_receptor/df$`Total strain number`
library(ggplot2)
library(gridExtra)
library(scales)   

make_plot <- function(df, i) {
  plot_data <- data.frame(
    variable = c("synthetase_ratio", "szelf_receptor_ratio"),
    value = c(df[i, "synthetase_ratio"], df[i, "szelf_receptor_ratio"])
  )
  
  ggplot(plot_data, aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = percent(value, accuracy = 1)),   
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("synthetase_ratio" = "#3498db",
                                 "szelf_receptor_ratio" = "#e74c3c")) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +  
    labs(title = rownames(df)[i], x = "", y = "") +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 8),
      axis.line.x = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

plots <- lapply(1:6, function(i) make_plot(df, i))
grid.arrange(grobs = plots, ncol = 6)
plot
#ggsave("../Figure/TF_ratio.pdf",
  #     arrangeGrob(grobs = plots, ncol = 6),
   #    width = 12, height = 3)

