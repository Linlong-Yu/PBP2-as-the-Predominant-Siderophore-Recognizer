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
sider_info <- read.csv('../input/sider_info.csv')
petrobactin_info <- sider_info[sider_info$Most.similar.known.cluster == 'petrobactin',]
sider_info$Most.similar.known.cluster[sider_info$Most.similar.known.cluster == 'desferrioxamin B/desferrioxamine E'| sider_info$Most.similar.known.cluster == 'desferrioxamin B'| sider_info$Most.similar.known.cluster == 'desferrioxamine E'| sider_info$Most.similar.known.cluster == 'legonoxamine A/desferrioxamine B/legonoxamine B' ] = 'desferrioxamine'
desferrioxamine_info <- sider_info[sider_info$Most.similar.known.cluster == 'desferrioxamine',]
desferrioxamine_info <- desferrioxamine_info[!duplicated(desferrioxamine_info$strainName),]
petrobactin_strain_name <- petrobactin_info$strainName[!duplicated(petrobactin_info$strainName)]
petrobactin_strain_name <- petrobactin_info[petrobactin_info$strainName == 'GCF_030718825.1',]
petrobactin_info <- petrobactin_info[!duplicated(petrobactin_info$strainName),]
A=intersect(petrobactin_info$strainName,valid_FatB$strain_name)
B = intersect(petrobactin_info$strainName,valid_FpuA$strain_name)
C = setdiff(B,A)

venn_fit1 <- 
  list(
    petrobactin = petrobactin_info$strainName,
    FatB = valid_FatB$strain_name,
    FpuA = valid_FpuA$strain_name
  )

#install.packages("ggvenn")
library(ggvenn)

p <- ggvenn(
  venn_fit1,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"), # 好看的调色板
  stroke_size = 0.5,
  set_name_size = 5
)
p

#ggsave("../Figure/Petrobactin_venn_plot.pdf", plot = p, width = 6, height = 6)

venn_fit2 <- 
  list(
    Desferrioxamine =desferrioxamine_info$strainName,
    DesE = valid_DesE$strain_name,
    CdtB = valid_CdtB$strain_name
  )

p <- ggvenn(
  venn_fit2,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"), # 好看的调色板
  stroke_size = 0.5,
  set_name_size = 5
)
p

#ggsave("../Figure/Desferrioxamine_venn_plot.pdf", plot = p, width = 6, height = 6)


