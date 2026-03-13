data <- read.csv("../input/genome_pbp2_summary.csv",
                 stringsAsFactors = FALSE)

strain_is_PBP <- unique(data$Strain_Name)

species_info <- read.csv("../input/NCBI_all_gtdb_meta_info.csv",
                         stringsAsFactors = FALSE)
gram_pos_strain <- read.csv("../input/genome_bama_screening.csv",
                           stringsAsFactors = FALSE)
gram_pos_strain <- gram_pos_strain[gram_pos_strain$Has_BamA =='False',]
gram_pos_strain$Genus <- species_info$Genus[match(gram_pos_strain$Strain_Name,species_info$Strain_name)]

test <- setdiff(gram_pos_strain$Strain_Name,strain_is_PBP)
species_noPBP <- species_info[species_info$Strain_name %in% test,]
species_noPBP <- species_noPBP[!duplicated(species_noPBP$Strain_name),]
A = as.data.frame(table(species_noPBP$Genus))
B = as.data.frame(table(gram_pos_strain$Genus))
colnames(A) <- c("Genus", "NoPBP_count")
colnames(B) <- c("Genus", "Total_count")

merged_df <- merge(
  A,
  B,
  by = "Genus",
  all.x = TRUE   
)
merged_df$non_PBP_ratio = merged_df$NoPBP_count/merged_df$Total_count
merged_df_high_ratio <- merged_df[merged_df$NoPBP_count > 10 & merged_df$non_PBP_ratio > 0.5,]
sum(merged_df_high_ratio$NoPBP_count)
library(dplyr)
library(ggplot2)
plot_data <- merged_df_high_ratio %>%
  mutate(Genus_Clean = gsub("g__", "", Genus)) %>%
  arrange(non_PBP_ratio) %>%
  mutate(Genus_Clean = factor(Genus_Clean, levels = Genus_Clean)) %>%
 
  mutate(Category = case_when(
    
    grepl("Bifido|Clostridium", Genus, ignore.case = TRUE) ~ "Anaerobes (Fe2+ rich)",
    
    grepl("Myco|Spiro|Urea|Phyto|Meso|Plasma", Genus, ignore.case = TRUE) ~ "Mollicutes (Parasitic)",
    
    grepl("Lacto|Lacti|Pedio|Leucono", Genus, ignore.case = TRUE) ~ "Lactobacillaceae (Mn-centric)",
    
    TRUE ~ "Other"
  ))

p <- ggplot(plot_data, aes(x = Genus_Clean, y = non_PBP_ratio, color = Category)) +
  
  geom_segment(aes(x = Genus_Clean, xend = Genus_Clean, y = 0, yend = non_PBP_ratio), 
               size = 1.2) +
 
  geom_point(aes(size = Total_count), alpha = 1) +
  
  geom_text(aes(label = scales::percent(non_PBP_ratio, accuracy = 0.1)), 
            hjust = -0.4, size = 3.5, color = "black", fontface = "bold") +
 
  coord_flip(ylim = c(0, 1.2)) + 
  scale_y_continuous(labels = scales::percent) +
  
  scale_color_manual(values = c(
    "Anaerobes (Fe2+ rich)" = "#E64B35",           
    "Lactobacillaceae (Mn-centric)" = "#3C5488",  
    "Mollicutes (Parasitic)" = "#00A087",          
    "Other" = "#999999"                            
  )) +
  

  scale_size_continuous(range = c(3, 10), breaks = c(20, 50, 100, 300)) +
 
  theme_minimal() +
  labs(#title = "Evolutionary Loss of PBP2 Genes in Specific Genera",
       #subtitle = "Categorized by ecological strategy and iron dependency",
       x = "",
       y = "Proportion of Strains Lacking PBP2",
       color = "Biological Category", 
       size = "Sample Size (N)") +    
  

  theme(
   
    legend.position = "right", 
    legend.box = "vertical",
    
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  guides(
    
    color = guide_legend(override.aes = list(size = 5)),
    
    size = guide_legend(override.aes = list(color = "gray30"))
  )


print(p)


ggsave("../Figure/PBP2_absence_lollipop_refined.pdf", width = 11, height = 8)

