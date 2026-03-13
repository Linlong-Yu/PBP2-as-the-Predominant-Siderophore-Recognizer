
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

domain_data <- read.table("../input/SBP_domains.dtbl", skip = 3, fill = TRUE, sep = "", 
                          comment.char = "#", header = FALSE)
domain_data <- domain_data[1:22]
domain_data <- domain_data[-5]

colnames(domain_data) <- c(
  "target_name", "target_accession", "tlen", "query_name",
  "qlen", "full_seq_evalue", "full_seq_score", "full_seq_bias", 
  "dom_N", "dom_of", "dom_c_evalue", "dom_i_evalue", "dom_score", "dom_bias", 
  "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",'accurancy')

domain_clean <- domain_data %>%
  group_by(query_name,target_accession) %>%
  slice_min(dom_i_evalue, n = 1) %>%
  ungroup()

domain_clean <- subset(domain_clean, !is.na(tlen))

domain_counts <- domain_clean %>%
  count(target_accession, target_name, sort = TRUE) %>%
  mutate(target_name = fct_reorder(target_name, n))
library(ggplot2)
p <- ggplot(domain_counts, aes(x = n, y = target_name)) +
  geom_col(fill = "#ffab91", width = 0.8) +
  geom_text(aes(label = n), hjust = -0.2,vjust = 0.4, size = 4) + 
  labs(
    x = NULL,  
    y = "Pfam Domain",
    title = "Shared Domains in 22 Iron Uptake SBP Receptors"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),     
    axis.ticks.x = element_blank(),    
    axis.line.x = element_blank(),     
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid = element_blank()      
  ) +
  xlim(0, max(domain_counts$n) * 1.1)  

print(p)
#ggsave("../Figure/shared_domains_plot.pdf", plot = p, device = "pdf", width = 8, height = 6)
