# README

## Input Files

- **blast_result/**  
  BLAST results for self-rec (self-receptor) genes.

- **fimo_result_DmdR1.csv** / **fimo_result_Fur.csv**  
  Results of motif scanning using **FIMO**, identifying genes with DmdR1 or Fur transcription factor (TF) binding sites in their upstream regions.

- **detect_all_BGC_all_genome_DmdR1.csv** / **detect_all_BGC_all_genome_Fur.csv**  
  Results showing whether siderophore BGCs contain DmdR1 or Fur TF binding sites within their genomic regions.

- **high_confidence_sider_info.csv**  
  High-confidence siderophore BGC information.

- **rec_in_5genera.csv** / **rec_with_accession_version.csv**  
  Information for PBP2 genes across the five genera analyzed.

- **sidero_info_accesion.csv** / **sider_info.csv**  
  Detailed information on siderophore BGCs.

## Scripts

- **number_distribution.R**  
  Analyzes the copy number distribution of PBP2 genes across all strain genomes.  
  The pipeline for identifying PBP2 genes in each genome:  
  (1) Extract all proteins with sequence lengths between 200–500 amino acids.  
  (2) Filter by HMM search using **PF01497**.  
  (Corresponds to **Figure 4a**.)

- **position_distribution.R**  
  Uses siderophore BGCs as reference points to calculate the genomic distance distribution between all PBP2 genes and siderophore BGCs.  
  Also computes the distance distribution of the **nearest** PBP2 to each siderophore BGC.  
  (Corresponds to **Figure 4b**, left two panels.)

- **self_position_distribution.R**  
  Calculates the distance between **self-rec** siderophore receptors and their corresponding siderophore BGCs,  
  as well as the proportion of self-rec genes that are also the **nearest** PBP2 to a siderophore BGC.  
  (Corresponds to **Figure 4b**, right two panels.)

- **Venn_plot.R**  
  Examines the overlap between strains containing siderophore BGCs and those containing self-rec receptors.  
  (Corresponds to **Figure 4c**.)

- **TF_analysis.R**  
  Investigates potential regulatory mechanisms between distant siderophore BGCs and self-rec receptors, focusing on TF binding site associations.  
  (Corresponds to **Figure 4d**.)
