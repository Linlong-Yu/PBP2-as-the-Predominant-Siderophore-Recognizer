# README

## Input Files

- **genus_gbk_summary.csv**  
  Extracted from antiSMASH results, recording the total number of secondary metabolite regions predicted for each strain (i.e., the number of predicted regions).

- **genus_info.csv**  
  Contains genus information for each strain, as well as the complete list of strain names analyzed from RefSeq.

- **RAxML_result.GCF_data_refined.tre**  
  Phylogenetic tree computed using **PhyloPhlAn**, representing evolutionary distances among different strains.

- **summ_genus.csv**  
  Contains information about each siderophore BGC.

- **syn_info_5_genera.csv**  
  Core biosynthetic genes were extracted from siderophore BGCs (as annotated by antiSMASH) to represent each siderophore’s synthetase genes.  
  Sequence distances among these biosynthetic genes were calculated, followed by hierarchical clustering to group siderophores.  
  For siderophore groups with sizes greater than 3 (in *Rhodococcus* and *Corynebacterium*) or greater than 5 (in *Streptomyces*, *Staphylococcus*, and *Bacillus*), product annotations were assigned by comparing with the **MIBiG** database.

## Scripts

- **phylogenetic_tree_total_5_genera.R**  
  Plots the phylogenetic tree of 5,940 strains.

- **SecondMetabolity_siderophorenumber.R**  
  Calculates the proportion of siderophores and secondary metabolites within each genus.

- **knowncluster.R**  
  Analyzes the similarity of siderophore BGCs to known BGCs in the **MIBiG** database.

- **sider_relative_abundance.R**  
  Computes the relative abundance of siderophore BGCs across genera.
