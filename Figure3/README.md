# README

## Input Files

- **PBP2_msa_889.mat**  
  Multiple sequence alignment (MSA) result containing 888 Bacillus receptor sequences (FatB, FpuA, FeuA, YfiY, YclQ, and YxeB) obtained through BLAST, plus one reference sequence (**FeuA**).

- **full_sec_struct.mat**  
  Secondary structure information of the reference FeuA protein (PDB ID: **2WHY**) corresponding to the MSA result.

- **GCF_fasta_concatenated.aln**  
  Phylogenetic tree computed using **PhyloPhlAn**, representing evolutionary relationships among different *Bacillus* strains.

## Scripts

- **Feature_sequence_MI.mlx**  
  Calculates the feature sites of PBP2 based on a **mutual information (MI)** approach.  
  Also visualizes the clustering of full-length PBP2 sequences and the identified feature sites.

- **bar_PBP2.R**  
  Generates a schematic representation of the secondary structure elements along the PBP2 sequence.

- **high_MI_complete_bacillibactin.pse**  
  PyMOL session file visualizing the positions of FeuA feature sites with high mutual information scores.

- **heatmap_phylogeny.mlx**  
  Visualizes clustering patterns of full-length PBP2 sequences, feature sites, and their relationships to phylogenetic distances.
