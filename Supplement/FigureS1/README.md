# README

## Input Files

- **genus_dataname.mat**  
  Contains the IDs of siderophore BGCs for each genus.

- **genus_distm.mat**  
  Pairwise distance matrix among siderophore BGCs within each genus.  
  The corresponding BGC IDs are recorded in **genus_dataname.mat**.

- **summ_genus.csv**  
  Detailed information of siderophore BGCs.

## Scripts

- **syn_cluster.mlx**  
  Evaluates clustering results under different threshold values and performs hierarchical clustering of siderophore BGCs based on sequence distance.
