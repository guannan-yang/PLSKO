This folder provides all reproducible codes for our article:

# PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection
biorXiv: https://www.biorxiv.org/content/10.1101/2024.08.06.606935v1

**Abstract**: The knockoff framework, combined with variable selection procedure, controls the false discovery rate (FDR) without the need for calculating p−values. Hence, it presents an attractive alternative to differential expression analysis of high-throughput biological data. However, current knockoff variable generators make strong assumptions or insufficient approximations that lead to FDR inflation when applied to biological data.

We propose Partial Least Squares Knockoff (PLSKO), an efficient and assumption-free knockoff generator that is robust to varying types of biological omics data. We compare PLSKO with a wide range of existing methods. In simulation studies, we show that PLSKO is the only method that controls FDR with sufficient statistical power in complex non-linear cases. In semi-simulation studies based on real data, we show that PLSKO generates valid knockoff variables for different types of biological data, including RNA-seq, proteomics, metabolomics and microbiome. In preeclampsia multi-omics case studies, we combined PLSKO with Aggregation Knockoff to address the randomness of knockoffs and improve power, and show that our method is able to select variables that are biologically relevant.

# Other resources
We also include the implementation of other knockoff-generating methods including IPAD (Fan et al., 2020), Sequential Knockoff (Kormaksson et al., 2021), Minimum variance-based reconstructability (MVR) knockoff and Maximum Entropy (ME) knockoff (Spector and Janson, 2022; implemented in Julia by [Chu](https://github.com/biona001/knockoffsr)). 
