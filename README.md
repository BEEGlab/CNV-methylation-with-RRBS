# Bioinformatic scripts used in "Epigenetic diversity of genes with copy number variations among natural populations of the three-spined stickleback"

## Description


This project aimed to assess DNA methylation of copy number variations (CNVs) in threespine sticklebacks using RRBS. It investigated the role of DNA methylation in the evolution and adaptive divergence of duplicated genes in threespine sticklebacks across salinity gradients. Despite the potential benefits of gene duplication for adaptation, dosage imbalance can destabilize cells. Our findings suggest that while DNA methylation is linked to young duplicate genes, it does not strongly interact with copy number variations (CNVs), challenging the hypothesis that methylation significantly contributes to dosage regulation following duplication.

This Git includes R scripts for analyzing genetic data, aiming to understand methylation patterns across CNVs.

The scripts published here have been used in the manuscript **"Epigenetic diversity of genes with copy number variations among natural populations of the three-spined stickleback"** by Frédéric J. J. Chain, Britta S. Meyer, Melanie J. Heckwolf, Sören Franzenburg, Christophe Eizaguirre and Thorsten B.H. Reusch. We will add the DOI/link as soon as the publication is fully accepted.


## Features
- RRBS data preprocessing and cleaning [scripts](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/1.Summaryscript_RRBS_CNV.md).
- Genome data processcing and cleaning [scripts](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/0.Summaryscript_WGS_CNV.md).
- CNV identification [scripts](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/1.CNVs_IdentifyCandidatesBasedOnDepth.md) and [here](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/2a.CNVs_ManualModBasedOnDepth.md) and [more](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/2b.CNVs_ManualModBasedOnDepth_Genes.md) 
- CNV and Vst [scripts](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/3.CNVs_Vst.md)
- CNV and Methylation assessment [scripts](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/4.CNVs_Methylation_AllSites.md).
- Statistical analysis for methylation comparison, [here](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/5.DMS.md).
- Visualization of methylation patterns and CNV distribution.
- Some metadata ([Samples96.out](https://github.com/BEEGlab/CNV-methylation-with-RRBS/blob/main/Samples96.out))



## Data
Scripts are designed for standard RRBS/Genome fastq. 
The WGS raw reads are accessible from the BioProject [PRJNA749309](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA749309) and the RRBS data in the BioProjects [PRJNA591803](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA591803) and [PRJNA743784](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA743784). 


## License
Licensed under the MIT License.

## Authors and Acknowledgments
- Contributors: Chain and Meyer

We express our gratitude for the financial support provided to the BAMBI project (grant agreement number: call 2012-76) by BONUS; the joint Baltic Sea research program, funded by the European Union; the Federal Ministry of Education and Research in Germany (to C.E. and T.B.H.R., reference number 03F0680A). The research conducted by B.S.M. received support from the Excellence Cluster “The Future Ocean” (EXC 80), and F.J.J.C. was supported by the National Science Foundation (grant number 2144259).

## Citations
Epigenetic diversity of genes with copy number variations among natural populations of the three-spined stickleback
Frédéric J. J. Chain, Britta S. Meyer, Melanie J. Heckwolf, Sören Franzenburg, Christophe Eizaguirre, Thorsten B. H. Reusch
First published: 14 July 2024 https://doi.org/10.1111/eva.13753

