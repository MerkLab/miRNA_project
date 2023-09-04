# miRNA_project
Pan-cancer functional annotation of miRNAs

This repository provides all codes used to analyze CRISPR-Cas9 screens using the lentiG-miR library in 47 human cancer cell lines.

Code includes low-level QC, correction for gene-independent cell responses to CRISPR-Cas9 targeting (CRISPRcleanR),  funcitonal validation of screen performance, and cell line-level estimation of gene essentiality using BAGEL2 and MAGeCK-RRA. Furthermore, miRNAs with a significant effect in the majority of or all cancer cell lines investigated were determined using the 90th percentile mehtod (common essential genes) and the ADaM method (core fitness genes).
