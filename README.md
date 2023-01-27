# Marker-assisted selection in a Pacific oyster population for an antiviral QTL conferring increased survival to OsHV-1 mortality events in Tomales Bay
This repository contains the data and code used in a [paper](https://doi.org/10.1016/j.aquaculture.2023.739291) in the journal Aquaculture.

Repository overview:
* SNP_pipeline.sh - Code to call and impute SNPs from genotyping-by-sequencing (GBS) data
* phenotypes.csv - Field trial survival phenotypes
* GWAS.R - Code to run GWAS
* A.csv - Additive pedigree relationship matrix of the MBP oyster families
* LMM.R - Code to run the linear mixed model (LMM) to obtain estimated breeding values (EBVs)
* genes_dCq.csv - Delta Cq values of antiviral genes and Chr8:9719736 SNP genotypes for oyster samples
* RT_qPCR.R - Code to analyze the gene expression data
