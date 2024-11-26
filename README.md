# Honours2024CodeKD
All scripts used in "Investigating gene expression changes caused by ZMYND8 knock-out in human brain organoids" by Kaya Dahlke in the Voineagu lab

## rnaseqOutputPreprocessing.R
Script taking in nf-core/RNA-seq output and filtering to create an RDS file with main input for future scripts. 

## Differential expression analyses:
### DEAcrossTimeWT.R & DESeparateTimepoints.R
![image](https://github.com/user-attachments/assets/aabc62ed-a6f9-4f1d-a966-17a1af694bda)
Two DESeq2 modules focusing on WT expression over time (DEAcrossTimeWT.R) and expression changes between WT and KO at each timepoint (DESeparateTimepoints.R)

## WTVsPublishedOrganoids.R
Comparison of DEAcrossTimeWT.R results to Pasca data (https://doi.org/10.1038/nmeth.3415 sup table 1), as well as validation experiments not included in thesis comparing to 2019 Sloan et al. (https://doi.org/10.1038%2Fs41596-018-0032-7 table 1) and 2023 Mulder et al. (https://doi.org/10.1186/s13287-023-03302-x sup file 5)

# WGCNAFromTPM.R
Weighted gene co-expression network analysis of all WT and KO samples.

## ZMYND8IsoformsAcrossTime.R
Isoform analysis from WT iPSCs to WT D100 using FANTOM5 annotated isoforms as input in addition to GTEx.  

## expViewShinyApp folder
Contains two part shinyApp (split into server.R and ui.R) for visualisation of TPM data: https://ivlab-unsw.shinyapps.io/organoidExpressionTPM/

## expViewShinyApp folder
Contains single part shinyApp (app.R only) for visualisation of WGCNA and timepoint DE results: https://ivlab-unsw.shinyapps.io/organoidExpressionTPM/
