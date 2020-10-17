library(TCGAbiolinks)
library("BiocManager")
library("biomaRt")
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#######################################################################
TCGAbiolinks:::getProjectSummary("TCGA-KIRC")
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Copy Number Variation",
                  data.type="Copy Number Segment")
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Copy Number Variation",
                  data.type="Masked Copy Number Segment")
GDCdownload(query,directory="C:/Users/syuan/Desktop/Master/kirc/Bioportal")
kirc_CNV <- GDCprepare(query,directory="C:/Users/syuan/Desktop/Master/kirc/Bioportal")
#######################################################################
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "DNA Methylation")
GDCdownload(query,directory="C:/Users/syuan/Desktop/Master/kirc/Bioportal")
kirc_Methylation <- GDCprepare(query,directory="C:/Users/syuan/Desktop/Master/kirc/Bioportal")
#######################################################################
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Simple Nucleotide Variation")
GDCdownload(query,directory="C:/Users/syuan/Desktop/Master/kirc/Bioportal")
kirc_Methylation <- GDCprepare(query,directory="C:/Users/syuan/Desktop/Master/kirc/Bioportal")