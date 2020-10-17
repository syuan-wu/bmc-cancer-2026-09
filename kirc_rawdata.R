library(TCGAbiolinks)
library("BiocManager")
library("biomaRt")
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#######################################################################
下載Clinical資料
#######################################################################
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type="patient" )
GDCdownload(query,directory="C:/Users/syuan/Desktop/Master/Clinical")
clinical <- GDCprepare(query,directory="C:/Users/syuan/Desktop/Master/Clinical")
write.csv(clinical[["clinical_patient_kirc"]],"C:/Users/syuan/Desktop/Master/kirc/clinical_patient_kirc.csv")
#######################################################################
建立m1,m0 patient NGS資料
#######################################################################
clinical_patient_kirc<-read.csv("C:/Users/syuan/Desktop/Master/kirc/clinical_patient_kirc.csv")
kirc_m1<-as.character(clinical_patient_kirc[which(clinical_patient_kirc$ajcc_metastasis_pathologic_pm=="M1"),3])
kirc_m0<-as.character(clinical_patient_kirc[which(clinical_patient_kirc$ajcc_metastasis_pathologic_pm=="M0"),3])
TCGAsetting <- GDCquery(project= "TCGA-KIRC", 
                        data.category = "Transcriptome Profiling", 
                        workflow.type = "HTSeq - FPKM",
                        access= "open",
                        experimental.strategy="RNA-Seq",
                        barcode =kirc_m0)
GDCdownload(TCGAsetting)
kirc_m0 <- GDCprepare(TCGAsetting, summarizedExperiment = F)
TCGAsetting <- GDCquery(project= "TCGA-KIRC", 
                        data.category = "Transcriptome Profiling", 
                        workflow.type = "HTSeq - FPKM",
                        access= "open",
                        experimental.strategy="RNA-Seq",
                        barcode =kirc_m1)
GDCdownload(TCGAsetting)
kirc_m1 <- GDCprepare(TCGAsetting, summarizedExperiment = F)
#######################################################################
ensg去小數點後的ensg_60483->protein_coding_gene
#######################################################################
kirc_m0<-as.data.frame(kirc_m0)
kirc_m0_split<-as.character(kirc_m0[,1])
split_data <- t(as.data.frame(strsplit(kirc_m0_split,"\\.")))
rownames(split_data) <- NULL
ensembl_gene_id<-split_data[,1]
gene_name <- getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype'),
                   filters = 'ensembl_gene_id', mart= ensembl, 
                   values = as.character(ensg_60483))
protein_coding_gene<-gene_name[which(gene_name$gene_biotype=="protein_coding"),]
#######################################################################
m0,m1 NGS資料儲存
#######################################################################
kirc_m0<-cbind(ensembl_gene_id,kirc_m0[2:ncol(kirc_m0)])
kirc_m1<-cbind(ensembl_gene_id,kirc_m1[2:ncol(kirc_m1)])
write.csv(merge( protein_coding_gene, kirc_m0,by = "ensembl_gene_id") ,file="C:/Users/syuan/Desktop/Master/kirc/kirc_m0_rawdata.csv")
write.csv(merge( protein_coding_gene, kirc_m1,by = "ensembl_gene_id") ,file="C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv")
#######################################################################
建立stage1-4 patient NGS資料
#######################################################################
clinical_patient_kirc<-read.csv("C:/Users/syuan/Desktop/Master/Clinical/clinical_patient_kirc.csv")
kirc_stage<-as.character(unique(clinical_patient_kirc$ajcc_pathologic_tumor_stage))
print(kirc_stage)
kirc_stage<-c(kirc_stage[3],kirc_stage[5],kirc_stage[4],kirc_stage[6])
stage_sep <- function(x) {
  stage_sep <- as.character(clinical_patient_kirc[which(clinical_patient_kirc$ajcc_pathologic_tumor_stage==x),3])
}
kirc_stage_list<-lapply(kirc_stage, stage_sep)
names(kirc_stage_list)<- kirc_stage
################################################################
TCGAsetting <- GDCquery( project= "TCGA-KIRC", 
                         data.category = "Transcriptome Profiling", 
                         workflow.type = "HTSeq - FPKM",
                         access= "open",
                         experimental.strategy="RNA-Seq",
                         barcode =kirc_stage_list[["Stage I"]] )
GDCdownload(TCGAsetting)
S1 <- GDCprepare(TCGAsetting, summarizedExperiment = F)
################################################################
TCGAsetting <- GDCquery( project= "TCGA-KIRC", 
                         data.category = "Transcriptome Profiling", 
                         workflow.type = "HTSeq - FPKM",
                         access= "open",
                         experimental.strategy="RNA-Seq",
                         barcode =kirc_stage_list[["Stage II"]] )
GDCdownload(TCGAsetting)
S2 <- GDCprepare(TCGAsetting, summarizedExperiment = F)
################################################################
TCGAsetting <- GDCquery( project= "TCGA-KIRC", 
                         data.category = "Transcriptome Profiling", 
                         workflow.type = "HTSeq - FPKM",
                         access= "open",
                         experimental.strategy="RNA-Seq",
                         barcode =kirc_stage_list[["Stage III"]] )
GDCdownload(TCGAsetting)
S3 <- GDCprepare(TCGAsetting, summarizedExperiment = F)
################################################################
TCGAsetting <- GDCquery( project= "TCGA-KIRC", 
                         data.category = "Transcriptome Profiling", 
                         workflow.type = "HTSeq - FPKM",
                         access= "open",
                         experimental.strategy="RNA-Seq",
                         barcode =kirc_stage_list[["Stage IV"]] )
GDCdownload(TCGAsetting)
S4 <- GDCprepare(TCGAsetting, summarizedExperiment = F)
#######################################################################
S1,S2,S3,S4 NGS資料儲存
#######################################################################
S1<-cbind(ensembl_gene_id,S1[2:ncol(S1)])
S2<-cbind(ensembl_gene_id,S2[2:ncol(S2)])
S3<-cbind(ensembl_gene_id,S3[2:ncol(S3)])
S4<-cbind(ensembl_gene_id,S4[2:ncol(S4)])
write.csv(merge( protein_coding_gene, S1 ,by = "ensembl_gene_id") ,file="C:/Users/syuan/Desktop/Master/kirc/kirc_s1_rawdata.csv")
write.csv(merge( protein_coding_gene, S2 ,by = "ensembl_gene_id") ,file="C:/Users/syuan/Desktop/Master/kirc/kirc_s2_rawdata.csv")
write.csv(merge( protein_coding_gene, S3 ,by = "ensembl_gene_id") ,file="C:/Users/syuan/Desktop/Master/kirc/kirc_s3_rawdata.csv")
write.csv(merge( protein_coding_gene, S4 ,by = "ensembl_gene_id") ,file="C:/Users/syuan/Desktop/Master/kirc/kirc_s4_rawdata.csv")
##############################################################
測試載TCGA其他資料
##############################################################
TCGAbiolinks:::getProjectSummary("TCGA-KIRC")
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Simple Nucleotide Variation",
                  data.type = "Annotated Somatic Mutation" 
                  #data.format = "BCR Biotab",
                  #file.type="patient"
)
GDCdownload(query,directory="C:/Users/syuan/Desktop/Master/kirc")
clinical <- GDCprepare(query,directory="C:/Users/syuan/Desktop/Master/Clinical")
getwd()
######
TCGAbiolinks:::getProjectSummary("TCGA-KIRC")
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Protein expression",
                  data.type = "Annotated Somatic Mutation" 
                  #data.format = "BCR Biotab",
                  #file.type="patient"
)
GDCdownload(query,directory="C:/Users/syuan/Desktop/Master/kirc")
clinical <- GDCprepare(query,directory="C:/Users/syuan/Desktop/Master/Clinical")
getwd()
Protein expression