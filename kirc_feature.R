#######################################################################
ANOVA
#######################################################################
kirc_stage_I_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s1_rawdata.csv")
kirc_stage_II_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s2_rawdata.csv")
kirc_stage_III_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_S3_rawdata.csv")
kirc_stage_IV_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s4_rawdata.csv")
result<-NULL
for (i in c(1:19615)){
  print(i)
  X <- as.numeric( c(kirc_stage_I_rawdata[i,5:ncol(kirc_stage_I_rawdata)],
                     kirc_stage_II_rawdata[i,5:ncol(kirc_stage_II_rawdata)],
                     kirc_stage_III_rawdata[i,5:ncol(kirc_stage_III_rawdata)],
                     kirc_stage_IV_rawdata[i,5:ncol(kirc_stage_IV_rawdata)]) )
  X <- X+1
  Y <- factor(rep(1:4, c(297,70,139,102)))
  df <- data.frame(X,Y)
  aov.df <- aov(X~Y, data=df)
  aov<-summary(aov.df)
  #if(aov[[1]]$`Pr(>F)`[1] <= 0.05)
  F_value<-aov[[1]]$`Pr(>F)`[1]
  Z<-cbind(kirc_stage_I_rawdata[i,1:3],F_value)
  result<-rbind(result,Z)
}
write.csv(result,"C:/Users/syuan/Desktop/Master/kirc/ANOVA_kirc.csv")
#ANOVA PC
X <- c(1:4)
Y <- c(1:4)
df <- data.frame(X,Y)
aov.df <- aov(X~Y, data=df)
aov<-summary(aov.df)
PC_result<-aov[[1]]$`Pr(>F)`[1]
#ANOVA NC
X <- c(1,1,1,1)
Y <- c(1:4)
df <- data.frame(X,Y)
aov.df <- aov(X~Y, data=df)
aov <-summary(aov.df)
NC_result<-aov[[1]]$`Pr(>F)`[1]
#ANOVA sample
X <- c(1,4,2,4)
Y <- c(1:4)
df <- data.frame(X,Y)
aov.df <- aov(X~Y, data=df)
aov<-summary(aov.df)
sample_result<-aov[[1]]$`Pr(>F)`[1]
######################################################################
Stage Correlation Coefficients
#######################################################################
kirc_stage_I_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s1_rawdata.csv")
kirc_stage_II_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s2_rawdata.csv")                                      
kirc_stage_III_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s3_rawdata.csv")                                      
kirc_stage_IV_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_s4_rawdata.csv")                                      
result<-NULL
for (i in c(1:19615)){
  print(i)
  kirc_I_mean<-mean( as.numeric(kirc_stage_I_rawdata[i,5:ncol(kirc_stage_I_rawdata)]))
  kirc_II_mean<-mean( as.numeric(kirc_stage_II_rawdata[i,5:ncol(kirc_stage_II_rawdata)]))
  kirc_III_mean<-mean( as.numeric(kirc_stage_III_rawdata[i,5:ncol(kirc_stage_III_rawdata)]))
  kirc_IV_mean<-mean( as.numeric(kirc_stage_IV_rawdata[i,5:ncol(kirc_stage_IV_rawdata)]))
  kirc_stage_mean<-cbind(kirc_stage_I_rawdata[i,1:3],kirc_I_mean,kirc_II_mean,kirc_III_mean,kirc_IV_mean)
  result<-rbind(result,kirc_stage_mean)
}
corresult<-NULL
for (i in c(1:19615)){
  print(i)
  kirc_cor<-cor( x=c(1:4),y=as.numeric(result[i,4:7]),method="spearman")
  corresult<-rbind(corresult,kirc_cor)}
kirc_cor<-cbind(result,corresult)
write.csv(kirc_cor,"C:/Users/syuan/Desktop/Master/kirc/kirc_cor.csv")
#######################################################################
m1/m0 FC
#######################################################################
kirc_ANOVA_correlation<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_ANOVA_correlation.csv")
kirc_m0_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m0_rawdata.csv")
kirc_m1_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv")
m0_mean<-as.data.frame( apply(kirc_m0_rawdata[,5:ncol(kirc_m0_rawdata)], FUN=mean,MARGIN=1) )
m1_mean<-as.data.frame( apply(kirc_m1_rawdata[,5:ncol(kirc_m1_rawdata)], FUN=mean,MARGIN=1) )
m1m0_FC<-m1_mean/m0_mean
m1m0_FC<-cbind(kirc_m0_rawdata$external_gene_name,m1m0_FC)
colnames(m1m0_FC)<-c("external_gene_name","m1m0_FC")
write.csv(m1m0_FC,"C:/Users/syuan/Desktop/Master/kirc/m1m0_FC.csv")
#######################################################################
T/N
#######################################################################
TCGAsetting <- GDCquery(project= "TCGA-KIRC", 
                        data.category = "Transcriptome Profiling", 
                        workflow.type = "HTSeq - FPKM",
                        access= "open",
                        experimental.strategy="RNA-Seq")
GDCdownload(TCGAsetting)
kirc_all <- GDCprepare(TCGAsetting, summarizedExperiment = F)
kirc_split<-as.data.frame(strsplit(colnames(kirc_all[2:612]),"-"))
kirc_split <- data.frame(lapply(kirc_split[4,], as.character), stringsAsFactors=FALSE)
kirc_normal<-kirc_all[,which(kirc_split=="11A")+1]
kirc_tumor<-kirc_all[,-c(1,(which(kirc_split=="11A")+1))]
kirc_all<-as.data.frame(kirc_all)
kirc_ensg<-as.character(kirc_all[,1])
split_data <- t(as.data.frame(strsplit(kirc_ensg,"\\.")))
rownames(split_data) <- NULL
ensembl_gene_id<-split_data[,1]
kirc_normal<-cbind(ensembl_gene_id, kirc_normal)
kirc_tumor<-cbind(ensembl_gene_id, kirc_tumor)
kirc_m1_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv")
ensg_external<-kirc_m1_rawdata[,2:3]
kirc_normal_rawdata<-merge(ensg_external,kirc_normal,by="ensembl_gene_id")
write.csv(kirc_normal_rawdata,"C:/Users/syuan/Desktop/Master/kirc/kirc_normal_rawdata.csv")
kirc_tumor_rawdata<-merge(ensg_external,kirc_tumor,by="ensembl_gene_id")
write.csv(kirc_tumor_rawdata,"C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_rawdata.csv")
colnames(kirc_normal_rawdata)
#######################################################################
?„²å­˜T/N
#######################################################################
kirc_normal_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_normal_rawdata.csv")
kirc_tumor_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_rawdata.csv")
kirc_TN_FC<-as.data.frame(apply(kirc_tumor_rawdata[,4:ncol(kirc_tumor_rawdata)],1,mean)/apply(kirc_normal_rawdata[,4:ncol(kirc_normal_rawdata)],1,mean))
kirc_TN_FC_merge<-cbind(kirc_tumor_rawdata$external_gene_name,kirc_TN_FC)
colnames(kirc_TN_FC_merge)<-c("external_gene_name","kirc_TN_FC")
write.csv(kirc_TN_FC_merge,"C:/Users/syuan/Desktop/Master/kirc/kirc_TN_FC_merge.csv")
#ANOVAå¯¦é??
X <- as.numeric( c(kirc_stage_I_rawdata[1,5:ncol(kirc_stage_I_rawdata)],
                   kirc_stage_II_rawdata[1,5:ncol(kirc_stage_II_rawdata)],
                   kirc_stage_III_rawdata[1,5:ncol(kirc_stage_III_rawdata)],
                   kirc_stage_IV_rawdata[1,5:ncol(kirc_stage_IV_rawdata)]) )
X <- X+1
Y <- factor(rep(1:4, c(297,70,139,102)))
df <- data.frame(X,Y)
aov.df <- aov(X~Y, data=df)
TukeyHSD(aov.df)
aov<-summary(aov.df)
#if(aov[[1]]$`Pr(>F)`[1] <= 0.05)
F_value<-aov[[1]]$`Pr(>F)`[1]
Z<-cbind(kirc_stage_I_rawdata[i,1:3],F_value)
result<-rbind(result,Z)
######################################################################
Marker Correlation Coefficients
#######################################################################
kirc_tumor_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_rawdata.csv")
marker_raw<-subset(kirc_tumor_rawdata,external_gene_name=="CDH2"|
                                          external_gene_name=="FN1"|
                                          external_gene_name=="VIM"|
                                          external_gene_name=="ZEB2"|
                                          external_gene_name=="FOXC2"|
                                          external_gene_name=="SNAI1"|
                                          external_gene_name=="SNAI2"|
                                          external_gene_name=="TWIST1"|
                                          external_gene_name=="TWIST2"|
                                          external_gene_name=="TGFBR1"|
                                          external_gene_name=="FGFR1"|
                                          external_gene_name=="CTNNB1"|
                                          external_gene_name=="MET"|
                                          external_gene_name=="EGFR"|
                     external_gene_name=="CDH1"|
                     external_gene_name=="BRMS1"|
                     external_gene_name=="MED23"|
                     external_gene_name=="CD82"|
                     external_gene_name=="KISS1"|
                     external_gene_name=="NME1"|
                     external_gene_name=="TP63"|
                     external_gene_name=="TXNIP"|
                     external_gene_name=="ARHGDIB"|
                     external_gene_name=="AKAP12"|
                     external_gene_name=="TIMP1"|
                     external_gene_name=="TIMP2"|
                     external_gene_name=="TIMP3"|
                     external_gene_name=="TIMP4"|
                     external_gene_name=="MAP2K4"|
                     external_gene_name=="DICER1")
marker_raw_numeric<-marker_raw[,4:ncol(marker_raw)]
result<-NULL
result2<-NULL
fun<-function(a){
    for (i in c(1:19615)){
r<-cor(x=as.numeric(a),y=as.numeric(kirc_tumor_rawdata[i,4:ncol(marker_raw)]),method="spearman")
result<-c(result,r) }
  result2<-rbind(result2,result)
  }
final_result<-apply( marker_raw_numeric,1,fun)
kirc_marker_coef<-final_result
colnames(kirc_marker_coef)<-marker_raw$external_gene_name
rownames(kirc_marker_coef)<-kirc_tumor_rawdata$external_gene_name
write.csv(kirc_marker_coef,"C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
#######################################################################
­«ºâT_high/normal
20200602
#######################################################################
kirc_tumor_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_order.csv")
kirc_tumor_order<-kirc_tumor_order[-seq(from=1, to=39229, by=2),]
kirc_tumor_order<-kirc_tumor_order[,-1]
kirc_tumor_order<-apply(kirc_tumor_order,2,as.numeric)
class(kirc_tumor_order[1,1])
mean_result<-NULL
y<-function(y){
  data<-y
  sd_result<-NULL
for(i in c(3:539)){
sd<-sd(data[1:i])
sd_result<-c(sd_result,sd)}
high_mean<-mean(as.numeric(data[1:which.max(sd_result)]))
mean_result<-c(mean_result,high_mean)
}
sd_result<-apply(kirc_tumor_order,1,y)
kirc_normal_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_normal_rawdata.csv")
normal_mean<-apply(kirc_normal_rawdata[,4:75],1,mean)
THN_FC<-cbind(as.character(kirc_normal_rawdata$external_gene_name),sd_result/normal_mean)
write.csv(THN_FC,"C:/Users/syuan/Desktop/Master/kirc/kirc_THN_FC.csv")
#######################################################################
¦X¨Öfeature
#######################################################################
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
THN_FC<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_THN_FC.csv")
kirc_surv_harz<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_surv_harz.csv")
kirc_merge_all<-cbind(kirc_merge_all,THN_FC[,3],kirc_surv_harz[,4])
write.csv(kirc_merge_all,"C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
