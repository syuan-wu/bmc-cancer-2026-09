library(survival)
library(survminer)
##########################################################
前處理下載好的clinical data
##########################################################
clinical <- t(read.table("C:/Users/syuan/Desktop/Master/kirc/KIRC.merged_only_clinical_clin_format.txt",header=T, row.names=1, sep='\t'))
clinical<-as.data.frame(clinical)
?toupper
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
?grep
ind_keep <- grep('days_to_new_tumor_event_after_initial_treatment',colnames(clinical))
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}
# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}
# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum(is.na(fl[i,]) ) < dim(fl)[2] ){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}
# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
clinical$days_to_new_tumor_event_after_initial_treatment
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')
##########################################################
# create vector with time to new tumor containing data to censor for new_tumor
all_clin$new_time <- c()
for (i in 1:length(as.numeric(as.character(all_clin$new_tumor_days)))){
  all_clin$new_time[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$new_tumor_days))[i]),
                                   as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$new_tumor_days))[i])
}
# create vector time to death containing values to censor for death
all_clin$new_death <- c()
for (i in 1:length(as.numeric(as.character(all_clin$death_days)))){
  all_clin$new_death[i] <- ifelse ( is.na(as.numeric(as.character(all_clin$death_days))[i]),
                                    as.numeric(as.character(all_clin$followUp_days))[i],as.numeric(as.character(all_clin$death_days))[i])
}
##########################################################
# create vector for death censoring
table(clinical$patient.vital_status)
# alive dead
# 372   161
all_clin$death_event <- ifelse(clinical$patient.vital_status == 'alive', 0,1)
#finally add row.names to clinical
rownames(all_clin) <- clinical$IDs
write.csv(all_clin,"C:/Users/syuan/Desktop/Master/Survival/Clinical/kirc_all_clin.csv")
######################################################
基因表現量從高排到低
######################################################
kirc_tumor_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_rawdata.csv")
case_id<-t(matrix(colnames(kirc_tumor_rawdata[4:ncol(kirc_tumor_rawdata)])))
result<-NULL
for(i in c(1:19615)){
print(i)
x<-t(as.matrix(as.numeric(kirc_tumor_rawdata[i,4:ncol(kirc_tumor_rawdata)])))
gene<-rbind(case_id,x)
colnames(gene)<-NULL
y<-gene[,order(as.numeric(gene[2,]),decreasing = T)]
result<-rbind(result,y)
}
kirc_tumor_order<-as.data.frame(result)    
write.csv(kirc_tumor_order,"C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_order.csv")
######################################################
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("TCGAutils")
#devtools::install_github("renkun-ken/rlist")
#install.packages("rlist")
library(rlist)
kirc_tumor_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_order.csv")
kirc_tumor_order<-kirc_tumor_order[-seq(2,19615*2,by=2),-1]
kirc_high<-kirc_tumor_order[,1:180]
kirc_low<-kirc_tumor_order[,360:539]
kirc_high<-as.data.frame(gsub("[.]", "-", as.matrix(kirc_high)))
kirc_low<-as.data.frame(gsub("[.]", "-", as.matrix(kirc_low)))
kirc_high <- data.frame(lapply(kirc_high, as.character), stringsAsFactors=FALSE)
kirc_low <- data.frame(lapply(kirc_low, as.character), stringsAsFactors=FALSE)
x<-function(x){
  b<-paste(a[1:3,x],collapse="-")
}
result<-data.frame(matrix(NA,180,1))
a<-NULL
for (i in c(1:19615)){
print(i)
a<-strsplit(as.character(kirc_high[i,1:180]),"-")
a<-list.cbind(a)
a<-a[1:3,1:180]
kirc_high_id<-as.data.frame(apply(as.data.frame(c(1:180)),1,x))
result<-cbind(result,kirc_high_id)
}
kirc_high_barcode<-result[,-1]
#
result<-data.frame(matrix(NA,180,1))
a<-NULL
for (i in c(1:19615)){
  print(i)
  a<-strsplit(as.character(kirc_low[i,1:180]),"-")
  a<-list.cbind(a)
  a<-a[1:3,1:180]
  kirc_low_id<-as.data.frame(apply(as.data.frame(c(1:180)),1,x))
  result<-cbind(result,kirc_low_id)
}
kirc_low_barcode<-result[,-1]
#
write.csv(kirc_high_barcode,"C:/Users/syuan/Desktop/Master/kirc/kirc_high_barcode.csv")
write.csv(kirc_low_barcode,"C:/Users/syuan/Desktop/Master/kirc/kirc_low_barcode.csv")
#######################################################################
算出所有survival curve p value
#######################################################################
kirc_high_barcode<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_high_barcode.csv")
kirc_low_barcode<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_low_barcode.csv")
kirc_all_clin<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_all_clin.csv")
colnames(kirc_all_clin)[1]<-"bcr_patient_barcode"
kirc_high_barcode<-kirc_high_barcode[,-1]
kirc_low_barcode<-kirc_low_barcode[,-1]
group<-as.data.frame(c(rep("H",180),rep("L",180)))
result<-NULL
for (i in c(1:19615)){
print(i)
kirc_sur<-cbind(c(as.character(kirc_high_barcode[,i]),as.character(kirc_low_barcode[,i])),group)
colnames(kirc_sur)<-c("bcr_patient_barcode","group")
kirc_sur<-merge(kirc_sur,kirc_all_clin,by="bcr_patient_barcode")
# run survival analysis
kirc_sur<-kirc_sur[order(kirc_sur$group),]
fit <- survfit(Surv(as.numeric(as.character(kirc_sur$new_death)),kirc_sur$death_event) ~ group, data = kirc_sur)
test<-surv_pvalue(fit, kirc_sur)$pval
result<-c(result,test)
}
pval<-result
#第一個引數是time，生存時間，對於右截尾資料，這是follow up time
#第二個引數是event, 即the status indicator, normally 0=alive,1=dead
#model <- survfit(Sur_obj~1) 如果用所有的資料，不進行分組，則後面引數用1
kirc_surv_pval<-cbind(kirc_tumor_rawdata[2:3],pval)
write.csv(kirc_surv_pval,"C:/Users/syuan/Desktop/Master/kirc/kirc_surv_pval.csv")
#write.csv(kirc_sur,"C:/Users/syuan/Desktop/Master/kirc/kirc_sur.csv")
#######################################################################
合併feature
#######################################################################
ANOVA_kirc<-read.csv("C:/Users/syuan/Desktop/Master/kirc/ANOVA_kirc.csv")
kirc_cor<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_cor.csv")
kirc_surv_pval<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_surv_pval.csv")
kirc_TN_FC_merge<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_TN_FC_merge.csv")
m1m0_FC<-read.csv("C:/Users/syuan/Desktop/Master/kirc/m1m0_FC.csv")
kirc_merge_all<-cbind( kirc_cor[,4:9],ANOVA_kirc[,4:5],m1m0_FC[,2:3],kirc_TN_FC_merge[,2:3],kirc_surv_pval[,3:4])
kirc_merge_all<-kirc_merge_all[,-c(7,9,11,13)]
write.csv(kirc_merge_all,"C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
#######################################################################
輸出marker survival curve
#######################################################################
kirc_high_barcode<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_high_barcode.csv")
kirc_low_barcode<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_low_barcode.csv")
kirc_positive_marker<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_positive_marker.csv")
kirc_negative_marker<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_negative_marker.csv")
kirc_sur<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_sur.csv")
kirc_all_clin<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_all_clin.csv")
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_result.csv")
group<-as.data.frame(c(rep("H",180),rep("L",180)))
#survival positive marker
colnames(kirc_all_clin)[1]<-"bcr_patient_barcode"
kirc_high_barcode<-kirc_high_barcode[,-1]
kirc_low_barcode<-kirc_low_barcode[,-1]
group<-as.data.frame(c(rep("H",180),rep("L",180)))
#positive marker
x<-as.numeric(kirc_positive_marker$X)
plot <- list()
for ( i in c(1:14) ){
print(i)
kirc_sur<- cbind( c( as.character(kirc_high_barcode[,x[i]]),as.character(kirc_low_barcode[,x[i]]) ),group  )
colnames(kirc_sur)<-c("bcr_patient_barcode","group")
kirc_sur<-merge(kirc_sur,kirc_all_clin,by="bcr_patient_barcode")
kirc_sur<-kirc_sur[order(kirc_sur$group),]
# run survival analysis
fit <- survfit(Surv(as.numeric(as.character(kirc_sur$new_death)),kirc_sur$death_event) ~ group, data = kirc_sur)
name<-as.character(kirc_positive_marker$gene_name[i] )
#繪製marker的生存曲線
plot[[i]]<-ggsurvplot(fit, title = name ,pval = TRUE,ylab="Survival Probability",xlab="Time (days)")
}
#Arrange multiple ggsurvplots and print the output
png("C:/Users/syuan/Desktop/Master/kirc/kirc_P_survplot.png",width = 1800, height = 1200, units = "px", pointsize = 24,
    bg = "white")
arrange_ggsurvplots(plot, print = TRUE,
                    ncol = 5, nrow = 3,title="KIRC marker")
dev.off()
#negative marker
x<-as.numeric(kirc_negative_marker$X)
plot <- list()
for ( i in c(1:16) ){
  print(i)
  kirc_sur<- cbind( c( as.character(kirc_high_barcode[,x[i]]),as.character(kirc_low_barcode[,x[i]]) ),group  )
  colnames(kirc_sur)<-c("bcr_patient_barcode","group")
  kirc_sur<-merge(kirc_sur,kirc_all_clin,by="bcr_patient_barcode")
  kirc_sur<-kirc_sur[order(kirc_sur$group),]
  # run survival analysis
  fit <- survfit(Surv(as.numeric(as.character(kirc_sur$new_death)),kirc_sur$death_event) ~ group, data = kirc_sur)
  name<-as.character(kirc_negative_marker$gene_name[i] )
  #繪製marker的生存曲線
  plot[[i]]<-ggsurvplot(fit, title = name ,pval = TRUE,ylab="Survival Probability",xlab="Time (days)")
}
#Arrange multiple ggsurvplots and print the output
png("C:/Users/syuan/Desktop/Master/kirc/kirc_N_survplot.png",width = 1800, height = 1200, units = "px", pointsize = 24,
    bg = "white")
arrange_ggsurvplots(plot, print = TRUE,
                    ncol = 4, nrow = 4,title="KIRC marker")
dev.off()
#######################################################################
hazard ratio
#######################################################################
kirc_tumor_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_tumor_rawdata.csv")
kirc_high_barcode<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_high_barcode.csv")
kirc_low_barcode<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_low_barcode.csv")
kirc_all_clin<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_all_clin.csv")
colnames(kirc_all_clin)[1]<-"bcr_patient_barcode"
kirc_high_barcode<-kirc_high_barcode[,-1]
kirc_low_barcode<-kirc_low_barcode[,-1]
group<-as.data.frame(c(rep("L",180),rep("H",180)))
result<-NULL
for (i in c(1:19615)){
kirc_sur<-cbind(c(as.character(kirc_low_barcode[,i]),as.character(kirc_high_barcode[,i])),group)
colnames(kirc_sur)<-c("bcr_patient_barcode","group")
kirc_sur<-merge(kirc_sur,kirc_all_clin,by="bcr_patient_barcode")
kirc_sur$group = relevel(kirc_sur$group, ref = "L")
fit.coxph <- coxph(Surv(new_death) ~ group, 
                   data = kirc_sur)
hr <- coef(summary(fit.coxph))[,2]
result<-c(result,hr)
}
kirc_surv_harz<-cbind(kirc_tumor_rawdata[2:3],result)
write.csv(kirc_surv_harz,"C:/Users/syuan/Desktop/Master/kirc/kirc_surv_harz.csv")
######################################################################
畫圖
######################################################################
which(kirc_tumor_rawdata$external_gene_name=="TP53")
#8021
kirc_sur<-cbind(c(as.character(kirc_low_barcode[,8021]),as.character(kirc_high_barcode[,8021])),group)
colnames(kirc_sur)<-c("bcr_patient_barcode","group")
kirc_sur<-merge(kirc_sur,kirc_all_clin,by="bcr_patient_barcode")
kirc_sur$group <- relevel(kirc_sur$group, ref = "L")
fit.coxph <- coxph(Surv(new_death) ~ group, 
                   data = kirc_sur)
hr <- coef(summary(fit.coxph))[,2]
png("C:/Users/syuan/Desktop/Master/kirc/kirc_coxplot_p53.png",width = 800, height = 400, units = "px", pointsize = 24,
    bg = "white")
ggforest(fit.coxph, data = kirc_sur,main="Hazard Ratio of TP53 expression", fontsize =0.9 )
dev.off()
head(result)

