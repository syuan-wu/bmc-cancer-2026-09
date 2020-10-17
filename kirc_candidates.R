library(rlist)
library(ROCR)
library(ggplot2)
library(dplyr)
##################################################################################
TCGA_data
##################################################################################
#19615
kirc_m1_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv",stringsAsFactors = F)
kirc_m0_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m0_rawdata.csv",stringsAsFactors = F)
kirc_m1<-as.data.frame(t(kirc_m1_rawdata[,-c(1:4)]))
colnames(kirc_m1)<-kirc_m1_rawdata$external_gene_name
kirc_m0<-as.data.frame(t(kirc_m0_rawdata[,-c(1:4)]))
colnames(kirc_m0)<-kirc_m0_rawdata$external_gene_name
tcga<-rbind( cbind(kirc_m1,response=rep(1,97)),cbind(kirc_m0,response=rep(0,481)))
#split data
EMT_top<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/predict_result_n36/EMT_top.csv",stringsAsFactors = F)           
tcga_397<-select(tcga, EMT_top$top)
tcga_397<-data.frame(tcga_397,tcga$response)
tcga_19214<-select(tcga, -c(EMT_top$top,"MATR3","TMSB15B", "PRAMEF7","PINX1"))  
##################################################################################
Combination of gene set and ROC of logistic model
##################################################################################
#19214 v.s. 397
#397
auc_result<-NULL
gene_result<-NULL
for(i in c(1:10000)){
n<-sample(397,20)
z<-paste( colnames(tcga_397)[n],collapse = '+' )
z<-paste0("response","~",z)
md<-glm(formula=z,data=tcga_397,family="binomial", na.action=na.exclude)
prediction<-predict(md,tcga_397)
pred <- prediction(prediction,tcga_397$response)
auc <- performance(pred, "auc")@y.values
auc_result<-c(auc_result,as.character(auc))
gene_result<-rbind(gene_result,paste(colnames(tcga_397)[n],collapse=","))
}
gene_auc<-data.frame(gene=gene_result,auc=as.numeric(as.character(auc_result)))
gene_auc<-arrange(gene_auc,desc(auc))
#19214
n1<-grep("\\.",colnames(tcga_19214))
colnames(tcga_19214)[n1]<-gsub("\\.","_",colnames(tcga_19214)[n1])
n2<-grep("-",colnames(tcga_19214))
colnames(tcga_19214)[n2]<-gsub("-","_",colnames(tcga_19214)[n2])
auc_result<-NULL
gene_result<-NULL
for(i in c(1:10000)){
  n<-sample(19214,20)
  z<-paste( colnames(tcga_19214)[n],collapse = '+' )
  z<-paste0("response","~",z)
  md<-glm(formula=z,data=tcga_19214,family="binomial", na.action=na.exclude)
  prediction<-predict(md,tcga_19214)
  pred <- prediction(prediction,tcga_19214$response)
  auc <- performance(pred, "auc")@y.values
  auc_result<-c(auc_result,as.character(auc))
  gene_result<-rbind(gene_result,paste(colnames(tcga_19214)[n],collapse=","))
}
gene_auc_19214<-data.frame(gene=gene_result,auc=as.numeric(as.character(auc_result)))
gene_auc_19214<-arrange(gene_auc_19214,desc(auc))
write.csv(gene_auc_19214,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/gene_auc_19214.csv")
#397 v.s. 19615
gene_auc<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/gene_auc.csv")
gene_auc_19214<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/gene_auc_19214.csv")
df_d<-rbind( cbind(gene_auc,group=rep("Candidate (397)",10000)),cbind(gene_auc_19214,group=rep("Others (19214)",10000)))
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/density_TCGA.png",width = 800, height = 500, units = "px", pointsize = 36,
    bg = "white")
ggplot(df_d, aes(x=auc, fill=group)) +  ggtitle("Randomly Selection of Combination of TCGA")+geom_density(alpha=0.4)+theme_bw()+theme(text=element_text(size=20))
dev.off()
##################################################################################
gene_auc<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/gene_auc_20.csv",stringsAsFactors = F)
N<-length(which(gene_auc$auc>0.7))
top_auc<-NULL
for (i in c(1:N)){
top_auc<-c(top_auc,unlist(strsplit(as.character(gene_auc$gene[i]), ",")) )
}
top_auc<-as.data.frame(table(top_auc))
top_auc<-arrange(top_auc,desc(Freq))
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/top_auc_20_0.7.png",width = 800, height = 500, units = "px", pointsize = 24,
    bg = "white")
ggplot(top_auc, aes( x=c(1:397),y=Freq ))+ geom_bar(stat = "identity")+theme_bw()
dev.off()
#frequency挑前?個當candidate
auc_result<-NULL
for (i in c(1:nrow(top_auc))){
  z<-paste( top_auc$top_auc[1:i],collapse = '+' )
  z<-paste0("response","~",z)
  #z
  md<-glm(formula=z,df_n,family="binomial", na.action=na.exclude)
  #summary(md)
  #coef(md)
  prediction<-predict(md,df_n)
  pred <- prediction(prediction,df_n$response)
  auc <- performance(pred, "auc")@y.values
  auc_result<-c(auc_result,auc)
}
auc_result<-as.data.frame(list.rbind(auc_result))
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/fre_auc_20.png",width = 800, height = 500, units = "px", pointsize = 24,
    bg = "white")
ggplot(auc_result, aes( x=c(1:397),y=V1 ))+ geom_point()+theme_bw()
dev.off()
auc_result<-as.data.frame(list.rbind(auc_result))
#####fre random
#frequency挑前?個當candidate
auc_r<-NULL
top_auc_r<-top_auc[sample(1:397),]
for (i in c(1:nrow(top_auc_r))){
  z<-paste( top_auc_r$top_auc[1:i],collapse = '+' )
  z<-paste0("response","~",z)
  #z
  md<-glm(formula=z,df_n,family="binomial", na.action=na.exclude)
  #summary(md)
  #coef(md)
  prediction<-predict(md,df_n)
  pred <- prediction(prediction,df_n$response)
  auc <- performance(pred, "auc")@y.values
  auc_r<-c(auc_r,auc)
}
auc_r<-as.data.frame(list.rbind(auc_r))
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/fre_auc_20_random.png",width = 800, height = 500, units = "px", pointsize = 24,
    bg = "white")
ggplot(auc_r, aes( x=c(1:397),y=V1 ))+ geom_point()+theme_bw()
dev.off()
#取k第i個組合雙重for loop
N2<-length(which(top_auc$Freq>150))
auc_total<-list()
n_result<-list()
for (k in c(1:N2)){
  print(k)
  y<-as.data.frame(combn(N2, k))
  auc_result<-list()
  for(i in c(1:ncol(y)) ){
    z<-paste(top_auc$top_auc[y[,i]],collapse = '+' )
    z<-paste0("response","~",z)
    #z
    md<-glm(formula=z,df_n,family="binomial", na.action=na.exclude)
    #summary(md)
    #coef(md)
    prediction<-predict(md,df_n)
    pred <- prediction(prediction,df_n$response)
    auc <- performance(pred, "auc")@y.values
    auc_result<-c(auc_result,auc)
  }
  auc_total<-c(auc_total,list(auc_result))
  n_result<-c(n_result,y)
}
#best combination
auc_c<-NULL
for (i in c(1:N2)){
  auc_c<-c(auc_c,auc_total[[i]])
}
auc_c<-list.rbind(auc_c)
write.csv(auc_c,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/auc_c.csv")
write.csv(auc_cc,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/auc_cc.csv")

auc_cc<-data.frame(num=c(1:nrow(auc_c)),auc=as.numeric(auc_c[,1]))
q<-function(q){c( paste0(q,collapse = ","),as.character(paste0(top_auc$top_auc[q],collapse = "+")) )}
result_f<-list.rbind( lapply(n_result,q) )
auc_c<-data.frame(cbind(auc_c,result_f))
auc_c<-arrange(auc_c,desc(X1))
which(auc_cc$auc==auc_c[1,1])
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/auc_from_N2.png",width = 1200, height = 500, units = "px", pointsize = 24,
    bg = "white")
ggplot(auc_cc, aes( x=num,y=auc )) + geom_point() +coord_cartesian(ylim = c(0.5, 0.8))+ geom_point() +theme_bw()
dev.off()
#df_r<-c( df[1,1], N, N2, paste( as.character(top[1:N2,1]),collapse = ',' ),auc_cc[nrow(auc_cc),2], as.character(auc_c[1,3]),as.character(auc_c[1,1])  )
#result<-cbind(result, c(df[1,1], N, N2, paste( as.character(top[1:N2,1]),collapse = ',' ), as.character(auc_c[1,1]),as.character(auc_c[1,2]),as.character(auc_c[1,3])))
write.csv(gene_auc,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/gene_auc.csv")
##########################################################
  AUC distribution 2020.10.14
##########################################################

##########################################################
Candidate genes-111
##########################################################
gsea_n36<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/predict_result_n36/gsea_n36.csv",stringsAsFactors = F)           
EMT_top<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/predict_result_n36/EMT_top.csv",stringsAsFactors = F)           
EMT_top<-EMT_top[1:111,]
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv",stringsAsFactors = F)
x<-data.frame(external_gene_name=as.character(EMT_top$top))
kirc_n_all<-merge(data.frame( external_gene_name=as.character(EMT_top$top) ),kirc_merge_all,by="external_gene_name")
#kirc_n_all<-subset(kirc_n_all,kirc_n_all$m1m0_FC>1&kirc_n_all$kirc_TN_FC>1&kirc_n_all$corresult>0.4)
kirc_m1_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv",stringsAsFactors = F)
kirc_m0_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m0_rawdata.csv",stringsAsFactors = F)
#n
kirc_m1<-t(merge(data.frame(external_gene_name=kirc_n_all$external_gene_name),kirc_m1_rawdata,by="external_gene_name"))
kirc_m0<-t(merge(data.frame(external_gene_name=kirc_n_all$external_gene_name),kirc_m0_rawdata,by="external_gene_name"))
x<-head( cbind(kirc_m1[5:nrow(kirc_m1),],rep(1,97)) )
df_n<-rbind( cbind(kirc_m1[5:nrow(kirc_m1),],rep(1,97)),cbind(kirc_m0[5:nrow(kirc_m0),],rep(0,481)))
colnames(df_n)<-c(as.character(kirc_n_all$external_gene_name),"response")
df_n<-data.frame(apply(df_n,2,as.numeric))
#######################################################
Combination of gene set and ROC of logistic model-111
#######################################################
#取k第i個組合雙重for loop
auc_total<-list()
n_result<-list()
for (k in c(1:111)){
  print(k)
  y<-as.data.frame(combn(111, k))
  auc_result<-list()
  for(i in c(1:length(y)) ){
    z<-paste( df_n[y[,i]],collapse = '+' )
    z<-paste0("response","~",z)
    #z
    md<-glm(formula=z,df_n,family="binomial", na.action=na.exclude)
    #summary(md)
    #coef(md)
    prediction<-predict(md,df_n)
    pred <- prediction(prediction,df_n$response)
    auc <- performance(pred, "auc")@y.values
    auc_result<-c(auc_result,auc)
  }
  auc_total<-c(auc_total,list(auc_result))
  n_result<-c(n_result,y)
}
#best combination
auc_c<-NULL
for (i in c(1:N2)){
  auc_c<-c(auc_c,auc_total[[i]])
}
auc_c<-list.rbind(auc_c)
auc_cc<-data.frame(num=c(1:nrow(auc_c)),auc=as.numeric(auc_c[,1]))
q<-function(q){c( paste0(q,collapse = ","),as.character(paste0(fre_top$top_auc[q],collapse = "+")) )}
result_f<-list.rbind( lapply(n_result,q) )
auc_c<-data.frame(cbind(auc_c,result_f))
auc_c<-arrange(auc_c,desc(X1))top_auc
which(auc_cc$auc==auc_c[1,1])
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/auc_from_N2.png",width = 1200, height = 500, units = "px", pointsize = 24,
    bg = "white")
ggplot(auc_cc, aes( x=num,y=auc )) + geom_point() +coord_cartesian(ylim = c(0.5, 0.75))+ geom_point() +theme_bw()
dev.off()
#df_r<-c( df[1,1], N, N2, paste( as.character(top[1:N2,1]),collapse = ',' ),auc_cc[nrow(auc_cc),2], as.character(auc_c[1,3]),as.character(auc_c[1,1])  )
#result<-cbind(result, c(df[1,1], N, N2, paste( as.character(top[1:N2,1]),collapse = ',' ), as.character(auc_c[1,1]),as.character(auc_c[1,2]),as.character(auc_c[1,3])))
write.csv(gene_auc,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/gene_auc.csv")