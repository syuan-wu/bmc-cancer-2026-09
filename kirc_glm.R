library(org.Hs.eg.db)
library(clusterProfiler)
library(ROCR)
library(dplyr)
#小鼠org.Mm.eg.db
library(DOSE)
#for enrichplot
library(msigdbr)
library(enrichplot)
#GSEA繪圖
library(pathview)
library(rlist)
#模型共線性檢查VIF(希望VIF值都在4以下)
library(car)
library(corrplot)
library(effects)
library(ggplot2)
library(ROCR)
library(pheatmap)
#######################################################################
self model-kirc_hc_2
#######################################################################
kirc_hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_1.csv")
data<-as.data.frame(kirc_hc_1[,3:8])
data[,c(2,4)] <- -log(data[,c(2,4)])
data<-scale(data)
colnames(data)[c(2,4,5,6)]<-c("-log(F_value)","-log(p-value)","kirc_THN_FC","kirc_Haza")
rownames(data)<-kirc_hc_1$kirc_merge_all.external_gene_name
data<-as.data.frame(data)
score<-3*data$m1m0_FC+3*data$kirc_THN_FC+2*data$`-log(p-value)`-data$kirc_Haza+(data$corresult+data$`-log(F_value)`)
data_score<-cbind(data,score)
data_score<-data_score[order(data_score$score,decreasing = T),] 
write.csv(data_score,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_221_selfmodel.csv",row.names =T)
#####################################################################
data_score<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_221_selfmodel.csv")
ranklist<-data_score$score
names(ranklist)<-data_score$X
head(ranklist)
sum(duplicated(ranklist))
sum(is.na(ranklist))
sum(is.nan(ranklist))
sum(is.infinite(ranklist))
#H
m_df <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_H.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:3, pvalue_table = T)
dev.off()
df_H<-as.data.frame(summary(gsea))
#C2
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_C2.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C2<-as.data.frame(summary(gsea))
#C5
m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_C5.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C5<-as.data.frame(summary(gsea))
#C1
m_df <- msigdbr(species = "Homo sapiens", category = "C1", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_C1.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1, pvalue_table = T)
dev.off()
df_C5<-as.data.frame(summary(gsea))
#C3
m_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_C3.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C3<-as.data.frame(summary(gsea))
#C4
m_df <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_C4.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C4<-as.data.frame(summary(gsea))
#C6
m_df <- msigdbr(species = "Homo sapiens", category = "C6", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea_C6.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C6<-as.data.frame(summary(gsea))
gsea<-rbind(df_H,df_C2,df_C5,df_C3,df_C4,df_C6)
gsea<-data.frame(gsea,dataset=c(rep("H",nrow(df_H)),rep("C2",nrow(df_C2)),rep("C5",nrow(df_C5)),rep("C3",nrow(df_C3)),rep("C4",nrow(df_C4)),rep("C6",nrow(df_C6)) ) )
write.csv(gsea,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_self_model/gsea.csv")
#######################################################################
logistic model
####################################################################### 
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
kirc_positive_marker<-subset(kirc_marker_coef,
                             kirc_marker_coef[,1]=="CDH2"|
                               kirc_marker_coef[,1]=="FN1"|
                               kirc_marker_coef[,1]=="VIM"|
                               kirc_marker_coef[,1]=="ZEB2"|
                               kirc_marker_coef[,1]=="FOXC2"|
                               kirc_marker_coef[,1]=="SNAI1"|
                               kirc_marker_coef[,1]=="SNAI2"|
                               kirc_marker_coef[,1]=="TWIST1"|
                               kirc_marker_coef[,1]=="TWIST2"|
                               kirc_marker_coef[,1]=="TGFBR1"|
                               kirc_marker_coef[,1]=="FGFR1"|
                               kirc_marker_coef[,1]=="CTNNB1"|
                               kirc_marker_coef[,1]=="MET"|
                               kirc_marker_coef[,1]=="EGFR")
kirc_negative_marker<-subset(kirc_marker_coef,
                             kirc_marker_coef[,1]=="CDH1"|
                               kirc_marker_coef[,1]=="BRMS1"|
                               kirc_marker_coef[,1]=="MED23"|
                               kirc_marker_coef[,1]=="CD82"|
                               kirc_marker_coef[,1]=="KISS1"|
                               kirc_marker_coef[,1]=="NME1"|
                               kirc_marker_coef[,1]=="TP63"|
                               kirc_marker_coef[,1]=="TXNIP"|
                               kirc_marker_coef[,1]=="ARHGDIB"|
                               kirc_marker_coef[,1]=="AKAP12"|
                               kirc_marker_coef[,1]=="TIMP1"|
                               kirc_marker_coef[,1]=="TIMP2"|
                               kirc_marker_coef[,1]=="TIMP3"|
                               kirc_marker_coef[,1]=="TIMP4"|
                               kirc_marker_coef[,1]=="MAP2K4"|
                               kirc_marker_coef[,1]=="DICER1")
premodel<-cbind(rbind(kirc_positive_marker,kirc_negative_marker),as.data.frame(c(rep(1,14),rep(0,16)) ))
colnames(premodel)[c(1,32)]<-c("gene","response")
###################################################################
#stepwise selection
full <- glm(formula=response~SNAI2+VIM+TIMP2+CDH1+MAP2K4+TP63+FGFR1+CD82+TIMP3+DICER1+TIMP1+MET+TGFBR1+ARHGDIB+MED23+
              FN1+TWIST1+SNAI1+AKAP12+EGFR+TIMP4+CTNNB1+ZEB2+KISS1+CDH2+BRMS1+FOXC2+TWIST2+NME1+TXNIP,data=premodel, family="binomial", na.action=na.exclude)
coef(full)
null <- glm(formula=response~1,data=premodel, family="binomial", na.action=na.exclude)
coef(null)
#從 null 開始或是從 full 開始都可以，只不過兩者的結果會不一樣(可以思考看看為什麼會這樣)
#從null
stepwise_model_n<-step(null, scope = list(upper=full), direction="both")
summary(stepwise_model_n)
coef(stepwise_model_n)
model_n30<-glm(formula = response ~ EGFR+TXNIP+TWIST1+AKAP12+CTNNB1 , 
               family = "binomial", data = premodel, na.action = na.exclude)
coef(model_n30)
#從full
stepwise_model_f<-step(full, scope = list(upper=full), direction="both")  
summary(stepwise_model_f)
coef(stepwise_model_f)
model_f30<-glm(formula = response ~CDH1+TIMP3+MED23+CTNNB1+CDH2+TWIST2+NME1 , family = "binomial", data = premodel, na.action = na.exclude)
#######################################################################
kirc-4-1-3277 prediction
#######################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_clustering_result.csv")
kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == "2")
predict_result_n30<-predict(model_n30,kirc_hc)
predict_result_f30<-predict(model_f30,kirc_hc)
full_result<-predict(full,kirc_hc)
kirc_hc_predict<-cbind(kirc_hc,predict_result_n30,predict_result_f30,full_result)
write.csv(kirc_hc_predict,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_hc_412_605.csv")
#######################################################################
依照response排序跑GSEA hc_41-n30
https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
#######################################################################
kirc_hc_412_605<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_hc_412_605.csv")
#
kirc_hc_412_605<-kirc_hc_412_605[order(kirc_hc_412_605$predict_result_n30,decreasing = T),]
ranklist<-kirc_hc_412_605$predict_result_n30
names(ranklist)<-kirc_hc_412_605$gene
head(ranklist)
sum(duplicated(ranklist))
sum(is.na(ranklist))
sum(is.nan(ranklist))
sum(is.infinite(ranklist))
#H
m_df <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_H_n.png",units="px",height=2000,width=3000,res=200)
gseaplot2(gsea,geneSetID=1, pvalue_table = F,base_size = 20)
dev.off()
df_H<-as.data.frame(summary(gsea))
#C2
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_C2.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C2<-as.data.frame(summary(gsea))
#C5
m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_C5.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C5<-as.data.frame(summary(gsea))
#C1
m_df <- msigdbr(species = "Homo sapiens", category = "C1", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_C1.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1, pvalue_table = T)
dev.off()
df_C1<-as.data.frame(summary(gsea))
#C3
m_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_C3.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C3<-as.data.frame(summary(gsea))
#C4
m_df <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_C4.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C4<-as.data.frame(summary(gsea))
#C6
m_df <- msigdbr(species = "Homo sapiens", category = "C6", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_C6.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C6<-as.data.frame(summary(gsea))
gsea_n30<-rbind(df_H,df_C2,df_C5,df_C3,df_C4,df_C6)
gsea_n30<-data.frame(gsea_n30,dataset=c(rep("H",nrow(df_H)),rep("C2",nrow(df_C2)),rep("C5",nrow(df_C5)),rep("C3",nrow(df_C3)),rep("C4",nrow(df_C4)),rep("C6",nrow(df_C6)) ) )
write.csv(gsea_n30,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30.csv")
#EMT
EMT_top<-data.frame(top=names(ranklist)[1:397])
write.csv(EMT_top,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/EMT_top.csv")
#NYCU seminar
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_n30_H_n.png",units="px",height=2000,width=3000,res=200)
gseaplot2(gsea,  title = "GSEA of 3277 gene candidates",geneSetID=1:5, pvalue_table = F,base_size = 20)
dev.off()
#######################################################################
依照response排序跑GSEA hc_41-f30
https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
#######################################################################
kirc_hc_412_605<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_hc_glm_predict.csv")
#
kirc_hc_412_605<-kirc_hc_412_605[order(kirc_hc_412_605$predict_result_f30,decreasing = T),]
ranklist<-kirc_hc_412_605$predict_result_f30
names(ranklist)<-kirc_hc_412_605$gene
head(ranklist)
sum(duplicated(ranklist))
sum(is.na(ranklist))
sum(is.nan(ranklist))
sum(is.infinite(ranklist))

#H
m_df <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_H.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:2, pvalue_table = T)
dev.off()
df_H<-as.data.frame(summary(gsea))
#C2
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_C2.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C2<-as.data.frame(summary(gsea))
#C5
m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_C5.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C5<-as.data.frame(summary(gsea))
#C1
m_df <- msigdbr(species = "Homo sapiens", category = "C1", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_C1.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1, pvalue_table = T)
dev.off()
df_C1<-as.data.frame(summary(gsea))
#C3
m_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_C3.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C3<-as.data.frame(summary(gsea))
#C4
m_df <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_C4.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C4<-as.data.frame(summary(gsea))
#C6
m_df <- msigdbr(species = "Homo sapiens", category = "C6", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30_C6.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C6<-as.data.frame(summary(gsea))
gsea_f30<-rbind(df_H,df_C2,df_C5,df_C3,df_C4,df_C6)
gsea_f30<-data.frame(gsea_f30,dataset=c(rep("H",nrow(df_H)),rep("C2",nrow(df_C2)),rep("C5",nrow(df_C5)),rep("C3",nrow(df_C3)),rep("C4",nrow(df_C4)),rep("C6",nrow(df_C6)) ) )
write.csv(gsea_f30,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_f30.csv")
#######################################################
#######################################################################
依照response排序跑GSEA hc_41-full
https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
#######################################################################
kirc_hc_412_605<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_hc_glm_predict.csv")
#
kirc_hc_412_605<-kirc_hc_412_605[order(kirc_hc_412_605$full_result,decreasing = T),]
ranklist<-kirc_hc_412_605$full_result
names(ranklist)<-kirc_hc_412_605$gene
head(ranklist)
sum(duplicated(ranklist))
sum(is.na(ranklist))
sum(is.nan(ranklist))
sum(is.infinite(ranklist))

#H
m_df <- msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_H.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1, pvalue_table = T)
dev.off()
df_H<-as.data.frame(summary(gsea))
#C2
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_C2.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C2<-as.data.frame(summary(gsea))
#C5
m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_C5.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C5<-as.data.frame(summary(gsea))
#C1
m_df <- msigdbr(species = "Homo sapiens", category = "C1", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_C1.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1, pvalue_table = T)
dev.off()
df_C1<-as.data.frame(summary(gsea))
#C3
m_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_C3.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C3<-as.data.frame(summary(gsea))
#C4
m_df <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_C4.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C4<-as.data.frame(summary(gsea))
#C6
m_df <- msigdbr(species = "Homo sapiens", category = "C6", subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
gsea <- GSEA(ranklist,nPerm = 1000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full_C6.png",units="px",height=3000,width=4000,res=250)
gseaplot2(gsea,  geneSetID=1:5, pvalue_table = T)
dev.off()
df_C6<-as.data.frame(summary(gsea))
gsea_full<-rbind(df_H,df_C2,df_C5,df_C3,df_C4,df_C6)
gsea_full<-data.frame(gsea_full,dataset=c(rep("H",nrow(df_H)),rep("C2",nrow(df_C2)),rep("C5",nrow(df_C5)),rep("C3",nrow(df_C3)),rep("C4",nrow(df_C4)),rep("C6",nrow(df_C6)) ) )
write.csv(gsea_full,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_full.csv")
#######################################################
PCA-meta marker
#######################################################
z<-paste(colnames(premodel)[2:37],collapse = '+')
z<-paste0("~",z)
z
pca<-prcomp(~corresult+F_value+m1m0_FC+pval+kirc_THN_FC+kirc_Haza+SNAI2+VIM+TIMP2+CDH1+MAP2K4+TP63+FGFR1+CD82+TIMP3+DICER1+TIMP1+MET+TGFBR1+ARHGDIB+MED23+FN1+TWIST1+SNAI1+AKAP12+EGFR+TIMP4+CTNNB1+ZEB2+KISS1+CDH2+BRMS1+FOXC2+TWIST2+NME1+TXNIP, premodel, scale=T)
pca
plot(pca,
     type="line",
     main="")
abline(h=1, col="blue")
vars <- (pca$sdev)^2
# 從pca中取出標準差(pca$sdev)後再平方，計算variance(特徵值)
vars
# 計算每個主成分的解釋比例 = 各個主成分的特徵值/總特徵值
props <- vars / sum(vars)    
props
cumulative.props <- cumsum(props)# 累加前n個元素的值
cumulative.props
plot(cumulative.props)
# 特徵向量(原變數的線性組合)
pca$rotation
top3.pca.eigenvector <- pca$rotation[, 1:3]
top3.pca.eigenvector
first.pca <- top3.pca.eigenvector[, 1]   #  第一主成份
second.pca <- top3.pca.eigenvector[, 2]  #  第二主成份
third.pca <- top3.pca.eigenvector[, 3]   #  第三主成份
first.pca[order(first.pca, decreasing=FALSE)]  
# 使用dotchart，繪製主成份負荷圖
dotchart(first.pca[order(first.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC1",                      # 主標題
         xlab="Variable Loadings") 
#######################################################
SVM-meta marker
#######################################################
library("e1071")
auc_result<-NULL
for( i in c(1:100)){
  print(i)
n <- floor(0.8*nrow(premodel)) 
s <- sample(nrow(premodel), n)
train <- premodel[s, ] # 80%
test <- premodel[-s, ] # 20%
model <- svm(train[,2:37],train$response,type="C-classification")
predict_result<-predict(model,test[,2:37])
print(predict_result)
pred <- prediction(as.numeric(predict_result),test$response)
auc <- performance(pred, "auc")@y.values
auc_result<-c(auc_result,auc)
}
auc_result<-list.rbind(auc_result)
mean(auc_result)
#######################################################
Decision Tree-meta marker
#######################################################
#install.packages("rpart")
premodel$response<-as.factor(premodel$response)
library(rpart)
n <- floor(0.8*nrow(premodel)) 
s <- sample(nrow(premodel), n)
train <- premodel[s, ] # 80%
test <- premodel[-s, ] # 20%
z<-paste(colnames(premodel)[2:37],collapse = '+')
z
cart.model<- rpart(response~corresult+F_value+m1m0_FC+pval+kirc_THN_FC+kirc_Haza+SNAI2+VIM+TIMP2+CDH1+MAP2K4+TP63+FGFR1+CD82+TIMP3+DICER1+TIMP1+MET+TGFBR1+ARHGDIB+MED23+FN1+TWIST1+SNAI1+AKAP12+EGFR+TIMP4+CTNNB1+ZEB2+KISS1+CDH2+BRMS1+FOXC2+TWIST2+NME1+TXNIP,data=train)
cart.model
#install.packages("rpart.plot")
library(rpart.plot) 
model<-prp(cart.model,  # 模型
    faclen=0,           # 呈現的變數不要縮寫
    fallen.leaves=TRUE, # 讓樹枝以垂直方式呈現
    shadow.col="gray",  # 最下面的節點塗上陰影
    # number of correct classifications / number of observations in that node
    extra=2)  
predict_result <- predict(cart.model, newdata=test, type="class")
table(real=test$response, predict=predict_result)
pred <- prediction(as.numeric(predict_result),test$response)
auc <- performance(pred, "auc")@y.values
printcp(cart.model) # 先觀察未修剪的樹，CP欄位代表樹的成本複雜度參數
#######################################################################
看model_f36, model_n36的acc
#######################################################################
test_acc_n36<-NULL
for( i in c(1:10)){
  print(i)
  n <- floor(0.8*nrow(premodel)) 
  s <- sample(nrow(premodel), n)
  train <- premodel[s, ] # 80%
  test <- premodel[-s, ] # 20%
  full <- glm(formula=response~corresult+F_value+m1m0_FC+kirc_THN_FC+kirc_Haza+pval+SNAI2+VIM+TIMP2+CDH1+MAP2K4+TP63+FGFR1+CD82+TIMP3+DICER1+TIMP1+MET+TGFBR1+ARHGDIB+MED23+
                FN1+TWIST1+SNAI1+AKAP12+EGFR+TIMP4+CTNNB1+ZEB2+KISS1+CDH2+BRMS1+FOXC2+TWIST2+NME1+TXNIP,data=train, family="binomial", na.action=na.exclude)
  null <- glm(formula=response~1,data=train, family="binomial", na.action=na.exclude)
  stepwise_model_n<-step(null, scope = list(upper=full), direction="both")
  z<-paste(names(coef(stepwise_model_n))[2:length(names(coef(stepwise_model_n)))],collapse = '+')
  z<-paste0("response","~",z)
  model<-glm(formula = z, family = "binomial", data = train, na.action = na.exclude)
  predict_result<-predict(model,test)
  predict_result[which(predict_result>0)]<-1
  predict_result[which(predict_result<0)]<-0
  print(predict_result)
  test$response
  #class(predict_result)
  #class(train$response)
  pred <- prediction(predict_result,test$response)
  auc <- performance(pred, "auc")@y.values
  test_acc_n36<-c(test_acc_n36,auc)
}
test_acc_f36<-NULL
for( i in c(1:10)){
  print(i)
  n <- floor(0.8*nrow(premodel)) 
  s <- sample(nrow(premodel), n)
  train <- premodel[s, ] # 80%
  test <- premodel[-s, ] # 20%
  full <- glm(formula=response~corresult+F_value+m1m0_FC+kirc_THN_FC+kirc_Haza+pval+SNAI2+VIM+TIMP2+CDH1+MAP2K4+TP63+FGFR1+CD82+TIMP3+DICER1+TIMP1+MET+TGFBR1+ARHGDIB+MED23+
                FN1+TWIST1+SNAI1+AKAP12+EGFR+TIMP4+CTNNB1+ZEB2+KISS1+CDH2+BRMS1+FOXC2+TWIST2+NME1+TXNIP,data=train, family="binomial", na.action=na.exclude)
  null <- glm(formula=response~1,data=train, family="binomial", na.action=na.exclude)
  stepwise_model_n<-step(full, scope = list(upper=full), direction="both")
  z<-paste(names(coef(stepwise_model_n))[2:length(names(coef(stepwise_model_n)))],collapse = '+')
  z<-paste0("response","~",z)
  model<-glm(formula = z, family = "binomial", data = train, na.action = na.exclude)
  predict_result<-predict(model,test)
  predict_result[which(predict_result>0)]<-1
  predict_result[which(predict_result<0)]<-0
  print(predict_result)
  test$response
  #class(predict_result)
  #class(train$response)
  pred <- prediction(predict_result,test$response)
  auc <- performance(pred, "auc")@y.values
  test_acc_f36<-c(test_acc_f36,auc)
}
##############################
df_gsea<-NULL
for ( i in c("H","C2","C5") ){
  m_df <- msigdbr(species = "Homo sapiens", category = i, subcategory = NULL)%>%dplyr::select(gs_name, gene_symbol)
  gsea <- GSEA(ranklist,nPerm = 10000,minGSSize = 10, maxGSSize = 1000,TERM2GENE = m_df,pvalueCutoff = 0.05,pAdjustMethod ="none")
  z<-paste0("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/gsea_",i)
  dic<-paste0(z,".png")
  png(dic,units="px",height=3000,width=4000,res=250)
  gseaplot2(gsea,  geneSetID=1:3, pvalue_table = T)
  dev.off()
  df<-as.data.frame(summary(gsea))
  df<-cbind(df,dataset=rep(i,nrow(df)))
  df_gsea<-rbind(df_gsea,df)
}
################################