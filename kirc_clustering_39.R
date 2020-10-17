library(ggplot2)
library(grid)
library(factoextra)
library(scatterplot3d)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)
library(enrichplot)
library(pathview)
library(ggpubr)
library(rlist)
library(pheatmap)
#######################################################################
feature data前處理
#######################################################################
#處理NA,NaN,Inf
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
data<-as.matrix(cbind(kirc_merge_all[,3:11],kirc_marker_coef[,2:31]))
sum(is.na(data))
sum(is.nan(data))
sum(is.infinite(data))
data[is.na(data)] <- 0
data[is.infinite(data)] <- 100
data<-as.data.frame(data)
#發現不一樣的ensg對到同個external_gene_name,所以用ensg當分群的rowname
#kirc_m1_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv")
#data<-cbind(kirc_m1_rawdata$external_gene_name,as.data.frame(data))
#dup_data<-data[data$`kirc_m1_rawdata$external_gene_name`=="MATR3"|data$`kirc_m1_rawdata$external_gene_name`=="PINX1"|data$`kirc_m1_rawdata$external_gene_name`=="PRAMEF7"|data$`kirc_m1_rawdata$external_gene_name`=="TMSB15B",]
#rownames(data)<-kirc_m1_rawdata$ensembl_gene_id
#data<-data[,-1]
data_scale<-scale(data)
#準備hcluster分群
#######################################################################
HC決定階層式分群最佳群聚數
可使用 Elbow 或 Silhouette method
#######################################################################
#install.packages("factoextra")
library(factoextra)
#Elbow Method
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC")
# Avg. Silhouett
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC") 
#######################################################################
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc, hang = -2, cex = 0.6)
hc_3<-cutree(hc,k=3)
#HC實驗
result<-NULL
for (i in c(1:10)){
  x<-cutree(hc,k=i)
  result<-cbind(result,x)
  colnames(result)[i]<-paste("hc",as.character(i),sep="_")
}
write.csv(result,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_cutree.csv")

png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC cluster Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
#######################################################################
HC 分群結果合併及summary
#######################################################################
#結果合併
kirc_clustering_result<-cbind(kirc_merge_all$external_gene_name,data,hc_3)
write.csv(kirc_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_result_39.csv")
#結果summary
  x<-c(length(which(kirc_clustering_result$hc_3==1)),
       length(which(kirc_clustering_result$hc_3==2)),
       length(which(kirc_clustering_result$hc_3==3)) )
  result<-c(x,sum(x))
write.csv(result,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_summary.csv")
#################################################################################
看PC NC在哪一群
#################################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_result_39.csv")
write.csv(subset(kirc_clustering_result,kirc_clustering_result$hc_3=="1"),"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_1.csv")
write.csv(subset(kirc_clustering_result,kirc_clustering_result$hc_3=="2"),"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_2.csv")
write.csv(subset(kirc_clustering_result,kirc_clustering_result$hc_3=="3"),"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_3.csv")
kirc_positive_marker<-subset(kirc_clustering_result,
                             kirc_clustering_result[,1]=="CDH1"|
                               kirc_clustering_result[,1]=="FN1"|
                               kirc_clustering_result[,1]=="VIM"|
                               kirc_clustering_result[,1]=="ZEB1"|
                               kirc_clustering_result[,1]=="FOXC1"|
                               kirc_clustering_result[,1]=="SNAI1"|
                               kirc_clustering_result[,1]=="SNAI1"|
                               kirc_clustering_result[,1]=="TWIST1"|
                               kirc_clustering_result[,1]=="TWIST1"|
                               kirc_clustering_result[,1]=="TGFBR1"|
                               kirc_clustering_result[,1]=="FGFR1"|
                               kirc_clustering_result[,1]=="CTNNB1"|
                               kirc_clustering_result[,1]=="MET"|
                               kirc_clustering_result[,1]=="EGFR")
kirc_negative_marker<-subset(kirc_clustering_result,
                             kirc_clustering_result[,1]=="CDH1"|
                               kirc_clustering_result[,1]=="BRMS1"|
                               kirc_clustering_result[,1]=="MED13"|
                               kirc_clustering_result[,1]=="CD81"|
                               kirc_clustering_result[,1]=="KISS1"|
                               kirc_clustering_result[,1]=="NME1"|
                               kirc_clustering_result[,1]=="TP63"|
                               kirc_clustering_result[,1]=="TXNIP"|
                               kirc_clustering_result[,1]=="ARHGDIB"|
                               kirc_clustering_result[,1]=="AKAP11"|
                               kirc_clustering_result[,1]=="TIMP1"|
                               kirc_clustering_result[,1]=="TIMP1"|
                               kirc_clustering_result[,1]=="TIMP3"|
                               kirc_clustering_result[,1]=="TIMP4"|
                               kirc_clustering_result[,1]=="MAP1K4"|
                               kirc_clustering_result[,1]=="DICER1")
write.csv(kirc_positive_marker,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_positive_marker.csv")
write.csv(kirc_negative_marker,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_negative_marker.csv")
#########################################################################
分群結果跑KEGG
#########################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_result_39.csv")
#3
kirc_hc_3<-subset(kirc_clustering_result,kirc_clustering_result$hc_3 == "3")
gene_df <- bitr(kirc_hc_3$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_3<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_3,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/KEGG_sum_hc_3.csv",row.names =F)
#2
kirc_hc_2<-subset(kirc_clustering_result,kirc_clustering_result$hc_3 == "2")
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_2<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/KEGG_sum_hc_2.csv",row.names =F)
#1
kirc_hc_1<-subset(kirc_clustering_result,kirc_clustering_result$hc_3 == "1")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_1<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/KEGG_sum_hc_1.csv",row.names =F)
###########################################################
H3細分
###########################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_result_39.csv")
kirc_hc_3<-subset(kirc_clustering_result,kirc_clustering_result$hc_3 == "3")
data<-kirc_hc_3[,3:41]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/E_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC HC3")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/A_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC HC3") 
dev.off()
#
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_hc_3.png",width = 1000, height = 500, units = "px", pointsize = 24,
    bg = "white")
plot(hc, hang = -2, cex = 0.6)
dev.off()
hc_2<-cutree(hc,k=2)
#
x<-c(length(which(hc_2=="1")),
     length(which(hc_2=="2")) )
sum_hc3_num<-c(x,sum(x))
write.csv(sum_hc3_num,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/sum_hc3_num.csv")
#
kirc_hc_3_clustering_result<-cbind(kirc_hc_3[,1:41],hc_2)
write.csv(kirc_hc_3_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_hc_3_clustering_result.csv",row.names =F)
#########################################################################
kirc_hc_3_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_hc_3_clustering_result.csv")
#2
kirc_hc_2<-subset(kirc_hc_3_clustering_result,kirc_hc_3_clustering_result$hc_2 == "2")
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_2<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/KEGG_sum_hc_32.csv",row.names =F)
#1
kirc_hc_1<-subset(kirc_hc_3_clustering_result,kirc_hc_3_clustering_result$hc_2 == "1")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_1<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/KEGG_sum_hc_31.csv",row.names =F)
###########################################################
heatmap
###########################################################
#35feature
kirc_clustering_39_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_39_clustering_result.csv")
kirc_hc_2<-subset(kirc_clustering_39_clustering_result,kirc_clustering_39_clustering_result$hc_2 == "2")
kirc_clustering_392_hp<-kirc_hc_2[,3:37]
rownames(kirc_clustering_392_hp)<-kirc_hc_2$kirc_merge_all.external_gene_name
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/heatmap35.png",units="px",height=3000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_clustering_392_hp,scale="row",cluster_col=F )
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/heatmap_detail35.png",units="px",height=14000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_clustering_392_hp,scale="row",cluster_col=F )
dev.off()
###########################################################
#39feature
kirc_clustering_39_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_clustering_39_clustering_result.csv")
kirc_hc_2<-subset(kirc_clustering_39_clustering_result,kirc_clustering_39_clustering_result$hc_2 == "2")
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
colnames(kirc_hc_2)[2]<-"external_gene_name"
kirc_clustering_392<-merge(kirc_hc_2[,2:37],kirc_merge_all[,2:6],by="external_gene_name")
kirc_clustering_392_hp<-kirc_clustering_392[,2:40]
rownames(kirc_clustering_392_hp)<-kirc_clustering_392$external_gene_name
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/heatmap.png",units="px",height=3000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_clustering_392_hp,scale="row",cluster_col=F )
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/heatmap_detail.png",units="px",height=14000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_clustering_392_hp,scale="row",cluster_col=F )
dev.off()
#########################################################
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
#
kirc_go_CC_hc2 <- enrichGO(gene = gene_df$ENTREZID, 
                           keyType = "ENTREZID", 
                           OrgDb = 'org.Hs.eg.db', 
                           ont = 'CC', 
                           pAdjustMethod = 'BH', 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.2, 
                           readable=T)
kirc_go_CC_hc2<-summary(kirc_go_CC_hc2)
ratio<-as.numeric(list.rbind(strsplit(as.character(kirc_go_CC_hc2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(kirc_go_CC_hc2[,4]),"/"))[,1])
kirc_go_CC_hc2<-cbind(kirc_go_CC_hc2,ratio)
write.csv(kirc_go_CC_hc2,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_go_CC_hc2.csv")
#
kirc_go_MF_hc2 <- enrichGO(gene = gene_df$ENTREZID, 
                           keyType = "ENTREZID", 
                           OrgDb = 'org.Hs.eg.db', 
                           ont = 'MF', 
                           pAdjustMethod = 'BH', 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.2, 
                           readable=T)
kirc_go_MF_hc2<-summary(kirc_go_MF_hc2)
ratio<-as.numeric(list.rbind(strsplit(as.character(kirc_go_MF_hc2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(kirc_go_MF_hc2[,4]),"/"))[,1])
kirc_go_MF_hc2<-cbind(kirc_go_MF_hc2,ratio)
write.csv(kirc_go_MF_hc2,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_go_MF_hc2.csv")
#
kirc_go_BP_hc2 <- enrichGO(gene = gene_df$ENTREZID, 
                           keyType = "ENTREZID", 
                           OrgDb = 'org.Hs.eg.db', 
                           ont = 'BP', 
                           pAdjustMethod = 'BH', 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.2, 
                           readable=T)
kirc_go_BP_hc2<-summary(kirc_go_BP_hc2)
ratio<-as.numeric(list.rbind(strsplit(as.character(kirc_go_BP_hc2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(kirc_go_BP_hc2[,4]),"/"))[,1])
kirc_go_BP_hc2<-cbind(kirc_go_BP_hc2,ratio)
write.csv(kirc_go_BP_hc2,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_39/kirc_go_BP_hc2.csv")