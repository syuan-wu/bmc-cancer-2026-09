#plot
library(ggplot2)
library(grid)
library(ggpubr)
library(scatterplot3d)
library(pheatmap)
library(easyGgplot2)

#data preproccessing
library(factoextra)
library(dplyr)
library(clusterProfiler)
library(DOSE)
library(rlist)
library(Rmisc)
#pathway
library(msigdbr)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
#######################################################################
feature data前處理
#######################################################################
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
data<-as.matrix(kirc_merge_all[,7:13])
#kirc_marker_coef[,2:31]
#remove TN_FC
data<-data[,-4]
sum(is.na(data))
sum(is.nan(data))
sum(is.infinite(data))
data[is.na(data)] <- 0
data[is.infinite(data)] <- 10000
sum(is.na(data))
sum(is.nan(data))
sum(is.infinite(data))
data<-as.data.frame(data)
data_scale<-scale(data)
#發現不一樣的ensg對到同個external_gene_name,所以用ensg當分群的rowname
#kirc_m1_rawdata<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_m1_rawdata.csv")
#data<-cbind(kirc_m1_rawdata$external_gene_name,as.data.frame(data))
#dup_data<-data[data$`kirc_m1_rawdata$external_gene_name`=="MATR3"|data$`kirc_m1_rawdata$external_gene_name`=="PINX1"|data$`kirc_m1_rawdata$external_gene_name`=="PRAMEF7"|data$`kirc_m1_rawdata$external_gene_name`=="TMSB15B",]
#rownames(data)<-kirc_m1_rawdata$ensembl_gene_id
#data<-data[,-1]
#準備hcluster分群
#######################################################################
1st Clustering
#######################################################################
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/E for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC")
dev.off()
# Avg. Silhouett define number of clustering
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/A for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC") 
dev.off()
#Clustering
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc, hang = -2, cex = 0.6)
hc_n<-cutree(hc,k=6)
#Clustering plot
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_clustering.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC cluster Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
#merge result
kirc_clustering_result<-cbind(kirc_merge_all$external_gene_name,data,hc_n)
write.csv(kirc_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_clustering_result.csv")
#######################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_clustering_result.csv")
for ( i in as.character(c(1:6)) ){
kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc_n == i)
gene_df <- bitr(kirc_hc$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc<-as.data.frame(kegg)
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,4]),"/"))[,1])
KEGG_sum_hc<-cbind(KEGG_sum_hc,ratio)
z<-paste0("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_",i)
dic<-paste0(z,".csv")
write.csv(KEGG_sum_hc,dic)
}
#Score
hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_1.csv")
hc_2<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_2.csv")
hc_3<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_3.csv")
hc_4<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_4.csv")
hc_5<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_5.csv")
hc_6<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_KEGG_sum_hc_6.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_clustering_result.csv")
hc_1<-merge(df,hc_1,by="Description")
hc_2<-merge(df,hc_2,by="Description")
hc_3<-merge(df,hc_3,by="Description")
hc_4<-merge(df,hc_4,by="Description")
hc_5<-merge(df,hc_5,by="Description")
hc_6<-merge(df,hc_6,by="Description")
score1<-c(sum(hc_1$Score),
          sum(hc_2$Score),
          sum(hc_3$Score),
          sum(hc_4$Score),
          sum(hc_5$Score),
          sum(hc_6$Score))
#sum_score*logp*ratio
score2<-c(sum(hc_1$Score*hc_1$ratio*-log(hc_1$pvalue)),
          sum(hc_2$Score*hc_2$ratio*-log(hc_2$pvalue)),
          sum(hc_3$Score*hc_3$ratio*-log(hc_3$pvalue)),
          sum(hc_4$Score*hc_4$ratio*-log(hc_4$pvalue)),
          sum(hc_5$Score*hc_5$ratio*-log(hc_5$pvalue)),
          sum(hc_6$Score*hc_6$ratio*-log(hc_6$pvalue))
          )
number<-c(length(which(kirc_clustering_result$hc_n==1)),
          length(which(kirc_clustering_result$hc_n==2)),
          length(which(kirc_clustering_result$hc_n==3)),
          length(which(kirc_clustering_result$hc_n==4)),
          length(which(kirc_clustering_result$hc_n==5)),
          length(which(kirc_clustering_result$hc_n==6)))
clustering_sum<-data.frame(cluster=c(c(1:6),"sum"),number=c(number,sum(number)),score1=c(score1,sum(score1)),score2=c(score2,sum(score2)) )
write.csv(clustering_sum,"C:/Users/syuan/Desktop/Master/kirc/Clustering/clustering_sum.csv")
#######################################################################
2nd Clustering: hc_4
#######################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_clustering_result.csv")
kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == "4")
data<-kirc_hc[,3:8]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/E for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC4 of KIRC")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/A for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC4 of KIRC") 
dev.off()
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/Clustering_result.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC HC4 cluster Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
hc_n<-cutree(hc,k=5)
kirc_clustering_result<-cbind(gene=as.character(kirc_hc$kirc_merge_all.external_gene_name),data,hc=hc_n)
write.csv(kirc_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_clustering_result.csv")
#definition
for ( i in as.character(c(1:5)) ){
  kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == i)
  gene_df <- bitr(kirc_hc$'gene', fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                   qvalueCutoff = 0.2)
  KEGG_sum_hc<-as.data.frame(kegg)
  ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,4]),"/"))[,1])
  KEGG_sum_hc<-cbind(KEGG_sum_hc,ratio)
  z<-paste0("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_KEGG_sum_hc_",i)
  dic<-paste0(z,".csv")
  write.csv(KEGG_sum_hc,dic)
}
#Score
hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_KEGG_sum_hc_1.csv")
hc_2<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_KEGG_sum_hc_2.csv")
hc_3<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_KEGG_sum_hc_3.csv")
hc_4<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_KEGG_sum_hc_4.csv")
hc_5<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_KEGG_sum_hc_5.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_clustering_result.csv")
hc_1<-merge(df,hc_1,by="Description")
hc_2<-merge(df,hc_2,by="Description")
hc_3<-merge(df,hc_3,by="Description")
hc_4<-merge(df,hc_4,by="Description")
hc_5<-merge(df,hc_5,by="Description")
score1<-c(sum(hc_1$Score),
          sum(hc_2$Score),
          sum(hc_3$Score),
          sum(hc_4$Score),
          sum(hc_5$Score))
#sum_score*logp*ratio
score2<-c(sum(hc_1$Score*hc_1$ratio*-log(hc_1$pvalue)),
          sum(hc_2$Score*hc_2$ratio*-log(hc_2$pvalue)),
          sum(hc_3$Score*hc_3$ratio*-log(hc_3$pvalue)),
          sum(hc_4$Score*hc_4$ratio*-log(hc_4$pvalue)),
          sum(hc_5$Score*hc_5$ratio*-log(hc_5$pvalue))
)
number<-c(length(which(kirc_clustering_result$hc==1)),
          length(which(kirc_clustering_result$hc==2)),
          length(which(kirc_clustering_result$hc==3)),
          length(which(kirc_clustering_result$hc==4)),
          length(which(kirc_clustering_result$hc==5))
          )
clustering_sum<-data.frame(cluster=c(c(1:5),"sum"),number=c(number,sum(number)),score1=c(score1,sum(score1)),score2=c(score2,sum(score2)) )
write.csv(clustering_sum,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/clustering_sum.csv")
#######################################################################
3rd Clustering: hc_41
#######################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_clustering_result.csv")
kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == "1")
data<-kirc_hc[,4:9]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/E for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC41 of KIRC")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/A for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC41 of KIRC") 
dev.off()
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/Clustering_result.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC HC 41 cluster Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
hc_n<-cutree(hc,k=4)
kirc_clustering_result<-cbind(gene=as.character(kirc_hc$gene),data,hc=hc_n)
write.csv(kirc_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_clustering_result.csv")
#########################
for ( i in as.character(c(1:4)) ){
  kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == i)
  gene_df <- bitr(kirc_hc$'gene', fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                   qvalueCutoff = 0.2)
  KEGG_sum_hc<-as.data.frame(kegg)
  ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,4]),"/"))[,1])
  KEGG_sum_hc<-cbind(KEGG_sum_hc,ratio)
  z<-paste0("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_",i)
  dic<-paste0(z,".csv")
  write.csv(KEGG_sum_hc,dic)
}
#Score
hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_1.csv")
hc_2<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_2.csv")
hc_3<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_3.csv")
hc_4<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_4.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_clustering_result.csv")
hc_1<-merge(df,hc_1,by="Description")
hc_2<-merge(df,hc_2,by="Description")
hc_3<-merge(df,hc_3,by="Description")
hc_4<-merge(df,hc_4,by="Description")
score1<-c(sum(hc_1$Score),
          sum(hc_2$Score),
          sum(hc_3$Score),
          sum(hc_4$Score))
#sum_score*logp*ratio
score2<-c(sum(hc_1$Score*hc_1$ratio*-log(hc_1$pvalue)),
          sum(hc_2$Score*hc_2$ratio*-log(hc_2$pvalue)),
          sum(hc_3$Score*hc_3$ratio*-log(hc_3$pvalue)),
          sum(hc_4$Score*hc_4$ratio*-log(hc_4$pvalue))
)
number<-c(length(which(kirc_clustering_result$hc==1)),
          length(which(kirc_clustering_result$hc==2)),
          length(which(kirc_clustering_result$hc==3)),
          length(which(kirc_clustering_result$hc==4))
)
clustering_sum<-data.frame(cluster=c(c(1:4),"sum"),number=c(number,sum(number)),score1=c(score1,sum(score1)),score2=c(score2,sum(score2)) )
write.csv(clustering_sum,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/clustering_sum.csv")
#######################################################################

#################################################################################
kirc_41_3277用30ft分
#################################################################################
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_4/kirc_clustering_result.csv")
kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == "1")
colnames(kirc_marker_coef)[1]<-"gene"
kirc_hc<-merge(kirc_marker_coef,kirc_hc,by="gene")
data<-kirc_hc[,2:31]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/E for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC41 of KIRC")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/A for HC.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC41 of KIRC") 
dev.off()
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/Clustering_result.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC HC 41 cluster Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
hc_n<-cutree(hc,k=3)
kirc_clustering_result<-cbind(gene=as.character(kirc_hc$gene),data,hc=hc_n)
write.csv(kirc_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_clustering_result.csv")
#########################
for ( i in as.character(c(1:3)) ){
  kirc_hc<-subset(kirc_clustering_result,kirc_clustering_result$hc == i)
  gene_df <- bitr(kirc_hc$'gene', fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                   pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                   qvalueCutoff = 0.2)
  KEGG_sum_hc<-as.data.frame(kegg)
  ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc[,4]),"/"))[,1])
  KEGG_sum_hc<-cbind(KEGG_sum_hc,ratio)
  z<-paste0("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_",i)
  dic<-paste0(z,".csv")
  write.csv(KEGG_sum_hc,dic)
}
#Score
hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_1.csv")
hc_2<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_2.csv")
hc_3<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_KEGG_sum_hc_3.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/kirc_clustering_result.csv")
hc_1<-merge(df,hc_1,by="Description")
hc_2<-merge(df,hc_2,by="Description")
hc_3<-merge(df,hc_3,by="Description")
score1<-c(sum(hc_1$Score),
          sum(hc_2$Score),
          sum(hc_3$Score))
#sum_score*logp*ratio
score2<-c(sum(hc_1$Score*hc_1$ratio*-log(hc_1$pvalue)),
          sum(hc_2$Score*hc_2$ratio*-log(hc_2$pvalue)),
          sum(hc_3$Score*hc_3$ratio*-log(hc_3$pvalue))
)
number<-c(length(which(kirc_clustering_result$hc==1)),
          length(which(kirc_clustering_result$hc==2)),
          length(which(kirc_clustering_result$hc==3))
)
clustering_sum<-data.frame(cluster=c(c(1:3),"sum"),number=c(number,sum(number)),score1=c(score1,sum(score1)),score2=c(score2,sum(score2)) )
write.csv(clustering_sum,"C:/Users/syuan/Desktop/Master/kirc/Clustering/hc_41/clustering_sum.csv")
########################################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_result.csv")
write.csv(subset(kirc_clustering_result,kirc_clustering_result$hc_2=="1"),"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_1.csv")
write.csv(subset(kirc_clustering_result,kirc_clustering_result$hc_2=="2"),"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_2.csv")
kirc_positive_marker<-subset(kirc_clustering_result,
                              kirc_clustering_result[,2]=="CDH2"|
                              kirc_clustering_result[,2]=="FN1"|
                              kirc_clustering_result[,2]=="VIM"|
                              kirc_clustering_result[,2]=="ZEB2"|
                              kirc_clustering_result[,2]=="FOXC2"|
                              kirc_clustering_result[,2]=="SNAI1"|
                              kirc_clustering_result[,2]=="SNAI2"|
                              kirc_clustering_result[,2]=="TWIST1"|
                              kirc_clustering_result[,2]=="TWIST2"|
                              kirc_clustering_result[,2]=="TGFBR1"|
                              kirc_clustering_result[,2]=="FGFR1"|
                              kirc_clustering_result[,2]=="CTNNB1"|
                              kirc_clustering_result[,2]=="MET"|
                              kirc_clustering_result[,2]=="EGFR")
kirc_negative_marker<-subset(kirc_clustering_result,
                             kirc_clustering_result[,2]=="CDH1"|
                              kirc_clustering_result[,2]=="BRMS1"|
                              kirc_clustering_result[,2]=="MED23"|
                              kirc_clustering_result[,2]=="CD82"|
                              kirc_clustering_result[,2]=="KISS1"|
                              kirc_clustering_result[,2]=="NME1"|
                              kirc_clustering_result[,2]=="TP63"|
                              kirc_clustering_result[,2]=="TXNIP"|
                              kirc_clustering_result[,2]=="ARHGDIB"|
                              kirc_clustering_result[,2]=="AKAP12"|
                               kirc_clustering_result[,2]=="TIMP1"|
                               kirc_clustering_result[,2]=="TIMP2"|
                               kirc_clustering_result[,2]=="TIMP3"|
                               kirc_clustering_result[,2]=="TIMP4"|
                               kirc_clustering_result[,2]=="MAP2K4"|
                               kirc_clustering_result[,2]=="DICER1")
write.csv(kirc_positive_marker,"C:/Users/syuan/Desktop/Master/kirc/kirc_positive_marker.csv")
write.csv(kirc_negative_marker,"C:/Users/syuan/Desktop/Master/kirc/kirc_negative_marker.csv")
#17883
#kirc_clustering_merge<-merge(subset(kirc_clustering_result,hc_3=="1"|hc_3=="3"),subset(kirc_clustering_result,km_3=="2"|km_3=="3"),by="kirc_merge_all.external_gene_name")
#write.csv(kirc_clustering_merge,"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_merge.csv")
#################################################################################
分群結果視覺化20200608
#################################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_result.csv")
my_comparisons <- list(c("1", "2"), c("2", "3"), c("1", "3"))
colnames(kirc_clustering_result)[3:8]<-c("Stage_correlation_coefficient","Stage_fval","m1m0_FC","Survival_curve_pval","THN_FC","Survival_Hazard_Ratio")
result<-list()
for (i in c(3:38)){
print(i)
  p<-ggboxplot(kirc_clustering_result,x="hc_3",y=colnames(kirc_clustering_result)[i],palette="jama",color="hc_3")+stat_compare_means(method = "wilcox.test",comparisons = my_comparisons)+annotate("text", x = -Inf, y = Inf, label = "Wilcox.test", hjust = -0.2,                                                                                                                                                                                                    vjust = 2)
  #add = "jitter"
  #add p-value
result<-c(result,list(p))
}
result[17]
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/clustering_feature2.png",width = 2000, height = 1800, units = "px", pointsize = 24,
    bg = "white")
multiplot(plotlist = result, cols = 6)
dev.off()
?grid.arrange
################################################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_result.csv")
p6<-ggplot()+geom_point(data=data,aes(x=SNAI1,y=SNAI2))+
  ggtitle("H Clustering")+
  facet_wrap(~ hc_3, ncol = 3)
p7<-ggplot()+geom_point(data=data,
                        aes(x=SNAI1,y=SNAI2,color=as.factor(hc_3)))+ggtitle("H Clustering")
p8<-ggplot()+geom_point(data=data,
                        aes(x=TWIST1,y=TWIST2,color=as.factor(hc_3)))+ggtitle("H Clustering")
p9<-ggplot()+geom_point(data=data,
                        aes(x=CDH1,y=CDH2,color=as.factor(hc_3)))+ggtitle("H Clustering")
p10<-ggplot()+geom_point(data=data,aes(x=m1m0_FC,y=kirc_TN_FC))+
  ggtitle("H Clustering")+
  facet_wrap(~ hc_3, ncol = 3)
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_plot10.png",width = 900, height = 600, units = "px", pointsize = 24,
    bg = "white")
plot(p10)
dev.off()
#3d plot P
png("C:/Users/syuan/Desktop/Master/kirc/kirc_3Dplot1.png",width = 1800, height = 1200, units = "px", pointsize = 24,
    bg = "white")
head(kirc_clustering_result)
colors <- c("#FF44AA", "#00BBFF", "#00DD00")
colors <- colors[as.numeric(as.character(kirc_clustering_result$hc_3))]
s3d<-scatterplot3d(kirc_clustering_result$TWIST1,
                   kirc_clustering_result$TWIST2,
                   kirc_clustering_result$SNAI1,
                   xlab="TWIST1", ylab="TWIST2", zlab="SNAI1",pch = 16, color=colors)
legend("right", legend = as.character(unique(kirc_clustering_result$hc_3)),
       col =  c("#FF44AA", "#00BBFF", "#00DD00"),pch = 16,title="clustering group")
dev.off()
#3d plot N
png("C:/Users/syuan/Desktop/Master/kirc/kirc_3Dplot3.png",width = 1800, height = 1200, units = "px", pointsize = 24,
    bg = "white")
head(kirc_clustering_result)
colors <- c("#FF44AA", "#00BBFF", "#00DD00")
colors <- colors[as.numeric(as.character(kirc_clustering_result$hc_3))]
s3d<-scatterplot3d(kirc_clustering_result$MAP2K4,
                   kirc_clustering_result$TIMP4,
                   kirc_clustering_result$TIMP3,
                   xlab="MAP2K4", ylab="TIMP4", zlab="TIMP3",pch = 16, color=colors)
legend("right", legend = as.character(unique(kirc_clustering_result$hc_3)),
       col =  c("#FF44AA", "#00BBFF", "#00DD00"),pch = 16,title="clustering group")
dev.off()
#boxplot
boxplotlist<-list()
for (i in c(3:37)){
print(i)
a<-list( ggplot(kirc_clustering_result,aes(x=as.factor(hc_3),y=kirc_clustering_result[,as.numerici]))+
  geom_boxplot(width=.5)+labs(x="Cluster",y=colnames(kirc_clustering_result)[i])+
  theme(plot.title=element_text(hjust = 0.5,face="bold",size=5)) )
class(a)
boxplotlist<-c(boxplotlist,a)
}
boxplotlist[[1]]
boxplotlist[[2]]
class(boxplotlist)
png("C:/Users/syuan/Desktop/Master/kirc/Clustering/kirc_feature_37.png",width = 600, height = 400, units = "px", pointsize = 24,
    bg = "white")
ggplot(kirc_clustering_result,aes(x=as.factor(hc_3),y=kirc_clustering_result[,37]))+
  geom_boxplot(width=.5)+labs(x="Cluster",y=colnames(kirc_clustering_result)[37])+
  theme(plot.title=element_text(hjust = 0.5,face="bold",size=5))+theme_bw()
dev.off()
#########################################################################
#2
kirc_hc_2<-subset(kirc_clustering_result,kirc_clustering_result$hc_2 == "2")
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_2<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_2.csv",row.names =F)
#1
kirc_hc_1<-subset(kirc_clustering_result,kirc_clustering_result$hc_2 == "1")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_1<-as.data.frame(kegg)
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_1.csv",row.names =F)
###########################################################
png("C:/Users/syuan/Desktop/Master/kirc/kirc_kegg_dotplot.png",units="px",height=3000,width=4000,res=250)
dotplot(kegg, showCategory=5)
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/kirc_kegg_network.png",units="px",height=1000,width=1600,res=200)
cnetplot(kegg,showCategory = 3,categorySize="pvalue")
dev.off()
###########################################################
kirc_hc3_KEGG_net_3<-subset(gene_df,gene_df$ENTREZID=="5154"|
         gene_df$ENTREZID=="5228"|
         gene_df$ENTREZID=="5156"|
         gene_df$ENTREZID=="2316"|
         gene_df$ENTREZID=="5159"|
         gene_df$ENTREZID=="7424"|
         gene_df$ENTREZID=="5595"|
         gene_df$ENTREZID=="207"|
         gene_df$ENTREZID=="3725"|
        gene_df$ENTREZID=="5895")
write.csv(kirc_hc3_KEGG_net_3,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc3_KEGG_net.csv")
#Focal adhesion
pathview(gene.data  = gene_df$ENTREZID,
         pathway.id = "hsa04510",
         species    = "hsa",
         limit      = list(gene=c(0,5), cpd=1) )
#MAPK signaling pathway
pathview(gene.data  = gene_df$ENTREZID,
         pathway.id = "hsa04010",
         species    = "hsa",
         limit      = list(gene=c(0,5), cpd=1) )
#why the color of gene on KEGG pathway plot are diff?
#the order of the genes?
which(gene_df$SYMBOL=="SHC1")
which(gene_df$SYMBOL=="SHC3")
which(gene_df$SYMBOL=="SHC4")
which(gene_df$SYMBOL=="ITGA11")
###########################################################
#第三群的每個feature
kirc_clustering_3<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_3.csv")
ggplot(kirc_clustering_3,aes(x=as.factor(hc_3),y=kirc_clustering_3[,as.numerici]))+
  geom_boxplot(width=.5)+labs(x="Cluster",y=colnames(kirc_clustering_3)[4])+
  theme(plot.title=element_text(hjust = 0.5,face="bold",size=5))
###########################################################
KEGG結果用ratio排序
###########################################################
KEGG_sum_hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_1.csv")
KEGG_sum_hc_2<-read.csv("C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_2.csv")
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,4]),"/"))[,1])
KEGG_sum_hc_1<-cbind(KEGG_sum_hc_1,ratio)
KEGG_sum_hc_1<-KEGG_sum_hc_1%>%arrange(desc(ratio))
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_1_order.csv")
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,4]),"/"))[,1])
KEGG_sum_hc_2<-cbind(KEGG_sum_hc_2,ratio)
KEGG_sum_hc_2<-KEGG_sum_hc_2%>%arrange(desc(ratio))
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_2_order.csv")
###########################################################
H2細分
###########################################################
kirc_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_result.csv")
kirc_hc_2<-subset(kirc_clustering_result,kirc_clustering_result$hc_2 == "2")
data<-kirc_hc_2[,3:38]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/E_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC HC2")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/A_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC HC2") 
dev.off()
#
dd <- dist(data_scale, method = "euclidean")
view(dd)
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/kirc_hc_2_clustering.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC cluster 2 Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
hc_2<-cutree(hc,k=2)
#
x<-c(length(which(hc_2=="1")),
     length(which(hc_2=="2"))
     )
clustering_summary<-c(x,sum(x))
write.csv(clustering_summary,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/clustering_summary.csv")
#
kirc_hc_2_clustering_result<-cbind(kirc_hc_2[,1:38],hc_2)
write.csv(kirc_hc_2_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/kirc_hc_2_clustering_result.csv",row.names =F)
#########################################################################
kirc_hc_2_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/kirc_hc_2_clustering_result.csv")
#2-1
kirc_hc_1<-subset(kirc_hc_2_clustering_result,kirc_hc_2_clustering_result$hc_2 == "1")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_1<-as.data.frame(kegg)
#2-2
kirc_hc_2<-subset(kirc_hc_2_clustering_result,kirc_hc_2_clustering_result$hc_2 == "2")
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_2<-as.data.frame(kegg)
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,4]),"/"))[,1])
KEGG_sum_hc_1<-cbind(KEGG_sum_hc_1,ratio)
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/KEGG_sum_hc_1_order.csv")
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,4]),"/"))[,1])
KEGG_sum_hc_2<-cbind(KEGG_sum_hc_2,ratio)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/KEGG_sum_hc_2_order.csv")
#
write.csv(kirc_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/kirc_hc_1.csv",row.names =F)
write.csv(kirc_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/kirc_hc_2.csv",row.names =F)
###########################################################
H2-2細分
###########################################################
kirc_hc_2<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/kirc_hc_2.csv")
data<-kirc_hc_2[,3:38]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/E_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC HC2-2")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/A_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC HC2-2") 
dev.off()
#
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_22_clustering.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC cluster 2-2 Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
hc_2<-cutree(hc,k=2)
#
x<-c(length(which(hc_2=="1")),
     length(which(hc_2=="2"))
)
clustering_summary<-c(x,sum(x))
write.csv(clustering_summary,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/clustering_summary.csv")
#
kirc_hc_22_clustering_result<-cbind(kirc_hc_2[,1:38],hc_2)
write.csv(kirc_hc_22_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_22_clustering_result.csv",row.names =F)
#########################################################################
kirc_hc_22_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_22_clustering_result.csv")
#2-2-1
kirc_hc_1<-subset(kirc_hc_22_clustering_result,kirc_hc_22_clustering_result$hc_2 == "1")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_1<-as.data.frame(kegg)
#2-2-2
kirc_hc_2<-subset(kirc_hc_22_clustering_result,kirc_hc_22_clustering_result$hc_2 == "2")
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_2<-as.data.frame(kegg)
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,4]),"/"))[,1])
KEGG_sum_hc_1<-cbind(KEGG_sum_hc_1,ratio)
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/KEGG_sum_hc_1_order.csv")
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,4]),"/"))[,1])
KEGG_sum_hc_2<-cbind(KEGG_sum_hc_2,ratio)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/KEGG_sum_hc_2_order.csv")
#
write.csv(kirc_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_1.csv",row.names =F)
write.csv(kirc_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_2.csv",row.names =F)
###########################################################
H2-2-1細分
###########################################################
kirc_hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_1.csv")
data<-kirc_hc_1[,3:38]
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/E_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC HC2-2-1")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/A_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC HC2-2-1") 
dev.off()
#
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/kirc_hc_221_clustering.png",width = 1500, height = 800, units = "px", pointsize = 24,
    bg = "white")
plot(main="KIRC cluster 2-2-1 Dendrogram",xlab="Ensemble Gene Id",hc, hang = -2, cex = 0.6)#+abline(h=300, col="red")
dev.off()
hc_3<-cutree(hc,k=3)
#
x<-c(length(which(hc_3=="1")),
     length(which(hc_3=="2")),
     length(which(hc_3=="3"))
)
clustering_summary<-c(x,sum(x))
write.csv(clustering_summary,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/clustering_summary.csv")
#
kirc_hc_221_clustering_result<-cbind(kirc_hc_1[,1:38],hc_3)
write.csv(kirc_hc_221_clustering_result,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/kirc_hc_221_clustering_result.csv",row.names =F)
#########################################################################
kirc_hc_221_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/kirc_hc_221_clustering_result.csv")
#2-2-1-1
kirc_hc_1<-subset(kirc_hc_221_clustering_result,kirc_hc_221_clustering_result$hc_3 == "1")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_1<-as.data.frame(kegg)
#2-2-1-2
kirc_hc_2<-subset(kirc_hc_221_clustering_result,kirc_hc_221_clustering_result$hc_3 == "2")
gene_df <- bitr(kirc_hc_2$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_2<-as.data.frame(kegg)
#2-2-1-3
kirc_hc_3<-subset(kirc_hc_221_clustering_result,kirc_hc_221_clustering_result$hc_3 == "3")
gene_df <- bitr(kirc_hc_1$'kirc_merge_all.external_gene_name', fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
kegg<-enrichKEGG(gene_df$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
KEGG_sum_hc_3<-as.data.frame(kegg)
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_1[,4]),"/"))[,1])
KEGG_sum_hc_1<-cbind(KEGG_sum_hc_1,ratio)
write.csv(KEGG_sum_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/KEGG_sum_hc_1_order.csv")
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_2[,4]),"/"))[,1])
KEGG_sum_hc_2<-cbind(KEGG_sum_hc_2,ratio)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/KEGG_sum_hc_2_order.csv")
#
ratio<-as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_3[,3]),"/"))[,1])/as.numeric(list.rbind(strsplit(as.character(KEGG_sum_hc_3[,4]),"/"))[,1])
KEGG_sum_hc_3<-cbind(KEGG_sum_hc_3,ratio)
write.csv(KEGG_sum_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/KEGG_sum_hc_3_order.csv")
#
write.csv(kirc_hc_1,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/kirc_hc_1.csv",row.names =F)
write.csv(kirc_hc_2,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/kirc_hc_2.csv",row.names =F)
write.csv(kirc_hc_3,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/kirc_hc_3.csv",row.names =F)
###########################################################
20200608 heatmap
###########################################################
#35feature
kirc_hc_3_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/kirc_hc_3_clustering_result.csv")
kirc_hc_2<-subset(kirc_hc_3_clustering_result,kirc_hc_3_clustering_result$hc_2 == "2")
kirc_hc_32_hp<-kirc_hc_2[,3:37]
rownames(kirc_hc_32_hp)<-kirc_hc_2$kirc_merge_all.external_gene_name
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/heatmap35.png",units="px",height=3000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_hc_32_hp,scale="row",cluster_col=F )
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/heatmap_detail35.png",units="px",height=14000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_hc_32_hp,scale="row",cluster_col=F )
dev.off()
###########################################################
#39feature
kirc_hc_3_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/kirc_hc_3_clustering_result.csv")
kirc_hc_2<-subset(kirc_hc_3_clustering_result,kirc_hc_3_clustering_result$hc_2 == "2")
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
colnames(kirc_hc_2)[2]<-"external_gene_name"
kirc_hc_32<-merge(kirc_hc_2[,2:37],kirc_merge_all[,2:6],by="external_gene_name")
kirc_hc_32_hp<-kirc_hc_32[,2:40]
rownames(kirc_hc_32_hp)<-kirc_hc_32$external_gene_name
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/heatmap.png",units="px",height=3000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_hc_32_hp,scale="row",cluster_col=F )
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_3/heatmap_detail.png",units="px",height=14000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",kirc_hc_32_hp,scale="row",cluster_col=F )
dev.off()
#########################################################
kirc_hc_2
data<-as.matrix(kirc_hc_2[,3:7])
sum(is.na(data))
sum(is.nan(data))
sum(is.infinite(data))
data[is.na(data)] <- 0
data[is.infinite(data)] <- 100
data<-as.data.frame(data)
data_scale<-scale(data)
#Elbow Method
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_32/E_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 15          # max number of clusters to consider
) + labs(title="Elbow Method for HC of KIRC3-2")
dev.off()
# Avg. Silhouett
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_32/A_method.png",width = 500, height = 300, units = "px", pointsize = 24,
    bg = "white")
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC3-2") 
dev.off()
#######################################################################
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_32/kirc_hc_32.png",width = 1000, height = 500, units = "px", pointsize = 24,
    bg = "white")
plot(hc, hang = -2, cex = 0.6)
dev.off()
hc_2<-cutree(hc,k=2)
##############################################################
量化分群-kirc
##############################################################
KEGG_sum_hc_1_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_1_order.csv")
KEGG_sum_hc_2_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/KEGG_sum_hc_2_order.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
KEGG_sum_hc_1_order<-merge(KEGG_sum_hc_1_order,df,by="Description")
KEGG_sum_hc_2_order<-merge(KEGG_sum_hc_2_order,df,by="Description")
#sum_score
x<-c(sum(KEGG_sum_hc_1_order$Score),
  sum(KEGG_sum_hc_2_order$Score))
score1<-c(x,sum(x))
#sum_score*logp*ratio
y<-c(sum(KEGG_sum_hc_1_order$Score*KEGG_sum_hc_1_order$ratio*-log(KEGG_sum_hc_1_order$pvalue)),
     sum(KEGG_sum_hc_2_order$Score*KEGG_sum_hc_2_order$ratio*-log(KEGG_sum_hc_2_order$pvalue)))
score2<-c(y,sum(y))
kirc_clustering_summary<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_summary.csv")
write.csv(cbind(kirc_clustering_summary,score1,score2),"C:/Users/syuan/Desktop/Master/kirc/kirc_clustering_summary.csv")
##############################################################
量化分群-kirc-2
##############################################################
KEGG_sum_hc_1_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/KEGG_sum_hc_1_order.csv")
KEGG_sum_hc_2_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/KEGG_sum_hc_2_order.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
KEGG_sum_hc_1_order<-merge(KEGG_sum_hc_1_order,df,by="Description")
KEGG_sum_hc_2_order<-merge(KEGG_sum_hc_2_order,df,by="Description")
#sum_score
x<-c(sum(KEGG_sum_hc_1_order$Score),
     sum(KEGG_sum_hc_2_order$Score))
score1<-c(x,sum(x))
#sum_score*logp*ratio
y<-c(sum(KEGG_sum_hc_1_order$Score*KEGG_sum_hc_1_order$ratio*-log(KEGG_sum_hc_1_order$pvalue)),
     sum(KEGG_sum_hc_2_order$Score*KEGG_sum_hc_2_order$ratio*-log(KEGG_sum_hc_2_order$pvalue)) )
score2<-c(y,sum(y))
clustering_summary<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/clustering_summary.csv")
write.csv(cbind(clustering_summary,score1,score2),"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_2/clustering_summary.csv")
##############################################################
量化分群-kirc-22
##############################################################
KEGG_sum_hc_1_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/KEGG_sum_hc_1_order.csv")
KEGG_sum_hc_2_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/KEGG_sum_hc_2_order.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
KEGG_sum_hc_1_order<-merge(KEGG_sum_hc_1_order,df,by="Description")
KEGG_sum_hc_2_order<-merge(KEGG_sum_hc_2_order,df,by="Description")
#sum_score
x<-c(sum(KEGG_sum_hc_1_order$Score),
     sum(KEGG_sum_hc_2_order$Score))
score1<-c(x,sum(x))
#sum_score*logp*ratio
y<-c(sum(KEGG_sum_hc_1_order$Score*KEGG_sum_hc_1_order$ratio*-log(KEGG_sum_hc_1_order$pvalue)),
     sum(KEGG_sum_hc_2_order$Score*KEGG_sum_hc_2_order$ratio*-log(KEGG_sum_hc_2_order$pvalue)) )
score2<-c(y,sum(y))
clustering_summary<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/clustering_summary.csv")
write.csv(cbind(clustering_summary,score1,score2),"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/clustering_summary.csv")
##############################################################
量化分群-kirc-221
##############################################################
KEGG_sum_hc_1_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/KEGG_sum_hc_1_order.csv")
KEGG_sum_hc_2_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/KEGG_sum_hc_2_order.csv")
KEGG_sum_hc_3_order<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/KEGG_sum_hc_3_order.csv")
df<-read.csv("C:/Users/syuan/Desktop/Master/KEGG/kegg_pathway.csv")
KEGG_sum_hc_1_order<-merge(KEGG_sum_hc_1_order,df,by="Description")
KEGG_sum_hc_2_order<-merge(KEGG_sum_hc_2_order,df,by="Description")
KEGG_sum_hc_3_order<-merge(KEGG_sum_hc_3_order,df,by="Description")
#sum_score
x<-c(sum(KEGG_sum_hc_1_order$Score),
     sum(KEGG_sum_hc_2_order$Score),
     sum(KEGG_sum_hc_3_order$Score)
     )
score1<-c(x,sum(x))
#sum_score*logp*ratio
y<-c(sum(KEGG_sum_hc_1_order$Score*KEGG_sum_hc_1_order$ratio*-log(KEGG_sum_hc_1_order$pvalue)),
     sum(KEGG_sum_hc_2_order$Score*KEGG_sum_hc_2_order$ratio*-log(KEGG_sum_hc_2_order$pvalue)),
     sum(KEGG_sum_hc_3_order$Score*KEGG_sum_hc_3_order$ratio*-log(KEGG_sum_hc_3_order$pvalue)) )
score2<-c(y,sum(y))
clustering_summary<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/clustering_summary.csv")
write.csv(cbind(clustering_summary,score1,score2),"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_221/clustering_summary.csv")
#######################################################################
heatmap-GLM_n36_EMT_top_397
#######################################################################
EMT_top<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/predict_result_n36/EMT_top.csv", stringsAsFactors = FALSE, row.names = NULL)
kirc_hc_1<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_1.csv", stringsAsFactors = FALSE, row.names = NULL)
colnames(EMT_top)[2]<-colnames(kirc_hc_1)[2]
df<-merge(EMT_top,kirc_hc_1,by=colnames(kirc_hc_1)[2])
df_1<-df[,4:39]
rownames(df_1)<-as.character(df$kirc_merge_all.external_gene_name)
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/heatmap_397_column.png",units="px",height=6000,width=2000,res=200)
pheatmap(main="KIRC HC 2-2-1-EMT-397",scale="column",df_1)
dev.off()
#
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_32/heatmap_652_column_d.png",units="px",height=14000,width=2000,res=200)
pheatmap(main="KIRC HC 3-2-652",data_scale,scale="column",cluster_col=F,color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-9,9,3),breaks=bk )
dev.off()
#######################################################################
397 v.s. 19615
#######################################################################
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
colnames(kirc_marker_coef)[1]<-"external_gene_name"
df<-cbind(kirc_marker_coef,kirc_merge_all)
EMT_top<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/predict_result_n36/EMT_top.csv", stringsAsFactors = FALSE, row.names = NULL)
result<-NULL
for(k in c(1:397)){
  print(k)
  n<-which(df$external_gene_name==EMT_top$top[k])
  result<-c(result,n)
}
group<-rep("Others",19615)
group[result]<-"Candidates"
df<-data.frame(df,group=group)
df<-df[,-c(32:37)]
result<-NULL
for (i in c(2:38)){
  name<-colnames(df)[i]
#  dic<-paste0("C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/",name,".png")
#  png(dic,width = 600, height = 500, units = "px", pointsize = 24,bg = "white")
  p<-ggboxplot( df ,x="group",y=name,palette="jama",color="group",title="")+stat_compare_means(method = "wilcox.test",comparisons = list(c("Others", "Candidates")) )+theme_bw(base_size = 18)
result<-c(result,list(p))
}
png("C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/multiplot_397_marker.png",width = 2400, height = 1800, units = "px", pointsize = 24,bg = "white")
multiplot(result[[1]],result[[2]],result[[3]],result[[4]],
          result[[5]],result[[6]],result[[7]],result[[8]],
          result[[9]],result[[10]],result[[11]],result[[12]],
          result[[13]],result[[14]],result[[15]],result[[16]],
          result[[17]],result[[18]],result[[19]],result[[20]],
          result[[21]],result[[22]],result[[23]],result[[24]],
          result[[25]],result[[26]],result[[27]],result[[28]],
          result[[29]],result[[30]],
          cols=6)
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/multiplot_397_parameter.png",width = 2400, height = 900, units = "px", pointsize = 24,bg = "white")
multiplot(result[[31]],result[[32]],result[[33]],result[[34]],result[[35]],result[[36]],result[[37]],cols=4 )
dev.off()
write.csv(df,"C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/group_397.csv")
#######################################################################
1858 v.s. 19615
#######################################################################
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
colnames(kirc_marker_coef)[1]<-"external_gene_name"
df<-cbind(kirc_marker_coef,kirc_merge_all)
kirc_hc_22_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_22/kirc_hc_22_clustering_result.csv", stringsAsFactors = FALSE, row.names = NULL)
kirc_1858<-subset(kirc_hc_22_clustering_result,hc_2=="1")
result<-NULL
for(k in c(1:1858)){
  print(k)
  n<-which(df$external_gene_name==kirc_1858$kirc_merge_all.external_gene_name[k])
  result<-c(result,n)
}
group<-rep("Others",19615)
group[result]<-"Candidates"
df<-data.frame(df,group=group)
df<-df[,-c(32:37)]
result<-NULL
for (i in c(2:38)){
  name<-colnames(df)[i]
  #  dic<-paste0("C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/",name,".png")
  #  png(dic,width = 600, height = 500, units = "px", pointsize = 24,bg = "white")
  p<-ggboxplot( df ,x="group",y=name,palette="jama",color="group",title="")+stat_compare_means(method = "wilcox.test",comparisons = list(c("Others", "Candidates")) )+theme_bw(base_size = 18)
  result<-c(result,list(p))
}
png("C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/multiplot_1858_marker.png",width = 2400, height = 1800, units = "px", pointsize = 24,bg = "white")
multiplot(result[[1]],result[[2]],result[[3]],result[[4]],
          result[[5]],result[[6]],result[[7]],result[[8]],
          result[[9]],result[[10]],result[[11]],result[[12]],
          result[[13]],result[[14]],result[[15]],result[[16]],
          result[[17]],result[[18]],result[[19]],result[[20]],
          result[[21]],result[[22]],result[[23]],result[[24]],
          result[[25]],result[[26]],result[[27]],result[[28]],
          result[[29]],result[[30]],
          cols=6)
dev.off()
png("C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/multiplot_1858_parameter.png",width = 2400, height = 900, units = "px", pointsize = 24,bg = "white")
multiplot(result[[31]],result[[32]],result[[33]],result[[34]],result[[35]],result[[36]],result[[37]],cols=4 )
dev.off()
write.csv(df,"C:/Users/syuan/Desktop/Master/kirc/clustering_boxplot/group_1858.csv")