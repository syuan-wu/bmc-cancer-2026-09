library(ggplot2)
library(grid)
library(scatterplot3d)
#############################################################################
kirc_merge_all<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_merge_all.csv")
kirc_filtered<-subset(kirc_merge_all,F_value<=0.05&pval<=0.05 )
result<-NULL
for (i in c(1:nrow(kirc_filtered))){
  print(i)
  value<-kirc_merge_all$m1m0_FC[i] * kirc_merge_all$kirc_TN_FC[i]
  result<-c(result,value)
}
FC_value<-result
kirc_filtered<-cbind(kirc_filtered,FC_value)
write.csv(kirc_filtered,"C:/Users/syuan/Desktop/Master/kirc/kirc_filtered.csv")
#############################################################################
kirc_filtered<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filtered.csv")
#5424
kirc_marker_coef<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_marker_coef.csv")
#19615
kirc_marker<-colnames(kirc_marker_coef)[2:31]
#30 markers
kirc_marker<-as.data.frame(kirc_marker)
colnames(kirc_marker)<-"external_gene_name"
kirc_filter_marker<-merge(kirc_filtered,kirc_marker,by="external_gene_name")
#知道filter後30個marker剩下11個
kirc_filter_marker<-as.character(kirc_filter_marker$external_gene_name)
#只挑出11個算coefficient
result<-NULL
for(i in c(1:11)){
x<-kirc_marker_coef[,which( colnames(kirc_marker_coef)==kirc_filter_marker[i] )]
result<-cbind(result,x)
}
colnames(result)<-kirc_filter_marker
kirc_filtered_marker_coef<-cbind(as.character(kirc_marker_coef[,1]),result)
colnames(kirc_filtered_marker_coef)[1]<-"external_gene_name"
#filtered+11個marker coefficient合併算clustering
data<-merge(kirc_filtered,as.data.frame(kirc_filtered_marker_coef),by="external_gene_name")
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
n_data<-t(apply(data[,8:24],1, as.numeric))
data_scale<-scale(n_data)
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
) + labs(title="Elbow Method for HC of KIRC filtered")
# Avg. Silhouett
fviz_nbclust(data_scale, 
             FUNcluster = hcut,   # hierarchical clustering
             method = "silhouette", # Avg. Silhouette
             k.max = 12             # max number of clusters
) +labs(title="Avg.Silhouette Method for HC of KIRC filtered") 
#######################################################################
dd <- dist(data_scale, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc, hang = -2, cex = 0.6)#abline(h=80, col="red")
hc_2<-cutree(hc,k=2)
hc_3<-cutree(hc,k=3)
#HC實驗
result<-NULL
for (i in c(1:10)){
  x<-cutree(hc,k=i)
  result<-cbind(result,x)
  colnames(result)[i]<-paste("hc",as.character(i),sep="_")
}
write.csv(result,"C:/Users/syuan/Desktop/Master/kirc/kirc_hc_cutree.csv")
#######################################################################
HC 分群結果合併及summary
#######################################################################
#結果合併
clustering_result<-cbind(data,hc_2,hc_3)
write.csv(clustering_result,"C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filtered_clustering_result.csv")
clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filtered_clustering_result.csv")
class(clustering_result$hc_2)
#結果summary
x1<-c(length(which(clustering_result$hc_2==1)),
      length(which(clustering_result$hc_2==2)),
      length(which(clustering_result$hc_2==3)))
  sum1<-c(x1,sum(x1))
  x2<-c(length(which(clustering_result$hc_3==1)),
        length(which(clustering_result$hc_3==2)),
        length(which(clustering_result$hc_3==3)))
  sum2<-c(x2,sum(x2))
  summary<-cbind(sum1,sum2)
  write.csv(summary,"C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filter_clustering_summary.csv")
#打卡打卡20200313
#################################################################################
分群結果視覺化
#################################################################################
kirc_filtered_clustering_result<-read.csv("C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filtered_clustering_result.csv")
p1<-ggplot()+geom_point(data=kirc_filtered_clustering_result,
                        aes(x=TWIST1,y=TWIST2,color=as.factor(hc_2)))+ggtitle("H Clustering")
p2<-ggplot()+geom_point(data=kirc_filtered_clustering_result,
                    aes(x=TWIST1,y=TWIST2,color=as.factor(hc_2)))+ggtitle("H Clustering")
p3<-ggplot()+geom_point(data=kirc_filtered_clustering_result,
                        aes(x=TIMP1,y=TIMP3,color=as.factor(hc_2)))+ggtitle("H Clustering")
par(mfrow=c(1,2))
plot(p1)
plot(p2)
#3d plot N
png("C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filter_3Dplot_N.png",width = 1800, height = 1200, units = "px", pointsize = 24,
    bg = "white")
head(kirc_filtered_clustering_result)
colors <- c("#1e90ff", "#ff4500")
colors <- colors[as.factor(kirc_filtered_clustering_result$hc_2)]
s3d<-scatterplot3d(kirc_filtered_clustering_result$TIMP1,
              kirc_filtered_clustering_result$TIMP3,
              kirc_filtered_clustering_result$TXNIP,
              xlab="TIMP1", ylab="TIMP3", zlab="TXNIP",pch = 16, color=colors)
legend("right", legend = as.character(unique(kirc_filtered_clustering_result$hc_2)),
       col =  c("#1e90ff", "#ff4500"),pch = 16,title="kirc filter clustering")
dev.off()
#3d plot P
png("C:/Users/syuan/Desktop/Master/kirc/Clustering_filtered/kirc_filter_3Dplot_2P1N.png",width = 1800, height = 1200, units = "px", pointsize = 24,
    bg = "white")
head(kirc_filtered_clustering_result)
colors <- c("#1e90ff", "#ff4500")
colors <- colors[as.factor(kirc_filtered_clustering_result$hc_2)]
s3d<-scatterplot3d(kirc_filtered_clustering_result$TWIST1,
                   kirc_filtered_clustering_result$TWIST2,
                   kirc_filtered_clustering_result$DICER1,
                   xlab="TWIST1", ylab="TWIST2", zlab="DICER1",pch = 16, color=colors)
legend("right", legend = as.character(unique(kirc_filtered_clustering_result$hc_2)),
       col =  c("#1e90ff", "#ff4500"),pch = 16,title="kirc filter clustering")
dev.off()
###########################################################
3D Plot example
###########################################################
data(iris)
head(iris)
colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(iris$Species)]
s3d <- scatterplot3d(iris[,1:3], pch = 16, color=colors)
legend("bottom", legend = levels(iris$Species),
       col =  c("#999999", "#E69F00", "#56B4E9"), 
       pch = c(16, 17, 18), 
       inset = -0.25, xpd = TRUE, horiz = TRUE)
###############################################################
heatmap?
###############################################################
heatmap(as.matrix(data[,8:43]),col=data$hc_2)
heatmap(as.matrix(data[,8:44]), scale = "none", col=hc_2, RowSideColors = rep(c("blue", "pink"), each=16), 
      #ColSideColors = c(rep("purple", 5), rep("orange", 6))
      )
