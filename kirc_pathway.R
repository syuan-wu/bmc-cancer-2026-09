library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(msigdbr)
library(enrichplot)
#############################################################################
kirc_filtered<-read.csv("C:/Users/syuan/Desktop/Master/kirc/kirc_filtered.csv")
kirc_filtered<-kirc_filtered[,-1]
gene_df <- bitr(kirc_filtered$external_gene_name, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
colnames(kirc_filtered)[2]<-"SYMBOL"
genelist<-merge(gene_df,kirc_filtered[,c(2,12)], by="SYMBOL")
genelist<-as.numeric(as.character(genelist$FC_value))
names(genelist)<-gene_df$ENTREZID
head(genelist)
genelist<-sort(genelist,decreasing = T)
head(genelist)
#############################################################################
kegg<-enrichKEGG(names(genelist), organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,
                 qvalueCutoff = 0.2)
write.csv(summary(kegg),"C:/Users/syuan/Desktop/Master/kirc/kirc_kegg_summary.csv",row.names =F)
dotplot(kegg,showCategory=20)
barplot(kegg,showCategory=20)
png("C:/Users/syuan/Desktop/Master/kirc/kirc_kegg_network.png",units="px",height=3000,width=4000,res=250)
cnetplot(kegg,categorySize="pvalue",foldChange = genelist)
dev.off()
which(names(genelist)=="96459")
gene_df[which(gene_df$ENTREZID == "96459"),]
FNIP1 - 頝idney tumor supression????
  #############################################################################
ggo<-groupGO(gene = names(genelist), OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
eng<-enrichGO(names(genelist), OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, 
              pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, 
              maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(summary(eng),"C:/Users/syuan/Desktop/Master/kirc/kirc_go_summary.csv",row.names =F)
barplot(eng,drop=TRUE,showCategory=20)
dotplot(eng,title="enrichGo CC dotplot")

#############################################################################
heatplot
#############################################################################
data(geneList, package="DOSE")
gene <- names(geneList)
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
de <- names(geneList)[1:60]
kk <- enrichKEGG(de, organism  ='hsa',  pvalueCutoff = 0.05)
final_productKEGG_Overrep <- heatplot(kegg, foldChange = genelist, showCategory= 20)+ ggplot2::coord_flip()
final_productKEGG_Overrep
