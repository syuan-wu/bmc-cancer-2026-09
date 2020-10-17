auc_result<-NULL
n_result<-NULL
for(i in c(1:1000)){
  n<-sample(397,9)
  z<-paste( colnames(df_n)[n],collapse = '+' )
  z<-paste0("response","~",z)
  #z
  md<-glm(formula=z,df_n,family="binomial", na.action=na.exclude)
  prediction<-predict(md,df_n)
  pred <- prediction(prediction,df_n$response)
  auc <- performance(pred, "auc")@y.values
  auc_result<-c(auc_result,auc)
  n_result<-c(n_result,list(n))
}
which.max(auc_result)
auc_result[which.max(auc_result)]
name_9<-colnames(df_n)[n_result[[which.max(auc_result)]]]
#取8
auc_result<-NULL
n_result<-NULL
for(i in c(1:1000)){
  n<-sample(397,8)
  z<-paste( colnames(df_n)[n],collapse = '+' )
  z<-paste0("response","~",z)
  #z
  md<-glm(formula=z,df_n,family="binomial", na.action=na.exclude)
  prediction<-predict(md,df_n)
  pred <- prediction(prediction,df_n$response)
  auc <- performance(pred, "auc")@y.values
  auc_result<-c(auc_result,auc)
  n_result<-c(n_result,list(n))
}
which.max(auc_result)
auc_result[which.max(auc_result)]
name_8<-colnames(df_n)[n_result[[which.max(auc_result)]]]

final<-merge(data.frame(name=name_10),data.frame(name=name_9),by="name")
final<-merge(final,data.frame(name=name_8),by="name")
final$name
#CCDC71L+CCL7+SNAI2
md<-glm(response~CCDC71L+CCL7+SNAI2,df_n,family="binomial", na.action=na.exclude)
prediction<-predict(md,df_n)
pred <- prediction(prediction,df_n$response)
performance(pred, "auc")@y.values
#######################################################
Combination of gene set and ROC of logistic model
#######################################################

##################################################################
best<-data.frame(result[10])
best[,635]
colnames( df[,as.numeric(best[,635])] )
#"SERPINE1""ADAM12""COL11A1""SCG2""THBS2""LUM""ECM1""WNT5A""FBN2"
z<-paste(colnames( df[,as.numeric(best[,635])] ),collapse = '+')
z<-paste0("response","~",z)
z
md<-glm(formula=z,df,family="binomial", na.action=na.exclude)
summary(md)
coef(md)
prediction<-predict(md,df)
pred <- prediction(prediction,df$response)
auc <- performance(pred, "auc")@y.values
auc
#plot
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
png("C:/Users/syuan/Desktop/Master/kirc/kirc_hc_31/roc.png",width = 600, height = 500, units = "px", pointsize = 24,
    bg = "white")
plot(perf, col = rainbow(7), main = "ROC curve", xlab = "1-Specificity(FPR)", ylab = "Sensitivity(TPR)")
#AUC = 0.5
abline(0, 1)
#實際AUC值
text(0.8, 0.5, auc)
dev.off()
#######################################################
Combination of gene set and ROC of logistic model-stepwise
#######################################################
df_n<-data.frame(apply(df_n,2,as.numeric))
#stepwise selection
z<-paste(colnames(df_n)[1:397],collapse = '+')
z<-paste0("response","~",z)
full <- glm(formula=z,data=df_n, family="binomial", na.action=na.exclude)
null <- glm(formula=response~1,data=df_n, family="binomial", na.action=na.exclude)
#從null
stepwise_model_n<-step(null, scope = list(upper=full), direction="both")
summary(stepwise_model_n)
coef(stepwise_model_n)
model_n<-glm(formula = response ~ EGFR+TXNIP+TWIST1+AKAP12+CTNNB1 , 
             family = "binomial", data = premodel, na.action = na.exclude)
#從full
stepwise_model_f<-step(full, scope = list(upper=full), direction="both")  
summary(stepwise_model_f)
coef(stepwise_model_f)
model_f<-glm(formula = response ~corresult+m1m0_FC+SNAI2+TIMP2+TP63+AKAP12+TIMP4+CTNNB1 , family = "binomial", data = premodel, na.action = na.exclude)
#######################################################
PCA
#######################################################
z<-paste(colnames(df_n)[1:397],collapse = '+')
z<-paste0("response","~",z)
z
pca<-prcomp(~ACKR3+ACSM4+ACTN1+ADA+ADAM12+ADAM18+ADAM19+ADAMTS12+ADAMTS14+ADAMTS2+ADAMTS3+ADAMTS4+ADORA2B+ADRA1B+ADRA2C+AEBP1+AHRR+ALDH1L2+ANGPT2+ANGPTL2+ANKLE2+ANKRD13B+AP000781.2+AP2A1+APCDD1+APCDD1L+APLN+APLP1+ARHGAP44+ARL4A+ARL4C+ARNTL2+ATP13A2+AVEN+AXL+B4GALNT1+BDKRB2+BEND6+BEST4+BHLHE22+BMP1+BYSL+C1orf21+C1orf216+C1QTNF1+C1QTNF6+C1R+C1S+C3orf80+C8orf34+CACNG8+CAMKK2+CAV1+CBLN4+CBX2+CCDC71L+CCIN+CCL7+CCN4+CCN5+CD163L1+CD248+CD276+CDC7+CDK2AP1+CERCAM+CFH+CFHR3+CHPF2+CHRDL2+CHST1+CHST14+CHST15+CHST3+CHSY3+CKAP4+CLDN11+CNN2+CNTNAP1+COL11A1+COL12A1+COL14A1+COL1A1+COL1A2+COL22A1+COL24A1+COL3A1+COL5A1+COL5A2+COL5A3+COL6A2+COL6A3+COLGALT1+CPE+CPXM1+CREG2+CRHR1+CTHRC1+CXCL5+CYP1B1+DAPK2+DCUN1D3+DENND2A+DIMT1+DIO2+DISP3+DLGAP4+DLX1+DLX2+DMRT3+DNAJB6+DPYSL3+ECM1+EDIL3+EFEMP1+EGFL6+EGFLAM+EGR2+EHD2+ELK1+ELOVL2+EN1+ENC1+EPGN+EPO+FADS2+FAM126A+FAM171A2+FAM78B+FAM83D+FAP+FBLN1+FBN1+FBN2+FGF7+FHL2+FJX1+FKBP10+FKBP1A+FLNC+FLRT2+FN1+FNDC1+FOLH1+FOSL1+FOXL1+FSCN1+FSTL1+FZD10+GAL3ST4+GALNT2+GAP43+GATA6+GBX2+GLP2R+GNA12+GPR85+GPX8+GRIK4+GRIN2A+GRIN2D+GRM2+HAS2+HIVEP3+HLX+HOXD11+HOXD12+HOXD13+HSD11B1+HTRA1+IBSP+IGFBP3+IGFL2+IL1R2+IL21R+IL2RA+IL6+IMPDH1+INHBA+INHBB+ITGA11+ITGA5+ITGBL1+JPH2+KCNB1+KCNG1+KCNJ2+KCNK17+KCNK2+KCNS3+KIAA1211+KIF4A+KIRREL1+KLF17+KLHDC8A+KPNA2+LEF1+LGI2+LGI4+LHX6+LIF+LINGO1+LMCD1+LMO1+LOXL2+LPAR4+LRFN1+LRRC15+LRRC4C+LTBP2+LUM+MAP1A+MAP3K7CL+MCHR1+MDFI+MEDAG+MEX3B+MFAP2+MGAT3+MICAL2+MLLT11+MLXIP+MMP10+MMP11+MMP13+MMP14+MMP16+MMP19+MMP2+MMP9+MMRN1+MPZL1+MSRB3+MTCP1+MXRA7+MYOZ3+NALCN+NAP1L1+NAV3+NCAM1+NCAM2+NCOR2+NPTX2+NRIP3+NT5DC2+NTM+NXPH4+OBP2A+OLFML2B+OPN1SW+OR1G1+OR51E1+OR51E2+OXTR+P3H1+P4HA3+PCBP3+PCDHB3+PCOLCE+PCSK1+PDE10A+PDGFRA+PDGFRL+PGF+PHACTR1+PI15+PIK3R6+PLCE1+PLEC+PLEKHG2+PLOD1+PLPPR2+PLPPR4+PLTP+PLXDC1+PLXNA4+PLXNB3+PMEPA1+PODNL1+POLR3G+POSTN+PPP2R2C+PREX2+PRRX1+PRRX2+PRSS35+PSD2+PTGIS+PTGS2+PTP4A3+PTPRZ1+PXN+PYGO1+RARG+RASD2+RASL10B+RFLNA+RGS4+RHBDL2+RIN1+RNASE10+RND3+RNF39+ROR2+RPE65+RTN2+RUNX1+RUNX2+S1PR2+SCG2+SDC3+SEMA3A+SEMA4B+SEMA7A+SEPTIN5+SEPTIN9+SERPINA12+SERPINA9+SERPINE1+SERPINF1+SERPINH1+SFRP2+SH2D4B+SH3GL1+SHC1+SHISAL1+SHOX2+SIX1+SLC12A5+SLC2A3+SLC30A3+SLC43A3+SLC52A3+SLC6A17+SLC6A9+SLFN11+SMIM3+SNAI1+SNAI2+SORCS3+SOX11+SPRED3+SPSB1+SRGAP1+SRPX2+SSC5D+ST3GAL2+STC1+STEAP1B+STMN2+SULF1+SULF2+SYBU+SYNPO+SYT4+TAFA3+TBC1D20+TBX15+TDG+TDO2+TEAD4+TENM2+TEX14+TGFB1+TGFB3+TGFBI+TGM2+THBS2+TMEM108+TMEM145+TMEM158+TMEM173+TMEM189+TMEM44+TPBGL+TPM4+TRAM2+TSHZ3+TUBB3+TUBB6+UCK2+ULBP2+UNC13A+UNC5A+VCAN+VOPP1+VWA3B+WNT5A+ZBTB7C+ZCCHC12+ZCCHC24+ZFHX4+ZNF114+ZNF365+ZNF469+ZNF730+ZP1, df_n, scale=T)
pca
plot(pca,
     type="line",
     main="")+abline(h=20, col="blue")
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
top100.pca.eigenvector <- pca$rotation[, 1:100]
top100.pca.eigenvector
first.pca <- top3.pca.eigenvector[, 1]   #  第一主成份
second.pca <- top3.pca.eigenvector[, 2]  #  第二主成份
third.pca <- top3.pca.eigenvector[, 3]   #  第三主成份
first.pca[order(first.pca, decreasing=FALSE)]  
# 使用dotchart，繪製主成份負荷圖
dotchart(first.pca[order(first.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC1",                      # 主標題
         xlab="Variable Loadings") 
version()
#########################################################
#c("KIF4A", "COL1A1", "VOPP1", "FNDC1", "MICAL2", "HOXD12", "SERPINE1", "HOXD13", "COL1A2", "C1orf21") 
ggboxplot( df_n ,x="response",y="C1orf21",palette="jama",color="response",title="")+stat_compare_means(method = "wilcox.test",comparisons = list(c("1", "0")) )
+theme_bw(base_size = 18)
#significance:KIF4A,HOXD12