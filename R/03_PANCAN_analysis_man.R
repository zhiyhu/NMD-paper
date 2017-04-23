#-----------------------#
# coding for manuscript #
#-----------------------#
# load packages-----
library(gplots)
library(scales)
library(ggplot2)
library(Rmisc)
library(plyr)
library(org.Hs.eg.db)
library(RColorBrewer)
library(biomaRt)
#--------------------------------------------------
# 1 NMD prevalant----------------
#--------------------------------------------------

# Distribution of mutations #----------------
#--------------------------------------------
df <-  data.frame(table(pancan_mut_all$cancer_abbr))
df$sample<-table(unique(pancan_mut_all[,c("sample","cancer_abbr")])[,"cancer_abbr"])
df$type="All mutations"
df2 <- data.frame(table(pancan_mut_all$cancer_abbr[pancan_mut_all$nmd.detail2==T]))
df2$sample<-table(unique(pancan_mut_all[,c("sample","cancer_abbr")])[,"cancer_abbr"])
df2$type ="NMD mutations"

df <- rbind(df,df2)
colnames(df)[1] <- "Cancers"
df <- df[order(df$Freq,decreasing = T),]
df$Cancers <- as.character(paste(df$Cancers))
df$Cancers <- paste(df$Cancers,"\nn=" ,df$sample,sep="")
df$Cancers <- factor(df$Cancers, levels = unique(df$Cancers))

# **Figure. Total number of mutations and the number of NMD mutation from 24 cancers. The top plot shows the number of all mutation and the bottom one shows NMD mutations.----
ggplot(data=df, aes(y = Freq, x = Cancers)) +
  geom_bar(position = "dodge", stat = "identity")+ facet_wrap(~type,ncol = 1,scales = "free") +#scale_fill_manual(values =colour_cancer)+
  xlab("")+ylab("Number of mutations")

ggplot(data=df[df$type=="NMD mutations",], aes(y = Freq, fill =Cancers,x="")) +
  geom_bar(width = 1, stat = "identity")+ coord_polar(theta = "y", start=0) +facet_grid(.~type) +#scale_fill_manual(values =colour_cancer)+
  xlab("")+ylab("Number of mutations")
ggplot(data=df[df$type=="All mutations",], aes(y = Freq, fill = Cancers,x=Cancers)) +geom_bar(width = 1, stat = "identity")

cancer_type$n_analysed_mut <- df$Freq[df$type=="All mutations"][match(cancer_type$Abbr,df$Cancers[df$type=="All mutations"])]
cancer_type$n_nmd_mut <- df$Freq[df$type=="NMD mutations"][match(cancer_type$Abbr,df$Cancers[df$type=="NMD mutations"])]
cancer_type$sample <- df$sample[df$type=="All mutations"][match(cancer_type$Abbr,df$Cancers[df$type=="All mutations"])]

# Prediction results #------------------------
#--------------------------------------------

# effect of prediction #-----------------------
# Convert your data into frequency matrix, then data.frame
#--------------------------------------------------
effect <- c("Silent","NM","MM","IFI","IFD","FSI","FSD")
##* all genes----
df2<-data.frame(table(pancan_mut_all[(pancan_mut_all$effect %in% effect),c("effect","nmd.detail2")]))
# Set columns ready for ggplot
df2$Freq<-as.factor(df2$Freq)
df2$effect <- as.character(df2$effect) #%>% as.numeric()
df2$nmd.detail2<- as.character(df2$nmd.detail2)# %>% as.numeric()

##* repeat with only tp53----
df3<-data.frame(table(pancan_mut_all[pancan_mut_all$gene=="TP53" & (pancan_mut_all$effect %in% effect),c("effect","nmd.detail2","cancer_abbr")]))
# Set columns ready for ggplot
df3$Freq<-as.factor(df3$Freq)
df3$effect <- as.character(df3$effect) #%>% as.numeric()
df3$nmd.detail2<- as.character(df3$nmd.detail2)# %>% as.numeric()
# **Figure. Novel annotation on NMD in TP53 across cancers.----
# plot a frequency heat plot using ggplot
ggplot(df3,aes(x=nmd.detail2,y=effect)) +geom_tile(aes(fill=as.numeric(Freq))) +
  geom_text(aes(label=Freq))+xlab(label = "NMD classification")+
  ylab(label = "TCGA variant classification")+ggtitle("TP53")+facet_grid(.~cancer_abbr)+ theme(legend.position="none")+scale_fill_gradient2()+ 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line =  element_blank(),panel.border = element_blank())+
  labs(fill="")

## table of cancers ##------
cancer_df <-data.frame (table(unique(pancan_mut_all[,c("sample","cancer_abbr")])[,"cancer_abbr"]))
colnames(cancer_df) <- c("cancer","n_sample")
cancer_df$mut <- table(pancan_mut_all$cancer_abbr)
cancer_df$NMD_mut <- table(pancan_mut_all$cancer_abbr[pancan_mut_all$nmd.detail2==T])
cancer_df$TP53_mut <- table(pancan_mut_all$cancer_abbr[pancan_mut_all$gene=="TP53"])
cancer_df$TP53_NMD_mut <- table(pancan_mut_all$cancer_abbr[pancan_mut_all$gene=="TP53"&pancan_mut_all$nmd.detail2==T])[match(cancer_df$cancer,names(table(pancan_mut_all$cancer_abbr[pancan_mut_all$gene=="TP53"&pancan_mut_all$nmd.detail2==T])))]

cancer_df$TP53_mut_sample <- table(unique(pancan_mut_all[pancan_mut_all$gene=="TP53",c("sample","cancer_abbr")])[,"cancer_abbr"])
cancer_df$TP53_NMD_mut_sample <- table(unique(pancan_mut_all[pancan_mut_all$gene=="TP53"&pancan_mut_all$nmd.detail2==T,c("sample","cancer_abbr")])[,"cancer_abbr"])[
  match(cancer_df$cancer,names(table(unique(pancan_mut_all[pancan_mut_all$gene=="TP53"&pancan_mut_all$nmd.detail2==T,c("sample","cancer_abbr")])[,"cancer_abbr"])))]
cancer_df$mut_in_TSG <-table(pancan_mut_all$cancer_abbr[grepl("TSG",pancan_mut_all$role_in_cancer)])[match(cancer_df$cancer,names(table(pancan_mut_all$cancer_abbr[grepl("TSG",pancan_mut_all$role_in_cancer)])))]
cancer_df$NMD_mut_in_TSG <-table(pancan_mut_all$cancer_abbr[grepl("TSG",pancan_mut_all$role_in_cancer)&pancan_mut_all$nmd.detail2==T])[match(cancer_df$cancer,names(table(pancan_mut_all$cancer_abbr[grepl("TSG",pancan_mut_all$role_in_cancer)&pancan_mut_all$nmd.detail2==T])))]
cancer_df$mut_in_TSG_sample <-table(unique(pancan_mut_all[grepl("TSG",pancan_mut_all$role_in_cancer),c("sample","cancer_abbr")])[,"cancer_abbr"])
cancer_df$NMD_mut_in_TSG_sample <-table(unique(pancan_mut_all[grepl("TSG",pancan_mut_all$role_in_cancer)&pancan_mut_all$nmd.detail2==T,c("sample","cancer_abbr")])[,"cancer_abbr"])[
  match(cancer_df$cancer,names(table(unique(pancan_mut_all[grepl("TSG",pancan_mut_all$role_in_cancer)&pancan_mut_all$nmd.detail2==T,c("sample","cancer_abbr")])[,"cancer_abbr"])))]

cancer_df$mut_in_cancer_gene <-table(pancan_mut_all$cancer[pancan_mut_all$driver==T])[match(cancer_df$cancer,names(table(pancan_mut_all$cancer[pancan_mut_all$driver==T])))]
cancer_df$NMD_mut_in_cancer_gene <-table(pancan_mut_all$cancer[pancan_mut_all$driver==T&pancan_mut_all$nmd.detail2==T])[match(cancer_df$cancer,names(table(pancan_mut_all$cancer[pancan_mut_all$driver==T&pancan_mut_all$nmd.detail2==T])))]
cancer_df$silent_mut <- table(pancan_mut_all$cancer[pancan_mut_all$effect=="Silent"])
cancer_df$silent_mut_in_TSG <- table(pancan_mut_all$cancer[grepl("TSG",pancan_mut_all$role_in_cancer)&pancan_mut_all$effect=="Silent"])[match(cancer_df$cancer,names(table(pancan_mut_all$cancer[grepl("TSG",pancan_mut_all$role_in_cancer)&pancan_mut_all$effect=="Silent"])))]

cancer_df$ratio_NMD_mut_mut <- cancer_df$NMD_mut/cancer_df$mut
cancer_df$ratio_TP53_NMD_mut_sample_TP53_mut_sample <- cancer_df$TP53_NMD_mut_sample/cancer_df$TP53_mut_sample
cancer_df$ratio_TP53_NMD_mut_sample_sample <- cancer_df$TP53_NMD_mut_sample/cancer_df$n_sample
cancer_df$ratio_NMD_mut_in_TSG_mut <- cancer_df$NMD_mut_in_TSG/cancer_df$mut
cancer_df$ratio_NMD_mut_in_TSG_sample_sample <-cancer_df$NMD_mut_in_TSG_sample/cancer_df$n_sample
cancer_df[is.na(cancer_df)] <- 0

ggplot(data=cancer_df,aes(x=as.numeric(mut),y=as.numeric(NMD_mut)))+geom_point()+
  geom_text(data=cancer_df, aes(x=as.numeric(mut)+10000,y=as.numeric(NMD_mut)+150, label = cancer), size=4) 

# **Percentage of samples with NMD-mutated TSG#------
cancer_df$cancer <- factor(cancer_df$cancer, levels = cancer_df[order(cancer_df$ratio_NMD_mut_in_TSG_sample_sample,decreasing = F),"cancer"])
sum(cancer_df$NMD_mut_in_TSG_sample)/sum(cancer_df$n_sample)
# [1] 0.2855663

cancer_df$cancer <- factor(cancer_df$cancer, levels = cancer_df[order(cancer_df$ratio_TP53_NMD_mut_sample_sample,decreasing = T),"cancer"])
ggplot(data=cancer_df[!is.na(cancer_df$ratio_TP53_NMD_mut_sample_sample),],
       aes(x=cancer,y =as.numeric(ratio_TP53_NMD_mut_sample_sample),fill=cancer))+
  geom_bar(stat="identity", position=position_dodge())+theme(legend.position="none")+
  labs(x="Cancers",y="Number of samples with NMD-mutated TP53/samples per cancer")

fisher.test(cancer_df[,c("mut","NMD_mut")])
#--------------------------------------------------

## table of samples ##------
sample_df <- data.frame(unique(pancan_mut_all[,c("sample","cancer_abbr")]))
colnames(sample_df)[2]<-"cancer"

tb <- table(pancan_mut_all$sample)
sample_df$mut <- tb[match(sample_df$sample,names(tb))]
tb <- table(pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T])
sample_df$NMD_mut <- tb[match(sample_df$sample,names(tb))]
tb <- table(pancan_mut_all$sample[pancan_mut_all$new.annotation=="NMD-escape"])
sample_df$NMD_mut <- tb[match(sample_df$sample,names(tb))]
tb <- table(pancan_mut_all$sample[pancan_mut_all$gene=="TP53"])
sample_df$TP53_mut <- tb[match(sample_df$sample,names(tb))]
tb <- table(pancan_mut_all$sample[pancan_mut_all$gene=="TP53"&pancan_mut_all$nmd.detail2==T])
sample_df$TP53_NMD_mut <- tb[match(sample_df$sample,names(tb))]

tb <- table(pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T&grepl("TSG",pancan_mut_all$role_in_cancer)])
sample_df$TSG_NMD_mut <- tb[match(sample_df$sample,names(tb))]

tb <- table(pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T&pancan_mut_all$driver==T])
sample_df$driver_NMD_mut <- tb[match(sample_df$sample,names(tb))]

sample_df[is.na(sample_df)] <- 0

sample_df$nmf_cluster <- nmf_clust$cluster[match(sample_df$sample,nmf_clust$SampleName)]

# number of substitution and indel----------
sample_df2 <- table(pancan_mut_all[,c("sample","effect")])
sample_df2 <-(matrix(sample_df2,ncol = 16,byrow = T))
rownames(sample_df2) <- rownames((table(pancan_mut_all[,c("sample","effect")])))
colnames(sample_df2) <- colnames((table(pancan_mut_all[,c("sample","effect")])))
sample_df2 <- as.data.frame(sample_df2)
sample_df2$sample <- rownames(sample_df2)
sample_df2$cancer <- sample_df$cancer[match(sample_df2$sample,sample_df$sample)]
sample_df2[rownames(sample_df2) %in% sample_df$sample[sample_df$cancer=="COAD"],]
cor(sample_df2[,1:16])
heatmap(cor(sample_df2[,1:16]))
sample_df2$mut_PLOD1 <- NA
# sample_df2$mut_PLOD1 <- pancan_mut_all[pancan_mut_all$gene=="POLD1","effect"]
sample_df2$TCGA.identifier <- substr(sample_df2$sample,start = 1,stop = 12)
sample_df2 <- merge(sample_df,sample_df2,by="sample")


# use mononucleotide_and_dinucleotide_marker_panel_analysis_status---------

sample_df2$MSS <- NA
sample_df2$MSS[match(clin_coad$TCGA.identifier,sample_df2$TCGA.identifier)] <-as.character (clin_coad$mononucleotide_and_dinucleotide_marker_panel_analysis_status)
sample_df2$MSS[match(clin_esca$TCGA.identifier,sample_df2$TCGA.identifier)] <- as.character (clin_esca$mononucleotide_and_dinucleotide_marker_panel_analysis_status)
sample_df2$MSS[match(clin_stad$TCGA.identifier,sample_df2$TCGA.identifier)] <- as.character (clin_stad$mononucleotide_and_dinucleotide_marker_panel_analysis_status)
sample_df2$MSS[match(clin_ucec$TCGA.identifier,sample_df2$TCGA.identifier)] <- as.character (clin_ucec$mononucleotide_and_dinucleotide_marker_panel_analysis_status)
table(sample_df2[,c("cancer.x","MSS")])
sample_mss<-sample_df2[!is.na(sample_df2$MSS),]

sample_mss$mut_PLOD1 <- F
sample_mss$mut_POLE <- F
sample_mss$mut_PLOD1[sample_mss$sample%in%pancan_mut_all$sample[pancan_mut_all$gene=="PLOD1"]] <- T
sample_mss$mut_POLE[sample_mss$sample%in%pancan_mut_all$sample[pancan_mut_all$gene=="POLE"]] <- T

# plot the relationship of microsatelite stability, number of mutations and NMD mutations------

ggplot(data=sample_mss,aes(x=as.numeric(mut),y=as.numeric(NMD_mut),col=as.factor(MSS),shape=paste(mut_PLOD1,mut_POLE)))+geom_point(alpha=0.6,size=2)+
  scale_x_log10() +
  scale_y_log10()+geom_abline(slope = 1,intercept = 0,col="grey",lty=2)+labs(x="Number of mutations in one sample",y="Number of NMD mutations in one sample")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+theme_bw()+  theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.key = element_rect( color = "white")) +theme(legend.position="bottom")+facet_wrap(~cancer.x,ncol = 2)

ggplot(data=sample_mss,aes(x=as.numeric(NMD_mut),y=as.numeric(IFI),col=as.factor(MSS)))+geom_point(alpha=0.6,size=1)+
  scale_x_log10() +
  scale_y_log10()+facet_wrap(~cancer.x,ncol = 2)
ggplot(data=sample_mss,aes(y=as.numeric(mut),x=cancer.y,fill=paste(as.factor(MSS),mut_PLOD1|mut_POLE)))+geom_boxplot(alpha=0.6)+scale_y_log10()
  scale_x_log10() +
  +geom_abline(slope = 1,intercept = 0,col="grey",lty=2)+labs(x="Number of mutations in one sample",y="Number of NMD mutations in one sample")+

# Ratio of samples with NMD-mutated TSGs across cancers-----
mean=1943/6924
df = data.frame(table(sample_df$cancer[sample_df$TSG_NMD_mut >0])/table(sample_df$cancer[sample_df$mut >5]))
df <- df[order(df$Freq,decreasing = T),]
df$Var1 <- factor(df$Var1,levels=df$Var1[order(df$Freq,decreasing = T)])    
ggplot(df, aes(x=Var1,y=Freq))+geom_bar(stat="identity")+coord_flip()+geom_hline(yintercept=mean,col="red",lty = 2)+labs(x="Cancers",y="Ratio of samples with NMD-mutated TSGs")
# linear regression resultls--------
lr <- lm(NMD_mut~mut+cancer,data = sample_df)
summary(lr) 
cor.rls.pvalue <- rep(0,24)
cor.rls.cor<- rep(0,24)
for(i in 1:24)
{
  cor.rls.pvalue[i] <- cor.test( x = sample_df[sample_df$cancer==cancer_type$Abbr[i],c("mut")],y= sample_df[sample_df$cancer==cancer_type$Abbr[i],c("NMD_mut")])$p.value
  cor.rls.cor[i] <- cor.test( x = sample_df[sample_df$cancer==cancer_type$Abbr[i],c("mut")],y= sample_df[sample_df$cancer==cancer_type$Abbr[i],c("NMD_mut")])$estimate  
}
View(data.frame(cancer=cancer_type$Abbr,pvalue=cor.rls.pvalue,pearson=cor.rls.cor,p.adj=p.adjust(cor.rls.pvalue,method = "BH")))
# Call:
#   glm(formula = NMD_mut ~ mut + cancer, data = sample_df)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -306.08    -2.72    -0.50     1.20   475.97  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.240e+00  2.104e+00   0.589  0.55566    
# mut          4.518e-02  3.484e-04 129.682  < 2e-16 ***
#   cancerBLCA  -3.126e-01  2.333e+00  -0.134  0.89343    
# cancerBRCA   5.794e-01  2.199e+00   0.264  0.79214    
# cancerCESC  -2.420e+00  2.549e+00  -0.949  0.34253    
# cancerCOAD   1.108e+01  2.434e+00   4.551 5.42e-06 ***
#   cancerDLBC  -9.912e+00  3.580e+00  -2.769  0.00564 ** 
#   cancerESCA  -6.094e-01  2.571e+00  -0.237  0.81267    
# cancerGBM   -1.835e+00  2.388e+00  -0.768  0.44226    
# cancerHNSC  -4.003e-01  2.282e+00  -0.175  0.86080    
# cancerKIRC   6.940e+00  2.513e+00   2.762  0.00575 ** 
#   cancerKIRP   7.954e-01  2.611e+00   0.305  0.76068    
# cancerLGG   -1.795e+00  2.281e+00  -0.787  0.43142    
# cancerLIHC   4.010e+00  2.346e+00   1.710  0.08737 .  
# cancerLUAD  -4.579e-01  2.274e+00  -0.201  0.84039    
# cancerLUSC  -3.181e+00  2.587e+00  -1.230  0.21881    
# cancerOV    -2.910e+00  2.346e+00  -1.240  0.21488    
# cancerPAAD  -3.576e+00  2.666e+00  -1.341  0.17981    
# cancerPRAD  -1.336e+00  2.374e+00  -0.563  0.57369    
# cancerREAD   2.951e+00  2.810e+00   1.050  0.29364    
# cancerSARC  -2.493e+00  2.445e+00  -1.020  0.30799    
# cancerSKCM  -1.583e+01  2.308e+00  -6.860 7.40e-12 ***
#   cancerSTAD   1.839e+01  2.415e+00   7.613 2.99e-14 ***
#   cancerTHCA  -1.279e+00  2.286e+00  -0.560  0.57570    
# cancerUCEC   1.203e+00  2.468e+00   0.487  0.62593    

p.adjust(summary(lr)$coefficients[,4],method = "BH")

# histological_type------
table(sample_df$histological_type)
# Kidney Clear Cell Renal Carcinoma 
# 212 
sample_df$histological_type <-KIRC_clinicalMatrix$neoplasm_histologic_grade[match(sample_df$sample,KIRC_clinicalMatrix$sampleID)]
table(sample_df$histological_type)
# G1 G2 G3 G4 GX 
# 1  8 95 80 24  4 
ggplot(data=sample_df[sample_df$cancer=="kidney clear cell carcinoma",],aes(x=as.numeric(mut),y=as.numeric(NMD_mut),col=as.factor(histological_type)))+geom_point(alpha=0.4)+
  scale_x_log10() +
  scale_y_log10()+geom_abline(slope = 1,intercept = 0,col="grey") +theme(legend.position="right")+labs(x="Number of mutations in one sample",y="Number of NMD mutations in one sample",col="histological type")

# good transfer from symbol to id------
SYMBOL2EG <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(SYMBOL2EG)
# Convert to a list
SYMBOL2EG <- as.list(SYMBOL2EG[mapped_genes])
gene.universe.id <- SYMBOL2EG[match(names(gene.universe),names(SYMBOL2EG))]
gene.universe.id <- sapply(gene.universe.id,function(x) return(x[1])) # map between gene id to symbor

gene_go <- names(gene.universe.id)[match(gene_go,gene.universe.id)]
gene_go <- na.omit(gene_go)

# gene.universe[match(gene_go,names(gene.universe))]

## which genes are more affected by NMD

# these genes have more than 1 mutations in 
gene.universe<- U_pan2$Z_score
names(gene.universe) <- U_pan2$gene
GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe, geneSel = function(x) x<0.2, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")

table <- GenTable(GOdata, classicFisher = resultFisher,weightfisher = resultweight, topNodes = 500 ,numChar=100)
# write.csv(table,file = "pan_cancer_results/twoBackgroundDiff_GOEA.csv")               
showSigOfNodes(GOdata, score(resultweight), firstSigNodes = 5, useInfo = 'all')

## Which genes are more frequently NMD mutated?------
write.csv(U_pan, file = "panCan12_results/MWW_table_pancancer20160917.csv",row.names = F)
ggplot(aes(x=n_nmd_mut, y = Z_score, col=as.factor(n_nmd_mut)),data=U_pan2)+geom_point()+
  theme(legend.position="none")+
  labs(x="Number of NMD mutations",y= "Z score (relative expression levels)")+
  geom_text(aes(label=as.character(paste(gene,cancer))),hjust=0, vjust=0,size=2,data = U_pan2[U_pan2$n_nmd_mut>10,])


## frequency of NMD mutations in 24 cancer types-----------
#--------------------------------------------------
pancan_mut_all$nmd.detail2 <- as.character(pancan_mut_all$nmd.detail2)
pancan_mut_all$nmd.detail2[is.na(pancan_mut_all$nmd.detail2)] <-  "unclassified"
pancan_mut_all$nmd.detail2 <- as.factor(pancan_mut_all$nmd.detail2)
charts.data=  data.frame(cancer=rep(cancer_type$Cancer,each=3),NMD = rep(c("non-NMD mutation","NMD mutation","unclassified"),nrow(cancer_type)),
                         counts = as.numeric(table(pancan_mut_all[,c("nmd.detail2","cancer")])),
                         percentage = as.numeric(table(pancan_mut_all[,c("nmd.detail2","cancer")]))/rep(table(pancan_mut_all$cancer),each=3))

n.samples <- unique(pancan_mut_all[,c("sample","cancer")])
n.samples <- table(n.samples$cancer) # number of samples per cancer type

charts.data$n.sample <- rep(n.samples,each=3) # number of samples per cancer type
charts.data$normalised.count <- charts.data$counts/charts.data$n.sample # normalised the number of mutations by sample numbers

# *plot before normalisation----
charts.data <- ddply(charts.data, .(cancer), 
                     transform, pos = cumsum(percentage) - (0.5 * percentage))

ggplot() + geom_bar(aes(y = percentage, x = cancer, fill = NMD), data = charts.data, stat="identity")+ 
  geom_text(data=charts.data, aes(x = cancer, y = pos, label = counts), size=4) + 
  theme(legend.position = "bottom",legend.direction="horizontal",legend.title = element_blank())

# *plot after normalisation----
charts.data <- ddply(charts.data, .(cancer), 
                     transform, pos = cumsum(normalised.count) - (0.5 * normalised.count))
ggplot() + geom_bar(aes(y = normalised.count, x = cancer, fill = NMD), data = charts.data, stat="identity")+ 
  geom_text(data=charts.data, aes(x = cancer, y = pos, label = round(normalised.count,1)), size=4) + 
  theme(legend.position = "bottom",legend.direction="horizontal",legend.title = element_blank()) +
  labs(x="Cancer",y="Number of mutations/number of samples")+ 
  scale_x_discrete(breaks = cancer_type$Cancer, labels = paste(cancer_type$Abbr,n.samples,sep="\nn = "))

chisq.test(as.matrix(table(pancan_mut_all[,c("nmd.detail2","cancer")])))

# **Cancer-specific signatures/heatmap-----
nmd_signature <- table(pancan_mut_all[pancan_mut_all$gene %in% cancergene$Gene.Symbol & pancan_mut_all$nmd.detail2==T,c("gene","sample")])
head(colnames(nmd_signature))
nmd_signature[nmd_signature==1] <- 5

nonnmd_signature <- table(pancan_mut_all[pancan_mut_all$gene %in% cancergene$Gene.Symbol & pancan_mut_all$nmd.detail2==F,c("gene","sample")])
nonnmd_signature <- nonnmd_signature[match(rownames(nmd_signature),rownames(nonnmd_signature)),
                                     match(colnames(nmd_signature),colnames(nonnmd_signature))]
nmd_signature[nonnmd_signature==1] <- 3
rm(nonnmd_signature)

nonnmd_signature <- table(pancan_mut_all[pancan_mut_all$gene %in% cancergene$Gene.Symbol & pancan_mut_all$nmd.detail2=="unclassified",c("gene","sample")])
nonnmd_signature <- nonnmd_signature[match(rownames(nmd_signature),rownames(nonnmd_signature)),
                                     match(colnames(nmd_signature),colnames(nonnmd_signature))]
nmd_signature[nonnmd_signature==1] <- 1

heatmap_cancer <- pancan_mut_all$cancer_abbr[match(colnames(nmd_signature),pancan_mut_all$sample)]
sum(rowSums(nmd_signature==5)>30)
nmd_signature2 <-  nmd_signature[,heatmap_cancer==cancer_type$Abbr[1]]
for(i in 2:nrow(cancer_type))
{
  nmd_signature2 <- cbind(nmd_signature2,nmd_signature[,heatmap_cancer==cancer_type$Abbr[i]])
}

nmd_signature2<- nmd_signature2[which(rowSums(nmd_signature2==5)>30),which(colSums(nmd_signature2==5)>0)]
heatmap_cancer <- pancan_mut_all$cancer_abbr[match(colnames(nmd_signature2),pancan_mut_all$sample)]
head(nmd_signature2)

colfunc <- colorRampPalette(c("white", "grey", "#00BFC4","#F8766D"))
colSeperate <- table(heatmap_cancer )[match(unique(heatmap_cancer),names(table(heatmap_cancer )))]
for(i in 24:2)
{
  colSeperate[i] <- sum(colSeperate[1:(i-1)])
}
# pancan_col <- c("#F8766D" "#ED813E" "#DE8C00" "#CD9600" "#B79F00" "#9DA700" "#7CAE00" "#49B500" "#00BA38" "#00BE67" "#00C08B"
# "#00C1A9" "#00BFC4" "#00BBDC" "#00B4F0" "#00A9FF" "#619CFF" "#9F8CFF" "#C77CFF" "#E36EF6" "#F564E3" "#FF61CC"
# "#FF64B0" "#FF6C91")

# figure 3b------
column_col <- pancan_col[match(heatmap_cancer,cancer_type$Abbr)]
heatmap.2(nmd_signature2 ,col=colfunc(15), trace = "none",dendrogram="row",Colv = F,
          labCol = FALSE, colsep=   colSeperate[-1],sepcolor="black")
barplot(column_col )


# p<-ggplot(data = data.frame(x=c("x","y"),y=c(10,20)),aes(x=x,y=y,fill=x))+geom_bar(stat = "identity")
# pancan_col=unique(ggplot_build(p)$data[[1]]$colour)
# [1] "#F8766D" "#ED813E" "#DE8C00" "#CD9600" "#B79F00" "#9DA700" "#7CAE00" "#49B500" "#00BA38" "#00BE67" "#00C08B"
# [12] "#00C1A9" "#00BFC4" "#00BBDC" "#00B4F0" "#00A9FF" "#619CFF" "#9F8CFF" "#C77CFF" "#E36EF6" "#F564E3" "#FF61CC"
# [23] "#FF64B0" "#FF6C91"

# fill x  y PANEL group ymin ymax xmin xmax colour size linetype alpha
# 1 #F8766D 1 10     1     1    0   10 0.55 1.45     NA  0.5        1    NA
# 2 #00BFC4 2 20     1     2    0   20 1.55 2.45     NA  0.5        1    NA


## frequency of NMD mutations across cancers: boxplot-----------
charts.data= table(pancan_mut_all[,c("sample","nmd.detail2")])
charts.data <- data.frame(sample= rownames(charts.data),non_NMD=charts.data[,1],NMD=charts.data[,2])
charts.data$cancer <- pancan_mut_all[match(rownames(charts.data),pancan_mut_all$sample),"cancer"]
head(charts.data)

ggplot(charts.data,aes(x=cancer,y=non_NMD))+geom_boxplot()

#### table of genes ###-------
gene_df_add <- unique(pancan_mut_all[pancan_mut_all$nmd.detail2==T,c("cancer_abbr","gene")])
gene_df_add  <- table(table(gene_df_add $gene))
gene_df_add <- data.frame(gene_df_add )
gene_df_add $Freq<- as.numeric(gene_df_add $Freq)

# Figure. Distribution of NMD-mutated genes in cancers-----
ggplot(data = gene_df_add, aes(x = Var1,y=Freq))+ geom_bar(stat="identity") + labs(x="Number of cancer",y="Number of genes with NMD mutations in the given cancers")
# ggplot(data = gene_df[gene_df$Freq>10,], aes(x = Var1,y=Freq))+ geom_bar(stat="identity") + labs(x="Number of cancer",y="Number of genes with NMD mutations in the given cancers")
gene_df[gene_df$Var1,]
# For transformed scales, binwidth applies to the transformed data.

# The bins have constant width on the transformed scale.
gene_df <- table(pancan_mut_all[,c("cancer","gene")]$gene)
gene_df <- cbind(data.frame(table(pancan_mut_all[pancan_mut_all$nmd.detail2==T,c("cancer","gene")]$gene)),data.frame(table(gene_df$gene)))
gene_df <- gene_df[,c(1,2,4)]
colnames(gene_df) <- c("gene","NMD.mutation","NMD.mutated.cancer")
# ggplot(data = share_gene_cancer_table[1:500, ],aes(x=NMD.mutation,y=NMD.mutated.cancer)) +geom_point()

# number of nmd mutations in each gene in each cancer------
head(gene_df)
gene_df <- table(pancan_mut_all[pancan_mut_all$nmd.detail2==T,c("gene","cancer_abbr")])
gene_df <- as.matrix(gene_df)
gene_df <- cbind(gene_df,rowSums(gene_df))
gene_df <- data.frame(gene_df)
gene_df <- data.frame(gene=rownames(gene_df),gene_df)
gene_df[order(gene_df$BLCA,decreasing = T),][1:20,]
gene_df$NMD.mutated.cancer <- rowSums(gene_df[,2:25]>0)
gene_df$average_NMD <- rowMeans(gene_df[,2:25])
ggplot(data =gene_df,aes(x=NMD.mutated.cancer,y=average_NMD))+geom_point(pch=20,col=alpha("black",0.6))+
  geom_text(aes(label=gene),hjust=0, vjust=0,size=2,data = gene_df[gene_df$V25>100,])
gene_df$cancer_driver <- F
gene_df$cancer_driver[gene_df$gene %in% cancergene$Gene.Symbol] <- T
gene_df$role.in.cancer <- cancergene$Role.in.Cancer[match(gene_df$gene,cancergene$Gene.Symbol)]
gene_df$related.cancer <- cancergene$Tumour.Types.Somatic.[match(gene_df$gene,cancergene$Gene.Symbol)]

gene_df[order(gene_df$BRCA,decreasing = T),][1:20,]

write.csv(gene_df,"data/table/gene_df.csv")

# plot the expression of frequently NMD-mutated tumor genes-----
gene_df$gene[which(gene_df$BRCA>5&gene_df$cancer_driver==T)]


leader_nmd_gene <- data.frame(table(unique(pancan_mut_all[,c("sample","cancer_abbr")])[,"cancer_abbr"]))
colnames(leader_nmd_gene) <- c("cancer","n_sample")
leader_nmd_gene$leader <- NA
for(i in 1:nrow(leader_nmd_gene))
{
  m <- gene_df[gene_df[,colnames(gene_df)==leader_nmd_gene$cancer[i]]>leader_nmd_gene$n_sample[i]*0.05,c(1,which(colnames(gene_df)==leader_nmd_gene$cancer[i]),29)]
  if(nrow(m)>0)
  {
    m <- m[order(m[,colnames(m)==leader_nmd_gene$cancer[i]] ,decreasing = T),]
    m$mark[m$cancer_driver==T] <- paste(m$gene[m$cancer_driver==T],"*",sep="")
    m$mark[m$cancer_driver==F] <- paste(m$gene[m$cancer_driver==F])
    m$mark <- paste(m$mark,"(",m[,colnames(m)==leader_nmd_gene$cancer[i]],")",sep="")
    leader_nmd_gene$leader[i] <- paste(m$mark,collapse = ", ")
  }
}



write.csv(leader_nmd_gene,"data/table/leader_nmd_gene_cancers.csv",row.names = F)

gene_id = c("KDM6A","ARID1A","APC","PBRM1","ATRX","SMAD4","RB1","NF1","RNF43","PTEN")
NMDtarget <- pancan_mut_all[pancan_mut_all$gene %in% gene_id,] 
NMDtarget <- NMDtarget[NMDtarget$cancer_abbr%in%names(table(NMDtarget$cancer_abbr))[table(NMDtarget$cancer_abbr[NMDtarget$nmd.detail2==T])>10],]
ggplot(aes(x = interaction(nmd.detail2,cancer_abbr),fill=nmd.detail2,y=as.numeric(expression)), data=NMDtarget[(NMDtarget$nmd.detail2!="unclassified") & !is.na(NMDtarget$cancer) ,] ) + geom_boxplot(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0,alpha=0.7))+
  labs(x="Cancers", y = "Expression level",fill = "NMD mutations",col =  "NMD mutations",alpha="") +facet_wrap(~ gene,ncol = 2)+theme(legend.position="bottom")+scale_x_discrete(labels=rep(names(table(NMDtarget$cancer_abbr)),each=2))

gene_id <- data.frame(gene =c("KDM6A","ARID1A","APC","PBRM1","ATRX","SMAD4","RB1","NF1","RNF43","PTEN"),ratio = NA,p_value=NA)
gene_id <- data.frame(gene =c("NF1","TP53"),ratio = NA,p_value=NA)

for(i in 1:2)
{
  m <- pancan_mut_expression[pancan_mut_expression$gene %in% gene_id$gene[i],] 
  rls <- wilcox.test(x = m$quantile_expression[m$nmd.detail2==T],y = m$quantile_expression[m$nmd.detail2==F],alternative = "less")
  gene_id$ratio[i] <- median(m$quantile_expression[m$nmd.detail2==T],na.rm = T)/median(m$quantile_expression[m$nmd.detail2==F],na.rm = T)
  gene_id$p_value[i] <- rls$p.value
}
gene_id$p_value.adj <- p.adjust(gene_id$p_value,method="BH")

#--------------------------------------------------
#--------------------------------------------------
##------ Mutation-level analysis ------##
#--------------------------------------------------

## remove low expressed genes-----
sum(unique(pancan_mut_all$gene) %in% rownames(pancan_seq))
# 17544
length(unique(pancan_mut_all$gene) )
# 17544
pancan_mut_all <- pancan_mut_all[pancan_mut_all$gene %in% rownames(pancan_seq),]
length(unique(pancan_mut_all$gene) )
# [1] 17544
length(unique(pancan_mut_all$sample))
# [1] 3152

# Expression analysis------
# read file----
load(file = "data/results/pancan_mut_classfied_with_quantileexpression.RData")

## rank the genes by the wt expression median-----
pancan_mut_expression$wt_expression_rank <- NA
table(unique(pancan_mut_expression[pancan_mut_expression$nmd.detail2==T,c("gene","cancer")])$cancer)

pancan_mut_expression$wt_expression_rank <-NA
order <- order( pancan_mut_expression$wt_expression_median,decreasing = F)
n<- round(length(order)/10)+1
pancan_mut_expression$wt_expression_rank[order] <- c(rep(1:9,each=n),rep(10,length(order)-9*n))
rm(order)

pancan_mut_expression$wt_expression_rank <- as.factor(pancan_mut_expression$wt_expression_rank)
table(pancan_mut_expression$wt_expression_rank)
write.csv(file = "panCan17_results/number_genes_each_group_10",x = table(pancan_mut_expression[pancan_mut_expression$nmd.detail2==T,c("cancer","wt_expression_rank")]))


#-------------------------------------------
## Dependence on wild type expression-------
#-------------------------------------------
ggplot(data=pancan_mut_expression[(!is.na(pancan_mut_expression$wt_expression_rank))&
                                      (pancan_mut_expression$nmd.detail2%in%c(T,F)),], 
       aes(x=wt_expression_rank, y=quantile_expression,fill=nmd.detail2)) +
  geom_boxplot( ) +
  labs(y="Relative expression of variants",x = "Rank by median of wildtype expression")+
  theme(legend.position = "bottom",legend.direction="horizontal",legend.title = element_blank()) +
  scale_fill_discrete(labels=c("non-NMD mutation", "NMD mutation")) +
  facet_wrap(~ cancer,ncol = 4)

mean_group <-data.frame(cancer = rep(cancer_type$Cancer,each = 20),group = rep(1:10,24),NMD = rep(c(rep(F,10),rep(T,10)),12), mean = rep(0,240))
for(ct in 1:nrow(cancer_type))
{
  for(i in 1:10)
  {
    mean_group$mean[mean_group$cancer==cancer_type$Cancer[ct]&mean_group$group==i&mean_group$NMD==T] <-  median(na.rm=T,pancan_mut_expression$quantile_expression[pancan_mut_expression$cancer==cancer_type$Cancer[ct]&pancan_mut_expression$wt_expression_rank==i&pancan_mut_expression$nmd.detail2==T])
    mean_group$mean[mean_group$cancer==cancer_type$Cancer[ct]& mean_group$group==i&mean_group$NMD==F] <- 
      median(na.rm=T,pancan_mut_expression$quantile_expression[pancan_mut_expression$cancer==cancer_type$Cancer[ct]&pancan_mut_expression$wt_expression_rank==i&pancan_mut_expression$nmd.detail2==F])
    
  }
}
mean_group$nan.index <- F
mean_group$nan.index[which(is.na(mean_group$mean))] <- T
mean_group$mean[is.na(mean_group$mean) ] <- 0
mean_group$cancer_abbr <- cancer_type$Abbr[match(mean_group$cancer,cancer_type$Cancer)]
write.csv(mean_group,file="PANCAN_results/figure2a_table.csv",row.names = F)

# ~Comparing two Regression Lines-----
# Test for the slope
anova(lm(mean~group+cancer,data=mean_group[mean_group$nan.index==F,]),lm(mean~group+NMD+cancer,data=mean_group[mean_group$nan.index==F,]))

## linear regression

## global expression difference

ggplot(aes(x= cancer, y = expression, fill = nmd.detail2),
       data= pancan_mut[pancan_mut_expression$gene %in% gene_list_nmd,]) +geom_boxplot()
ggplot(aes(y=Z_score,x=cancer_type,fill=cancer_type),data=meta_U)+geom_boxplot()


##~ Figure 8. Predicted NMD mutations have lower expression than non-NMD mutations in frameshift indels and nonsense mutations-----
## difference in NMD FSI/D and non-NMD FSI/D
mut<-pancan_mut_expression[!is.na(pancan_mut_expression$nmd.detail2)&pancan_mut_expression$cancer_abbr!="UCEC",]


test.list <- data.frame(variant = rep(c("FSD","FSI","NM"),each=3))
test.list$pvalue <- NA

for(i in 1:3)
{
  
  rls <- wilcox.test(quantile_expression ~ have.ptc, 
                     data = mut[mut$effect==test.list$variant[i*3-2] &(mut$nmd.detail2==F),],
                     alternative = "greater")
  test.list$pvalue[i*3-2] <- rls$p.value
  rls <- wilcox.test(quantile_expression ~ nmd.detail2, 
                     data = rbind(mut[mut$effect==test.list$variant[i*3-1] &(mut$nmd.detail2==F)&(mut$have.ptc==F),],mut[mut$effect==test.list$variant[i] &(mut$nmd.detail2==T),]),
                     alternative = "greater")
  test.list$pvalue[i*3-1] <- rls$p.value
  rls <- wilcox.test(quantile_expression ~ nmd.detail2, 
                     data = mut[mut$effect==test.list$variant[i*3] &(mut$have.ptc==T),],
                     alternative = "greater")
  test.list$pvalue[i*3] <- rls$p.value
}
test.list$corrected.pvalue[!is.na(test.list$pvalue)] <- p.adjust(test.list$pvalue[!is.na(test.list$pvalue)],method = "BH")

pp<-ggplot(data = mut[c(grep("FS",mut$effect),grep("NM",mut$effect)),],aes(x = effect, y = quantile_expression)) + geom_boxplot(data = mut[c(grep("FS",mut$effect),grep("NM",mut$effect)),],aes( fill = interaction(have.ptc,nmd.detail2)))+
  labs(y="Relative expression of variants",x = "TCGA variant classification")+
  theme(legend.position = "bottom",legend.direction ="vertical") +
  scale_fill_discrete(name = "",labels=c("non-NMD mutation without PTC", "non-NMD mutation with PTC","NMD mutation with PTC")) +
  scale_x_discrete(labels=c("Frameshift deletion","Frameshift insertion","Nonsense mutation"))+ylim(c(-0.1,1.2))+
  
  # add aster
  label.df <- data.frame(Group = c("FSI", "FSD","NM"),
                         Value = c(1.1, 1.1,1.1))

pp+ geom_text(data = label.df, label = "***")

df1 <- data.frame(a = c(0.5, 0.5:1,1), b = c(1.02,1.1,1.02))
df2 <- data.frame(a = c(1, 1,2, 2), b = c(35, 36, 36, 35))
df3 <- data.frame(a = c(2, 2, 3, 3), b = c(24, 25, 25, 24))

pp+ geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1, y = 1.13, label = "*", size = 8) 
+
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 38, label = "**", size = 8) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 27, label = "n.s.", size = 8)

write.csv(x = test.list,file = "panCan12_results/wilcoxTest_oneside_diff_FS_NM.csv")
# test.list$corrected.pvalue <- round(test.list$corrected.pvalue)


#-------------------------------------------
# linear regression analysis---------
#-------------------------------------------
# gene expression=tumout type+gene-specific effect+nmd effect+noise
require(MASS)
require(robustbase)
library(caret)
set.seed(123)
inTrain <- createDataPartition(y =pancan_mut_expression$expression,
                               ## the outcome data are needed
                                 p = .75,
                                ## The percentage of data in the training set
                                 list = FALSE)
                               ## The format of the results
training <- pancan_mut_expression[ inTrain,]
testing  <- pancan_mut_expression[-as.numeric(inTrain),]
testing <- testing[testing$cancer_abbr!="UCEC",]
nmd.mod1 = lmrob(quantile_expression ~  nmd.detail2 + wt_expression_median +cancer_abbr + cnv,
                    data = training,fast.s.large.n = Inf)
normalised_expression <- predict(object = nmd.robmod1,newdata = testing[testing$cancer_abbr!="UCEC",])
cor(na.omit(cbind(normalised_expression,testing$expression[testing$cancer_abbr!="UCEC"])))
cor(na.omit(cbind(testing$wt_expression_median[testing$cancer_abbr!="UCEC"],testing$expression[testing$cancer_abbr!="UCEC"])))

cancer_coeff <- summary(nmd.mod1)$coefficient

testing$quantile_expression-testing$wt_expression_median-testing$cancer_abbr

#--------------------------------------------------
# TP53 ---------------------------------------------
#--------------------------------------------------

table(pancan_mut_expression$nmd.detail2)
# FALSE   TRUE 
# 796972  43788
gene_id = "TP53" # TP53
tp53 <- pancan_mut[pancan_mut_expression$cancer_abbr == "OV"& pancan_mut_expression$gene == gene_id,] #ovarian cancer
tp53 <- as.data.frame(tp53[tp53$gene==gene_id,c("sample","effect","nmd.detail2","expression","cancer_abbr","nmfc")])
tp53_wt <- unlist(c(get.wt.expression(gene = "TP53",cancer = "OV",nmf = 1),
                    get.wt.expression(gene = "TP53",cancer = "OV",nmf = 2),
                    get.wt.expression(gene = "TP53",cancer = "OV",nmf = 3)))
tp53_wt <- data.frame( sample = names( tp53_wt ),
                       effect= rep("WT",length(tp53_wt )),
                       nmd.detail2 = rep(F,length(tp53_wt )),
                       expression = unlist(tp53_wt) ,
                       cancer_abbr = rep("OV", length(tp53_wt )),
                       nmfc = c(rep(1,length(get.wt.expression(gene = "TP53",cancer = "OV",nmf = 1))),
                                rep(2,length(get.wt.expression(gene = "TP53",cancer = "OV",nmf = 2))),
                                rep(3,length(get.wt.expression(gene = "TP53",cancer = "OV",nmf = 3)))))
tp53 <- rbind(tp53,tp53_wt )
tp53$expression <- as.numeric(tp53$expression) 
save(tp53,file="PANCAN/results/Figure_TP53_OV_expression.RData")

# TP53 in all cancer types  
tp53 <- pancan_mut_expression[pancan_mut_expression$gene == "TP53",]
tp53$effect[grep("FS",tp53$effect)] <- "Frameshift indels"
tp53$effect[grep("IF",tp53$effect)] <- "Inframe indels"
tp53$cancer <- as.factor(tp53$cancer)
colnames(tp53)[colnames(tp53)=="expression"] <- "Expression"
tp53$Expression <- as.numeric(tp53$Expression)
tp53$effect[tp53$effect%in%c("NM","Frameshift indels")] <-"FS/NM"

wilcox.test(Expression~nmd.detail2,data = tp53, alternative = "greater")

#---------------------------------------
# ---------NF---------------------------
#---------------------------------------

table(pancan_mut_expression$nmd.detail2)
#   FALSE    TRUE 
# 1202212   73855 
gene_id = "NF1" # TP53

# TP53 in all cancer types  
tp53 <- pancan_mut_expression[pancan_mut_expression$gene == "NF1",]
tp53$effect[grep("FS",tp53$effect)] <- "Frameshift indels"
tp53$effect[grep("IF",tp53$effect)] <- "Inframe indels"
tp53$cancer <- as.factor(tp53$cancer)
colnames(tp53)[colnames(tp53)=="expression"] <- "Expression"
tp53$Expression <- as.numeric(tp53$Expression)
tp53$effect[tp53$effect%in%c("NM","Frameshift indels")] <-"FS/NM"

wilcox.test(Expression~nmd.detail2,data = tp53, alternative = "greater")


gene_df <- 
  top_tsg <- c("KDM6A","ARID1A","TP53","RB1")
NMDtarget <- pancan_mut[pancan_mut_expression$gene %in% top_tsg & pancan_mut_expression$cancer=="BLCA",] 
ggplot(aes(x = nmd.detail2,fill=(nmd.detail2),y=as.numeric(expression)), data=NMDtarget[(!is.na(NMDtarget$expression))&(!is.na(NMDtarget$nmd.detail2)),] ) + geom_violin(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0,alpha=0.7))+
  labs(x="Cancers", y = "Expression level",fill = "NMD mutations",col =  "NMD mutations",alpha="") +facet_wrap(~ gene,ncol = 1)+theme(legend.position="none")

# p-value < 2.2e-16

#--------------------------------------------------
# Calculation the test (U) statistic per gene-----
#--------------------------------------------------
# rather than using a normal approximation
# calculate the U/z score and MWW p-value for each gene

# two background
# ~calculate for pan cancer------
mut<- pancan_mut_expression
U_pan <- matrix(NA, nrow=nrow(unique(mut[mut$nmd.detail2==T,c("gene","cancer")])),ncol = 10) # calculate the U/z score and MWW p-value for each gene
colnames(U_pan) <- c("gene","cancer","U","Z_score","nmd_mut_expression_mean","bg_expression_mean","p_value","n_nmd_mut","n_bg","bg_normal")
U_pan  <- as.data.frame(U_pan )
U_pan[,c("gene","cancer")] <- unique(mut[mut$nmd.detail2==T,c("gene","cancer_abbr")])
U_pan <- U_pan[(!is.na(U_pan$gene)&(!is.na(U_pan$cancer))),]
U_pan2 <- U_pan # calculate for pan cancer: compare NMD mutation and non-NMD mutation
U_pan3 <- U_pan ## compare non-NMD and WT expression
pancan_mut_expression$expression<-as.numeric(pancan_mut_expression$expression)
for(i in 28648:nrow(U_pan))
{
  gene_id <- U_pan$gene[i]
  cancer <- U_pan$cancer[i]
  m <- mut[mut$gene==gene_id&mut$cancer_abbr==cancer,]
  m <- m[!is.na(m$nmd.detail2),]
  m$expression<-as.numeric(m$expression)
  if(table(m$nmd.detail2)["TRUE"]>0)
  {
    expression_wt <- NULL
    
    expression_wt <- as.numeric(get.wt.expression(gene_id,cancer,nmf = NA))
    
    if(length(expression_wt)>0)
    {
      
      
      results = wilcox.test(as.numeric(m$expression[m$nmd.detail2==T]), expression_wt ,alternative = "less")
      
      U_pan[i,1:6] <- c( gene_id,cancer,results$statistic, results$statistic/(length(m$expression[m$nmd.detail2==T])* length(expression_wt)),
                         mean(as.numeric(m$expression[m$nmd.detail2==T])), mean(expression_wt))
      U_pan[i,7] <-results$p.value
      U_pan[i,8:9] <-c (length(m$expression[m$nmd.detail2==T]), length(expression_wt))
      if(any(expression_wt>0))
      {
        if(length(expression_wt)>2)
        {
          rl <- shapiro.test(expression_wt)
          U_pan[i,10] <- rl$p.value
        }
        
      }
      if(table(m$nmd.detail2)["TRUE"]>0 & table(m$nmd.detail2)["FALSE"]>0)
      {
        results = wilcox.test(as.numeric(m$expression[m$nmd.detail2==F]), expression_wt)
        U_pan3[i,1:6] <- c( gene_id,cancer,results$statistic, results$statistic/(length(m$expression[m$nmd.detail2==F])* length(expression_wt)),
                            mean(as.numeric(m$expression[m$nmd.detail2==F])), mean(expression_wt))
        U_pan3[i,7] <-results$p.value
        U_pan3[i,8:9] <-c (length(m$expression[m$nmd.detail2==F]), length(expression_wt))
        if(any(expression_wt>0))
        {
          if(length(expression_wt)>2){
            rl <- shapiro.test(expression_wt)
            U_pan3[i,10] <- rl$p.value
          }
          
        }
      }
      
      
    }
    
    if(table(m$nmd.detail2)["TRUE"]>0 & table(m$nmd.detail2)["FALSE"]>0)
    {
      
      results = wilcox.test(as.numeric(m$expression[m$nmd.detail2==T]), as.numeric(m$expression[m$nmd.detail2==F]) ,alternative = "less")
      U_pan2[i,1:6] <- c( gene_id,cancer,results$statistic, results$statistic/(length(m$expression[m$nmd.detail2==T])* length(m$expression[m$nmd.detail2==F])),
                          mean(as.numeric(m$expression[m$nmd.detail2==T])), mean(as.numeric(m$expression[m$nmd.detail2==F])))
      U_pan2[i,7] <-results$p.value
      U_pan2[i,8:9] <-c (length(m$expression[m$nmd.detail2==T]), length(m$expression[m$nmd.detail2==F]))
      
      
    }
    
  }
}

for(i in 3:10)
{
  U_pan[,i] <- as.numeric(U_pan[,i])
}
for(i in 3:10)
{
  U_pan2[,i] <- as.numeric(U_pan2[,i])
}
for(i in 3:10)
{
  U_pan3[,i] <- as.numeric(U_pan3[,i])
}

# ~write csv of U-----
write.csv(U_pan,file = "PANCAN/results/U_pan.RData",row.names = F)
write.csv(U_pan2,file = "PANCAN/results/U_pan2.RData",row.names = F)
write.csv(U_pan3,file = "PANCAN/results/U_pan3.RData",row.names = F)

# these genes have more than 1 mutations in 
gene.universe<- U_pan2$Z_score
names(gene.universe) <- U_pan2$gene
gene.universe <- gene.universe[!is.na(gene.universe)]
GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe, geneSel = function(x) x<0.2, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")

table <- GenTable(GOdata, classicFisher = resultFisher,weightfisher = resultweight, topNodes = 500 ,numChar=100)
# write.csv(table,file = "pan_cancer_results/twoBackgroundDiff_GOEA.csv")               
showSigOfNodes(GOdata, score(resultweight), firstSigNodes = 5, useInfo = 'all')

## Which genes are more frequently NMD mutated?------
write.csv(U_pan, file = "panCan12_results/MWW_table_pancancer20160917.csv",row.names = F)


U_pan$cancer_gene <-NA

View(cancergene[cancergene$Gene.Symbol %in% names(table(U_pan$gene))[table(U_pan$gene)>7 & names(table(U_pan$gene)) %in% U_pan$gene[U_pan$p_value < 0.05]] ,])

U_pan$cancer_gene[U_pan$gene %in% cancergene$Gene.Symbol] <- "unclassified"
U_pan$cancer_gene[U_pan$gene %in% cancergene$Gene.Symbol[cancergene$Role.in.Cancer=="oncogene"]] <- "oncogene"
U_pan$cancer_gene[U_pan$gene %in% cancergene$Gene.Symbol[cancergene$Role.in.Cancer=="TSG"]] <- "TSG"
U_pan$cancer_gene[U_pan$gene %in% cancergene$Gene.Symbol[cancergene$Role.in.Cancer=="oncogene/TSG"]] <- "oncogene/TSG"

ggplot(aes(x=n_nmd_mut, y = Z_score, col=as.factor(cancer_gene)),data=U_pan)+geom_point(alpha=0.7)+
  theme(legend.position="bottom")+
  labs(x="Number of NMD mutations",y= "Z score (NMD efficiency)")+
  geom_text(aes(label=as.character(paste(gene,cancer,sep="\n"))),hjust=-0.1, vjust=0,size=1,data = U_pan[(U_pan$n_nmd_mut>20&U_pan$Z_score<0.3)|(U_pan$n_nmd_mut>40),])

### gene features------

gene_list
gene_list <- U_pan$gene
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
hg19_txlens <- transcriptLengths(txdb, with.cds_len=TRUE,
                                 with.utr5_len=TRUE,
                                 with.utr3_len=TRUE)
hg19_txlens <- hg19_txlens[which(hg19_txlens$gene_id %in% gene_list),]
pancan_mut_expression$pos.ptc <- as.numeric(pancan_mut_expression$pos.ptc)
ggplot(data=pancan_mut[(pancan_mut_expression$nmd.detail2==T),],aes(x=paste(wt_expression_rank,as.numeric(pos.ptc)<200),y=quantile_expression,fill=as.numeric(pos.ptc)<200))+
  geom_boxplot(alpha=0.3,position = "dodge")

## How do they differ?-----

## Plot NMD affected genes-----

# ~ Supplementary Figure 9. (A) . Expression levels of 25 genes with z scores (NMD vs wildtype) < 0.1. -----
gene_id = paste(U_pan$gene,U_pan$cancer)[U_pan$n_nmd_mut>2&U_pan$Z_score<0.1]
gene_id <- gene_id[!is.na(gene_id )]
NMDtarget <- pancan_mut_expression[paste(pancan_mut_expression$gene,pancan_mut_expression$cancer_abbr) %in% gene_id[1:25],] 
NMDtarget$z_score <-round (U_pan$Z_score[match(paste(NMDtarget$gene,NMDtarget$cancer_abbr),paste(U_pan$gene,U_pan$cancer))],2)

NMDtarget$gene <- paste(NMDtarget$gene,NMDtarget$cancer_abbr,"z score=",NMDtarget$z_score)
NMDtarget <- NMDtarget[!is.na(NMDtarget$z_score)&!is.na(NMDtarget$nmd.detail2),]
ggplot(aes(x = nmd.detail2,fill=nmd.detail2,y=expression), data=NMDtarget ) + geom_boxplot(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0,alpha=0.7))+
  labs(x="Cancers", y = "Expression level",fill = "NMD mutations",col =  "NMD mutations",alpha="") +facet_wrap(~ gene,ncol = 5)+
  scale_alpha(labels=NULL)+theme(legend.position="bottom")
# scale_x_discrete(breaks=unique(paste(pancan_mut_expression$cancer,pancan_mut_expression$nmd.detail2)),
#                  labels=rep(c(cancer_type$Cancer),each=2))

# ~ Supplementary Figure 9. (B) . Expression levels of 25 genes with z scores (NMD vs non-NMD) < 0.1. -----
gene_id = paste(U_pan2$gene,U_pan2$cancer)[U_pan2$n_nmd_mut>2&U_pan2$Z_score<0.1]
gene_id <- gene_id[!is.na(gene_id )]
NMDtarget <- pancan_mut_expression[paste(pancan_mut_expression$gene,pancan_mut_expression$cancer_abbr) %in% gene_id[1:25],] 
NMDtarget$z_score <-round (U_pan2$Z_score[match(paste(NMDtarget$gene,NMDtarget$cancer_abbr),paste(U_pan2$gene,U_pan2$cancer))],2)

# NMDtarget$z_score <-round (U_pan$Z_score[match(paste(NMDtarget$gene,NMDtarget$cancer),paste(U_pan$gene,U_pan$cancer))],2)
NMDtarget$gene <- paste(NMDtarget$gene,NMDtarget$cancer_abbr,"z score=",NMDtarget$z_score)
NMDtarget <- NMDtarget[!is.na(NMDtarget$z_score)&!is.na(NMDtarget$nmd.detail2),]

# NMDtarget$gene <- paste(NMDtarget$gene,NMDtarget$cancer,"z score=",NMDtarget$z_score)
# NMDtarget <- NMDtarget[!is.na(NMDtarget$z_score),]
ggplot(aes(x = nmd.detail2,fill=nmd.detail2,y=expression), data=NMDtarget ) + geom_boxplot(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0,alpha=0.7))+
  labs(x="Cancers", y = "Expression level",fill = "NMD mutations",col =  "NMD mutations",alpha="") +facet_wrap(~ gene,ncol = 5)+
  scale_alpha(labels=NULL)+theme(legend.position="bottom")
# scale_x_discrete(breaks=unique(paste(pancan_mut_expression$cancer,pancan_mut_expression$nmd.detail2)),
# labels=rep(unique(pancan_mut_expression$cancer),each=2))

# Tumour-suppressor effect by NMD------------
#--------------------------------------------
#--------------------------------------------
cancergene<- read.csv("pan_cancer_data/Census_allWed Jul 13 11_38_59 2016.csv",header=T,as.is = T)
U_pan$cancer_gene[is.na(U_pan$cancer_gene)] <- "unrelated"
# z score distribution of COSMIC
ggplot(U_pan[U_pan$cancer_gene %in% c("TSG","oncogene","unrelated"),], aes(x = Z_score,col = cancer_gene,fill=cancer_gene)) +
  geom_histogram(aes(y=..density..),binwidth = 0.02,position = "dodge")+ labs (x = "z score (NMD vs wildtype)", fill = "Cancer-related genes",col = "Cancer-related genes")+
  theme(legend.position="bottom")+facet_wrap(~cancer_gene,ncol=1)

t.test(U_pan$Z_score[U_pan$cancer_gene=="TSG"],U_pan$Z_score[U_pan$cancer_gene=="unrelated"],alternative = "less")
t.test(U_pan$Z_score[U_pan$cancer_gene=="oncogene"],U_pan$Z_score[U_pan$cancer_gene=="unrelated"],alternative = "less")

sum(U_pan$n_nmd_mut[which(U_pan$cancer=="unclassified")])

table(cancergene$Role.in.Cancer)
table(U_pan$cancer)/sum(table(U_pan$cancer))
c(sum(U_pan$n_nmd_mut[which(U_pan$cancer=="unclassified")]),
  sum(U_pan$n_nmd_mut[which(U_pan$cancer=="oncogene")]),
  sum(U_pan$n_nmd_mut[which(U_pan$cancer=="oncogene/TSG")]),
  sum(U_pan$n_nmd_mut[which(U_pan$cancer=="TSG")]))/sum(c(sum(U_pan$n_nmd_mut[which(U_pan$cancer=="unclassified")]),
                                                          sum(U_pan$n_nmd_mut[which(U_pan$cancer=="oncogene")]),
                                                          sum(U_pan$n_nmd_mut[which(U_pan$cancer=="oncogene/TSG")]),
                                                          sum(U_pan$n_nmd_mut[which(U_pan$cancer=="TSG")])))
fisher.test(matrix(c(84,55,32,40),byrow = T,ncol = 2),alternative = "greater")

ggplot(U_pan2, aes(x = Z_score,fill  = smg, col  = smg)) +
  geom_density(alpha=0.5) + labs (x = "z score (NMD vs non-NMD)", fill = "Significantly mutated genes",col = "Significantly mutated genes")+
  theme(legend.position="bottom")


## some example 
pancan_mut[pancan_mut_expression$gene==1029,]
ggplot(aes(x=cancer,fill=nmd.detail2,y=expression),data=pancan_mut[pancan_mut_expression$gene==1029,])+
  geom_boxplot() + geom_abline(intercept = U_pan$wt_expression_mean[U_pan$gene_id==1029],slope=0)

# another gene list from nature 2013 pan-cancer analysis, significantly mutated genes------
smg <- read.csv("pan_cancer_data/SMG_cancer_Kandoth_2013.csv")
sng.list <- smg$Gene
sum(sng.list %in% pancan_mut_expression$gene) #119
length(sng.list)
sng.list.geneId <- pancan_mut_expression$gene[match(sng.list,pancan_mut_expression$gene)]

ggplot(aes(x=n_nmd_mut, y = Z_score),data=U_pan[U_pan$gene %in% sng.list.geneId,])+geom_point(alpha = 5/10)+
  theme(legend.position="bottom")+
  labs(x="Number of NMD mutations",y= "Z score (NMD vs wildtype)",col="Role in cancer")+
  geom_text(aes(label=as.character(gene)),hjust=0, vjust=0,data = U_pan[U_pan$gene %in% sng.list.geneId ,])

U_pan$smg <- F
U_pan$smg[U_pan$gene %in% sng.list.geneId ] <- T
t.test(U_pan$Z_score[U_pan$smg==T],U_pan$Z_score[U_pan$smg==F],alternative = "less")
library(plyr)
cdat <- ddply(U_pan, "smg", summarise, rating.mean=median(Z_score,na.rm = T))
cdat

# plot: z score distribution of SMGs------
ggplot(U_pan, aes(x = Z_score,fill  = smg, col  = smg)) +
  geom_density(alpha=0.3) + labs (x = "z score (NMD vs wildtype)", fill = "Significantly mutated genes",col = "Significantly mutated genes")+
  theme(legend.position="bottom")+
  geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=smg),
             linetype="dashed", size=1)

gene_id = sng.list.geneId[sng.list.geneId %in%  U_pan2$gene_id]
NMDtarget <- pancan_mut[pancan_mut_expression$gene %in% gene_id[1:25],] 
NMDtarget$z_score <-round (U_pan2$Z_score[match(NMDtarget$gene,U_pan$gene_id)],2)
NMDtarget$gene <- paste(NMDtarget$gene,NMDtarget$Hugo_Symbol,"z score=",NMDtarget$z_score)

ggplot(aes(x = paste(cancer,nmd.detail2),fill=nmd.detail2,y=expression), data=NMDtarget ) + geom_boxplot(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0),alpha=0.5)+ 
  scale_colour_discrete(guide = FALSE) +scale_fill_discrete(guide = FALSE)+
  labs(x="Cancers", y = "Expression level") +
  scale_x_discrete(breaks=unique(paste(pancan_mut_expression$cancer,pancan_mut_expression$nmd.detail2)),
                   labels=rep(c("BRCA","GBM","LUSC","OV","PRAD"),each=2))+facet_wrap(~ gene,ncol = 5) 

gene_id = sng.list.geneId[sng.list.geneId %in%  U_pan2$gene_id]
NMDtarget <- pancan_mut[pancan_mut_expression$gene %in% gene_id[26:51],] 
NMDtarget$z_score <-round (U_pan2$Z_score[match(NMDtarget$gene,U_pan$gene_id)],2)
NMDtarget$gene <- paste(NMDtarget$gene,NMDtarget$Hugo_Symbol,"z score=",NMDtarget$z_score)

ggplot(aes(x = paste(cancer,nmd.detail2),fill=nmd.detail2,y=expression), data=NMDtarget ) + geom_boxplot(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0),alpha=0.5)+ 
  scale_colour_discrete(guide = FALSE) +scale_fill_discrete(guide = FALSE)+
  labs(x="Cancers", y = "Expression level") +
  scale_x_discrete(breaks=unique(paste(pancan_mut_expression$cancer,pancan_mut_expression$nmd.detail2)),
                   labels=rep(c("BRCA","GBM","LUSC","OV","PRAD"),each=2))+facet_wrap(~ gene,ncol = 5) 

gene_id = sng.list.geneId[sng.list.geneId %in%  U_pan2$gene_id]
NMDtarget <- pancan_mut[pancan_mut_expression$gene %in% gene_id[52:84],] 
NMDtarget$z_score <-round (U_pan2$Z_score[match(NMDtarget$gene,U_pan$gene_id)],2)
NMDtarget$gene <- paste(NMDtarget$gene,NMDtarget$Hugo_Symbol,"z score=",NMDtarget$z_score)

ggplot(aes(x = paste(cancer,nmd.detail2),fill=nmd.detail2,y=expression), data=NMDtarget ) + geom_boxplot(alpha=0.4) + 
  geom_jitter(aes(col = nmd.detail2,height = 0),alpha=0.5)+ 
  theme(legend.position = "bottom")+
  labs(x="Cancers", y = "Expression level",fill="NMD mutation",col="NMD mutation") +
  scale_x_discrete(breaks=unique(paste(pancan_mut_expression$cancer,pancan_mut_expression$nmd.detail2)),
                   labels=rep(c("BRCA","GBM","LUSC","OV","PRAD"),each=2))+facet_wrap(~ gene,ncol = 5) 


# detect possible wrong mutations------
#  get flanking sequencing-----
library('BSgenome.Hsapiens.UCSC.hg19')

cutoff <- 3  # at least 3 repeat in length
getflan
chr <- 'chr10'
position <- 890939
alleles <- '[T/-]'
offset <- 10

# Function to pull flanking sequence. Defaults to +/- 10 bp
# adapted from http://www.gettinggeneticsdone.com/2011/04/using-rstats-bioconductor-to-get.html [20/11/2016]
getflank <- function(position, alleles="N", chr="chr12", offset=3) {
    leftflank  <- getSeq(Hsapiens,chr,position-offset,position-1)
    rightflank <- getSeq(Hsapiens,chr,position+1,position+offset)
    paste(leftflank,alleles,rightflank,sep="")
}
indel_list <- pancan_mut_all[pancan_mut_all$effect%in%c("FSD","FSI","IFD","IFI"),]
dim(indel_list) #76041  4.3%
table(indel_list$nmd.detail2)
sum(nchar(indel_list$reference)==1&nchar(indel_list$alt)==1)
# [1] 53778
table(pancan_mut_all$nmd.detail2)
# FALSE         TRUE unclassified 
# 1429163        88037       261610 
table(indel_list$nmd.detail2)
# FALSE         TRUE unclassified 
# 30383        37346         8312 
indel_list$flank_seq <- NA
indel_list$sn_repeat <- NA # is it a single nucleotide repeat
for(i in 1:nrow(indel_list))
{
    if(nchar(indel_list$reference[i])==1&nchar(indel_list$alt[i])==1&indel_list$chr[i]!="chrMT")
    {
        allele <- c(indel_list$reference[i],indel_list$alt[i])
        allele <- allele[allele!="-"]
        indel_list$flank_seq[i] <- getflank(position = indel_list$start[i], # get the flanking sequence
                                            allele=allele,
                                            chr = indel_list$chr[i])
        reads  <- c(substr(indel_list$flank_seq[i],start = 1,stop = 4),
                    substr(indel_list$flank_seq[i],start = 2,stop = 5),
                    substr(indel_list$flank_seq[i],start = 3,stop = 6),
                    substr(indel_list$flank_seq[i],start = 4,stop = 7))
        indel_list$sn_repeat[i] <- any(reads %in% c("AAAA","TTTT","CCCC","GGGG"))
    }
    
}

write.table(indel_list,"PANCAN_results/20161120repeats_indels.txt",sep="\t",
            row.names = F)
# indel_list <- read.delim(file = "PANCAN_results/20161120repeats_indels.txt")

ggplot(data.frame(table(indel_list[,c("MSS","sn_repeat","nmd.detail2","cancer_abbr")]))[data.frame(table(indel_list[,c("MSS","sn_repeat","nmd.detail2","cancer_abbr")]))$Freq>0,],
       aes(x=paste(cancer_abbr,MSS),y=Freq,fill=paste(nmd.detail2,sn_repeat)))+
    geom_bar(stat = "identity",position = "stack")+coord_flip()+
    scale_fill_discrete(name="Mutations/flanking nucleotides",
                        breaks=c("FALSE FALSE", "FALSE TRUE", "TRUE FALSE",
                                 "TRUE TRUE","unclassified FALSE","unclassified TRUE"),
                        labels=c("nonNMD-elicit unrepeated", 
                                 "nonNMD-elicit repeated",
                                 "NMD-elicit unrepeated", 
                                 "NMD-elicit repeated",
                                 "unclassified unrepeated",
                                 "unclassified repeated"))+
    labs(y="counts",x="cancer/MS instability")

table(indel_list$MSS[indel_list$sn_repeat==T])

pancan_mut_all <- cbind(pancan_mut_all,
                        indel_list[match(paste(pancan_mut_all$sample,pancan_mut_all$gene,pancan_mut_all$start,pancan_mut_all$alt),
                                         paste(indel_list$sample,indel_list$gene,indel_list$start,indel_list$alt)),
                                   c("flank_seq","sn_repeat")])

pancan_mut_expression <- cbind(pancan_mut_expression,
                        indel_list[match(paste(pancan_mut_expression$sample,pancan_mut_expression$gene,pancan_mut_expression$start,pancan_mut_expression$alt),
                                         paste(indel_list$sample,indel_list$gene,indel_list$start,indel_list$alt)),
                                   c("flank_seq","sn_repeat")])


repeat_mut_loci <- unique(indel_list[,c("chr","start","gene")])
repeat_mut_loci$chr <- as.character(repeat_mut_loci$chr)
nrow(repeat_mut_loci)
# [1] 58295
for(i in 1:nrow(repeat_mut_loci))
{
    
    repeat_mut_loci$freq[i] <- sum(indel_list$chr==repeat_mut_loci$chr[i]&indel_list$start==repeat_mut_loci$start[i],na.rm = T)
}

max(repeat_mut_loci$freq)

repeat_mut_loci$flank_seq <- indel_list$flank_seq[!is.na(indel_list$flank_seq)][match(paste(repeat_mut_loci$chr,repeat_mut_loci$start),
                                                        paste(indel_list$chr,indel_list$start)[!is.na(indel_list$flank_seq)])]
repeat_mut_loci$sn_repeat <- indel_list$sn_repeat[!is.na(indel_list$flank_seq)][match(paste(repeat_mut_loci$chr,repeat_mut_loci$start),
                                                                                      paste(indel_list$chr,indel_list$start)[!is.na(indel_list$flank_seq)])]
repeat_mut_loci$nmd <- indel_list$nmd.detail2[!is.na(indel_list$flank_seq)][match(paste(repeat_mut_loci$chr,repeat_mut_loci$start),
                                                                                      paste(indel_list$chr,indel_list$start)[!is.na(indel_list$flank_seq)])]

write.table(repeat_mut_loci,"PANCAN_results/20161120repeat_mut_loci.txt",sep="\t",
            row.names = F)

# Hypermutation---------

# PTEN----

# NMD-elicit mutations in PTEN were associated with a profound and significant 
# reduction in PTEN expression (P = 5.6e-09, MWW test). 

pten_a <- pancan_mut_all$expression[pancan_mut_all$cancer_abbr=="STAD"&pancan_mut_all$gene=="PTEN"&pancan_mut_all$nmd.detail2==T]
pten_b <- pancan_mut_all$expression[pancan_mut_all$cancer_abbr=="STAD"&pancan_mut_all$gene=="PTEN"&pancan_mut_all$nmd.detail2==F]
