# load packages-----
library(gplots)
library(scales)
library(ggplot2)
library(ggrepel)
library(Rmisc)
library(plyr)
library(org.Hs.eg.db)
library(RColorBrewer)
require(grid)
library(gridExtra)
require(cowplot)
# Create a transparent theme object
transparent_theme <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))

get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
# A helper function to define a region on the layout------
define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
} 
# read files #---------------------------------
#--------------------------------------------------

cancer_type <- read.csv("data/cancer_list.csv",header = T,as.is = T)

cancergene<- read.csv("data/Census_allWed Jul 13 11_38_59 2016.csv",header=T,as.is = T)
tsg <- cancergene$Gene.Symbol[grep("TSG",cancergene$Role.in.Cancer)]
pancan_mut_all$driver <- pancan_mut_all$gene %in% cancergene$Gene.Symbol
pancan_mut_all$role_in_cancer <- cancergene$Role.in.Cancer[match(pancan_mut_all$gene,cancergene$Gene.Symbol)]

# plot figure 1 c d-------
# Figure1c----

effect <- c("Silent","NM","MM","IFI","IFD","FSI","FSD")
pancan_mut_all$nmd.detail.new <- NA
pancan_mut_all$nmd.detail.new[pancan_mut_all$nmd.detail2==T] <- "NMD-elicit"
pancan_mut_all$nmd.detail.new[pancan_mut_all$nmd.detail2==F&pancan_mut_all$have.ptc==T] <- "NMD-escape"
pancan_mut_all$nmd.detail.new[pancan_mut_all$have.ptc==F] <- "non-PTC"
table(pancan_mut_all$nmd.detail.new)
##* all genes----
df2<-data.frame(table(pancan_mut_all[(pancan_mut_all$effect %in% effect),c("effect","nmd.detail.new")]))
# Set columns ready for ggplot
df2$Freq<-as.numeric(df2$Freq)
df2$effect <- as.character(df2$effect) #%>% as.numeric()
df2$nmd.detail.new<- as.character(df2$nmd.detail.new)# %>% as.numeric()

# Figure 1c----
ggplot(df2,aes(x=nmd.detail.new,y=effect)) +geom_tile(aes(fill=log10(as.numeric(Freq))))+
  geom_text(aes(label=Freq))+xlab(label = "classification")+scale_x_discrete(labels = c("NMD-elicit","NMD-escape","non-PTC"))+
  ylab(label = "TCGA variant classification")+theme(legend.position="none")+ 
    scale_fill_gradient2(high = "#c6e5fc", low = "#ecf6fe", space = "Lab", guide = "colourbar")+
  theme_bw()+theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line =  element_blank(),panel.border = element_blank(),
                   legend.position="none")


# Figure 1d-----
ggplot(data=sample_df,aes(x=as.numeric(mut),y=as.numeric(NMD_mut)))+geom_point(col="grey50",alpha=0.2,size=0.5)+
  scale_x_log10() +
  scale_y_log10()+geom_abline(slope = 1,intercept = 0,col="grey",lty=2)+
    facet_wrap(~cancer,ncol = 6)+
    labs(x="Number of mutations in one sample",y="Number of NMD-elicit mutations in one sample")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+theme_bw()+  theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.key = element_rect( color = "white")) +theme(legend.position="none")


# figure S1------
ggplot(aes(x=n_nmd_mut, y = Z_score, col=as.factor(cancer_gene)),data=U_pan[U_pan$gene!="TTN",])+geom_point(alpha=0.5)+
  theme(legend.position="bottom")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  labs(x="Number of NMD-elicit mutations",y= "Z score",colour="")+scale_color_discrete(labels=c("oncogene","oncogene/TSG","TSG","other cancer-related","cancer-unrelated"))+
  geom_text_repel(aes(label=as.character(paste(gene,cancer,sep="\n"))),size=2,data = U_pan[((U_pan$n_nmd_mut>30&U_pan$Z_score<0.3)|(U_pan$n_nmd_mut>40))&U_pan$gene!="TTN",])+
  theme_bw()+  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     legend.key = element_rect( color = "white")) +
  theme(legend.position = "bottom",legend.direction="horizontal",legend.title = element_blank()) 

# figure 3c------

# table of mutation
gene_list_signitaure <- unique(c("TP53","MEN1",
                          "KDM6A","ARID1A","RB1","STAG2","EP300","FAT1",
                          "MAP3K1","CDH1",
                          "APC","RNF43","ATM","ARID1A",
                          "TNFAIP3",
                          "NF1",
                          "NSD1","CDKN2A","NOTCH1",
                          "PNRM1","VHL","FANCE",
                          "ATRX",
                          "STK11",
                          "SMAD4",
                          "FBXW7",
                          "RB1",
                          "PIK3R1","CTCF","ZFHX3"))
nmd_signature <- table(pancan_mut_all[pancan_mut_all$gene %in% gene_list_signitaure& pancan_mut_all$new.annotation!="Silent",c("gene","nmd.detail2","sample")])
nmd_signature <- data.frame(nmd_signature ,stringsAsFactors = F)
nmd_signature$Freq[nmd_signature$Freq>0] <-1
nmd_signature <- nmd_signature[nmd_signature$Freq>0,]
nmd_signature$gene <- factor(as.character(nmd_signature$gene),levels = gene_list_signitaure[28:1])

nmd_signature2_cancer_panel <- data.frame(sample = unique(nmd_signature$sample),
                                          cancer = sample_df$cancer[match(unique(nmd_signature$sample),sample_df$sample)])
nmd_signature2_cancer_panel$n_mut  = sample_df$mut[match(unique(nmd_signature$sample),sample_df$sample)]
sample_level = nmd_signature2_cancer_panel$sample[order(paste(nmd_signature2_cancer_panel$cancer,
                                                              (100000000-nmd_signature2_cancer_panel$n_mut)),
                                                        decreasing = F)]
nmd_signature2_cancer_panel$sample <- factor(as.character(nmd_signature2_cancer_panel$sample),
                                             levels=sample_level)
nmd_signature2_cancer_panel$cancer <- factor(as.character(nmd_signature2_cancer_panel$cancer ),
                                             levels = c('ACC','BLCA','BRCA','CESC','COAD','DLBC','ESCA','GBM','HNSC','KIRC','KIRP','LGG','LIHC','LUAD','LUSC','OV','PAAD','PRAD','READ','SARC','SKCM','STAD','THCA','UCEC') )
# figure 2c-bar----
ggplot(nmd_signature2_cancer_panel,aes(x=sample,y=1,fill=cancer))+
    geom_tile()+
    scale_fill_manual(values=c(brewer.pal(name = "Set3",12),brewer.pal(name = "Set2",8),brewer.pal(name = "Set1",4)))+
    transparent_theme+theme(legend.position="bottom")+labs(fill="")+
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 1,ncol = 12, 
                               byrow = TRUE,
                               label.theme = element_text(size = 7,angle=0)))



nmd_signature$sample <- factor(as.character(nmd_signature$sample),levels=sample_level)
nmd_signature$nmd.detail2 <- as.character(nmd_signature$nmd.detail2)
nmd_signature$nmd.detail2[nmd_signature$nmd.detail2==T] <-"NMD-elicit mutation"
nmd_signature$nmd.detail2[nmd_signature$nmd.detail2==F] <-"other nonsilent mutations"

# Figure2c-heatmap----
ggplot(data = nmd_signature ,aes(x=sample,y=gene,fill=nmd.detail2))+geom_tile()+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
          axis.line=element_blank(),
          legend.position="top")+labs(x="",y="",fill="")+
   scale_fill_manual(values=c("#FF7F00",alpha("#4DAF4A",alpha = 0.7),alpha("#999999",alpha = 0.7)))

# figure 2a TP53----------
labels <- c("NMD","MM","non-NMD","MM")

tp53 <- pancan_mut_expression[pancan_mut_expression$gene == "TP53",]
tp53$effect[grep("FS",tp53$effect)] <- "Frameshift indels"
tp53$effect[grep("IF",tp53$effect)] <- "Inframe indels"
tp53$cancer <- as.factor(tp53$cancer)
colnames(tp53)[colnames(tp53)=="expression"] <- "Expression"
tp53$Expression <- as.numeric(tp53$Expression)
tp53$effect[tp53$effect%in%c("NM","Frameshift indels")] <-"FS/NM"

# Figure 2a-1 TP53------
ggplot(aes(x = paste(effect,(nmd.detail2)),fill=nmd.detail2,y=Expression,col=nmd.detail2), 
       data=tp53[!is.na(tp53$nmd.detail2)&(tp53$effect %in% c("FS/NM","MM")),]) + geom_boxplot(alpha=0.4)+ 
    labs(x="", y = "TP53 expression level",fill = "",col="") +
    scale_x_discrete(breaks=c("FS/NM TRUE","FS/NM FALSE","MM FALSE"),
                     labels=labels)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90, vjust=1))+theme(axis.ticks.length = unit(0.1, "cm"))+
    theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
          legend.position="none")


# figure 2b NF1----

labels <- c("NMD-elicit","non-NMD","missense")
tp53 <- pancan_mut_expression[pancan_mut_expression$gene == "NF1",]
tp53$effect[grep("FS",tp53$effect)] <- "Frameshift indels"
tp53$effect[grep("IF",tp53$effect)] <- "Inframe indels"
tp53$cancer <- as.factor(tp53$cancer)
colnames(tp53)[colnames(tp53)=="expression"] <- "Expression"
tp53$Expression <- as.numeric(tp53$Expression)
tp53$effect[tp53$effect%in%c("NM","Frameshift indels")] <-"FS/NM"

# Figure 2a-1 NF1----
ggplot(aes(x = paste(effect,(nmd.detail2)),fill=nmd.detail2,y=Expression,col=nmd.detail2), 
       data=tp53[!is.na(tp53$nmd.detail2)&(tp53$effect %in% c("FS/NM","MM")),]) + geom_boxplot(alpha=0.4)+ 
    labs(x="", y = "NF1 expression level",fill = "",col="") +
    scale_x_discrete(breaks=c("FS/NM TRUE","FS/NM FALSE","MM FALSE"),
                     labels=labels)+theme_bw()+  
    theme(legend.position="none")+
    theme(axis.text.x = element_text( angle=90, vjust=1))+
    theme(axis.ticks.length = unit(0.1, "cm"))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_rect( color = "white"))


# figure 2New: compare MM NMD-elicit and NMD-escape------


theme.ggplot <- theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      legend.key = element_rect( color = "white")) +theme(legend.position="bottom")
pancan_mut_expression$new.annotation2 <- pancan_mut_expression$new.annotation
pancan_mut_expression$new.annotation2[grep("non-PTC",pancan_mut_expression$new.annotation)] <- "non-PTC"

# figure 1d-----

ggplot(aes(x = new.annotation2,fill=nmd.detail2,y=quantile_expression), 
       data=pancan_mut_expression[as.numeric(pancan_mut_expression$wt_expression_rank)>4&
                                      pancan_mut_expression$nmd.detail2%in%c(T,F),]) + 
    geom_boxplot()+ geom_hline(yintercept = 0.5,lty=2)+
    labs(x="Reannotation of mutations", y = "Relative expression of variants",fill="")+
    scale_fill_discrete(labels=c("non-NMD-elicit","NMD-elicit")) +
    theme(legend.position="bottom",axis.text.x=element_text(angle=30,vjust=0.6))


t.test(x=pancan_mut_expression[as.numeric(pancan_mut_expression$wt_expression_rank)>4&pancan_mut_expression$nmd.detail2%in%c(T,F)&pancan_mut_expression$new.annotation2=="NMD-escape","quantile_expression"],
       y=pancan_mut_expression[as.numeric(pancan_mut_expression$wt_expression_rank)>4&pancan_mut_expression$nmd.detail2%in%c(T,F)&pancan_mut_expression$new.annotation2=="non-PTC","quantile_expression"],alternative = "less")


## Figure 4a-----
ggplot(sample_df[!is.na(sample_df$hyper_mutated),],
       aes(x=as.numeric(mut),y=as.numeric(NMD_mut),col=factor(TSG_NMD_mut),shape=factor(hyper_mutated)))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_point(alpha=0.3)+
  facet_wrap(~cancer,ncol=3)+scale_x_log10()+scale_y_log10()+geom_abline(slope = 1,intercept = 0,col="grey",lty=2)+
  labs(x="# mutations in one sample",y="# NMD mutations in one sample",col="# NMD-mutated TSGs",shape="Hypermutated sample")+
  theme_bw()+  theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.key = element_rect( color = "white")) +
  theme(legend.position="bottom")

# figure 4b------
ggplot(sample_df[!is.na(sample_df$hyper_mutated),],
       aes(x=factor(cancer), y=as.numeric(TSG_NMD_mut),fill=factor(hyper_mutated)))+
  geom_boxplot()+
  geom_hline(yintercept = -1)+ylim(-1,20)+
  labs(x="",y="# NMD-mutated TSGs per sample",fill="Hypermutated samples")+
  theme_bw()+  
  theme(panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(),
                                                panel.border = element_blank(),
                                                panel.background = element_blank(),
                                                legend.key = element_rect( color = "white")) +
  theme(legend.position="bottom")

t.test(x=sample_df$TSG_NMD_mut[sample_df$hyper_mutated==T],y=sample_df$TSG_NMD_mut[sample_df$hyper_mutated==F],alternative = "greater")

# figure 3a-----
ggplot(go_mut_tb3,aes(x=sample,y=GO,fill=NMD_mutated))+geom_tile(col="white")+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=8),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    theme(legend.position="bottom")+
    labs(title = "")+scale_fill_manual(values=color_heatmap)+
    guides(fill = guide_legend(nrow=3,title ="",byrow = F))


ggplot(go_mut_sample,aes(x=sample,y=Freq,fill=new.annotation))+
    geom_bar(stat="identity",position = "stack")+labs(x="",fill="")+
    theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),legend.position="top")
dev.off()

# figure 3b-----
ggplot(go_mut_tb3[go_mut_tb3$sample %in% sample_df$sample[sample_df$cancer=="STAD"],],aes(x=sample,y=GO,fill=NMD_mutated))+geom_tile(col="white")+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(size=8),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    theme(legend.position="bottom")+
    labs(title = "")+scale_fill_manual(values=color_heatmap)+
    guides(fill = guide_legend(nrow=5,title ="",byrow = F))

