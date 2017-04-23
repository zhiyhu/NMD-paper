source("00_PANCAN_function.R")

# Classify mutations into NMD and non-NMD

# Load package----------------
#-------------------------------
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# biocLite('BSgenome.Hsapiens.UCSC.hg19')
# biocLite("GenomicFeatures")

library(TxDb.Hsapiens.UCSC.hg19.knownGene);library('BSgenome.Hsapiens.UCSC.hg19') # NCBI build 37
library("Biostrings")
library("GenomicFeatures")
library(stringr)
library(scales)
txdb37 <- TxDb.Hsapiens.UCSC.hg19.knownGene# NCBI build 37
genome37 <- BSgenome.Hsapiens.UCSC.hg19# NCBI build 37
library("biomaRt") # ensembl genome annotation

# prepare pos info of transcript-----
#------------------------------------
gene_list <- unique(pancan_mut_all$geneId)

tx37 <- select(txdb37, keys = as.character(unique(gene_list)),
               columns=c("GENEID","TXID","TXSTART","TXEND","EXONID","EXONCHROM",
                         "EXONSTRAND","EXONSTART","EXONEND","CDSEND" ,"CDSID","CDSSTART"), keytype="GENEID")
tx37 <- tx37[!is.na(tx37$CDSSTART),]

gene2tx37 <- select(txdb37, keys = as.character(unique(gene_list)), columns="TXID", keytype="GENEID")
cds <- cdsBy(txdb37, by="tx" )
cds_seqs37 <- extractTranscriptSeqs(genome37, cds) # all CDS in NCBI build 37
cds_seqs37

# database of build 37, have all CDS
cds37 <- data.frame(txid = names(cds_seqs37),geneid = gene2tx37$GENEID[match(names(cds_seqs37),gene2tx37$TXID)],seq = as.character(cds_seqs37),width = width(cds_seqs37))
cds37 <- cds37[!is.na(cds37$geneid),]
cds37 <- cds37[!duplicated(paste(cds37$geneid,cds37$seq)),]

tx37 <- tx37[tx37$TXID %in% cds37$txid,]

# Classify the mutations into NMD and non-NMD----------
#------------------------------------------------------
nrow(pancan_mut)

# exclude mutations outside coding sequence
pancan_mut<- pancan_mut[!(pancan_mut$effect %in% c("Splice_Site","3'UTR","5'UTR","IGR","Intron","RNA","targeted_region")),] 
nrow(pancan_mut)
# [1] 1344589

# filter the mutations without annotations in cds db
has.anno <- which(pancan_mut$geneId %in% cds37$geneid)
pancan_mut <- pancan_mut[has.anno,]; rm(has.anno)
nrow(pancan_mut) 
# [1] 1343298
save(pancan_mut,file="data/results/pancan_mut_before_assess")
# see if the gene has correct CDS annotation: start with start codon, end to stop codon
pancan_mut$build  <- 37

# find dominate transcript for each gene in each cancer
assess <- data.frame(gene=as.numeric(unique(pancan_mut$geneId))) # extract all genes

assess$assess <- NA # 17280
for(i in 1:nrow(assess)) # see if the gene has correct CDS annotation: start with start codon, end to stop codon
  assess$assess[i] <-  assess.annotation(assess$gene[i])

write(unique(pancan_mut$gene[!pancan_mut$geneId %in% assess$gene[assess$assess==0]]),
      file="data/results/filtered_gene_list_incorrect_CDS.txt")
pancan_mut <- pancan_mut[pancan_mut$geneId %in% assess$gene[assess$assess==0],] 
# remove the genes are not correct annotated
nrow(pancan_mut)
# [1]  1338241

save(pancan_mut,file="data/results/pancan_mut_for_classify")
write.csv(assess,"data/results/gene_assess_tb.csv",row.names = F)
rm(assess)

# Classify the mutated sequences into NMD T or F----------
nmd.rl <- matrix(ncol=8, nrow = nrow(pancan_mut))
colnames(nmd.rl) <- names(classify.nmd(pancan_mut[1,]))
cancer_type$Cancer[order(cancer_type$Number_total_mutations,decreasing = T)]

sample <- which(pancan_mut$cancer %in% 
                    cancer_type$Cancer[order(cancer_type$Number_total_mutations,decreasing = T)][4:24])
for(i in sample)
{
  if( (pancan_mut$geneId[i] %in% tx37$GENEID & pancan_mut$build[i]==37))
    nmd.rl[i,] <- classify.nmd(pancan_mut[i,])
  else nmd.rl[,"note"] <- "can't find geneID in txdb"
}
save(nmd.rl,file="data/data/nmd_classified")

pancan_mut <- cbind(pancan_mut,nmd.rl[,-3])
pancan_mut <- pancan_mut[pancan_mut$cancer %in% cancer_type$Cancer[order(cancer_type$Number_total_mutations,decreasing = T)][1:24],]
save(pancan_mut,file="data/data/pancan_mut_classified")

colnames(pancan_mut)[which(colnames(pancan_mut)=="mut_nmd")] <- "nmd.detail2"
pancan_mut$nmd.detail2[grep("splice site mutated",pancan_mut$note)] <- NA
sum(is.na(pancan_mut$nmd.detail2))
# [1] 28924
pancan_mut$nmd.detail<- NA
pancan_mut$nmd.detail[which(pancan_mut$nmd.detail2==T)] <- "NMD"
pancan_mut$nmd.detail[pancan_mut$nmd.detail2==F] <- "nonNMD"
pancan_mut$nmd.detail[(grepl("cannot find ORF",pancan_mut$note) | 
                           is.na(pancan_mut$PTC.Stop)) &
                          !is.na(pancan_mut$nmd.detail2)] <- "nonstop"

# Mergy the mutation and RNAseq data --------------------
#-------------------------------------------------------

pancan_seq <- pancan_seq[, colnames(pancan_seq) %in% pancan_mut$sample]
# find the mRNA level of gene in sample
pancan_mut$expression  <- mapply(function(x,y) return(as.numeric(as.character(pancan_seq[rownames(pancan_seq)==x, 
                                                                                         colnames(pancan_seq)==y ] ))),
                                 pancan_mut$gene,
                                 pancan_mut$sample) 
table(pancan_mut$cancer)
# BLCA   BRCA   CESC   COAD    GBM   HNSC   KIRC   LAML    LGG   LUAD   LUSC 
# 36103  80282  38010 102983  10375  47663  24181   1892   9185  65938  60590 
# OV   PRAD   READ   SKCM   UCEC    UCS 
# 6281  11269  23575 383797   2191   8652 

pancan_mut[pancan_mut$effect=="Frame_Shift_Del","effect"] <- "FSD"
pancan_mut$effect[pancan_mut$effect=="Frame_Shift_Ins"] <- "FSI"
pancan_mut$effect[pancan_mut$effect=="In_Frame_Del"] <- "IFD"
pancan_mut$effect[pancan_mut$effect=="In_Frame_Ins"] <- "IFI"
pancan_mut$effect[pancan_mut$effect=="Missense_Mutation"] <- "MM"
pancan_mut$effect[pancan_mut$effect=="Nonsense_Mutation"] <- "NM"
sum(grepl("single exon",pancan_mut$note)&pancan_mut$nmd.detail2==T,na.rm = T)
# [1] 0
pancan_mut$nmd.detail2[grepl("single exon",pancan_mut$note)&pancan_mut$nmd.detail2==T] <- F
save(pancan_mut,file = "data/results/pancan_mut_expre.RData")

table(pancan_mut$nmd.detail2)
# FALSE    TRUE 
# 1202212   73855  
prop.table(table(pancan_mut$nmd.detail2))
# FALSE       TRUE
# 0.94925368 0.05074632
prop.table(table(pancan_mut$have.ptc))
# FALSE       TRUE 
# 0.91716954 0.08283046

# Calculate the wildtype expression----------------------
#---------------------------------------------------------
load("data/pancan_mut_matrix.RData")
load("data/pancan_rnaseq.RData")
load("data/pan_cancer_mut_anno_table.RData")
mut_tb <- mut_tb[rownames(mut_tb) %in%  pancan_mut$sample, colnames(mut_tb) %in% pancan_mut$gene]
dim(mut_tb)
#  6924 17166
dim(pancan_mut)
# [1] 912967     24

pancan_mut$num_nmd_mut <- table(pancan_mut$gene[pancan_mut$nmd.detail2==T])[
    match(pancan_mut$gene,names(table(pancan_mut$gene[pancan_mut$nmd.detail2==T])))]

pancan_mut <- pancan_mut[pancan_mut$gene %in% rownames(pancan_seq),]
dim(pancan_mut)
# [1] 1304991      29
mut_tb <- mut_tb[rownames(mut_tb) %in%  pancan_mut$sample, 
                 colnames(mut_tb) %in% pancan_mut$gene]
dim(mut_tb)
# [1] 4191 24992
save(pancan_mut,file = "data/results/pancan_mut_matrix.RData")

table(pancan_mut$cancer_abbr)
# ACC   BLCA   BRCA   CESC   COAD   DLBC   ESCA    GBM   HNSC   KIRC   KIRP    LGG 
# 7595 108367  71554  34921 100113  15330  30476   9764  80078  14706  10678  25108 
# LIHC   LUAD   LUSC     OV   PAAD   PRAD   READ   SARC   SKCM   STAD   THCA   UCEC 
# 47016 135653  59173   9759  28615  10784  22871  22567 323201 126146   8378   2138 

table(unique(pancan_mut[,c("cancer_abbr","sample")])$cancer_abbr)
# ACC BLCA BRCA CESC COAD DLBC ESCA  GBM HNSC KIRC KIRP  LGG LIHC LUAD LUSC   OV PAAD 
# 78  390  981  192  210   48  182  152  504  212  168  516  367  488  178  188  150 
# PRAD READ SARC SKCM STAD THCA UCEC 
# 332   81  257  470  272  501    7 

# firehose NMF cluster data
nmf_clust <- read.delim(file = "data/firehoseNMF2016_01_28/pancan_bestclus.txt",
                        header = T, as.is = T)
nmf_clust <- nmf_clust[nmf_clust$SampleName!="SampleName",]
sum(!unique(pancan_mut$sample)%in%nmf_clust$SampleName)
# [1] 133
sum(unique(pancan_mut$sample)%in%nmf_clust$SampleName)
# [1] 6791
table(pancan_mut$cancer_abbr[!pancan_mut$sample%in%nmf_clust$SampleName])
# BRCA  COAD  ESCA   GBM  HNSC  LUAD    OV  SARC  SKCM  THCA 
# 260   116   110    44   284   401    19   196 53310   236 
pancan_mut$nmfc <- nmf_clust$cluster[match(pancan_mut$sample,nmf_clust$SampleName)]
table(pancan_mut$nmfc )
sum(is.na(pancan_mut$nmfc ))
# [1] 54976

# CNV data
pancan_cnv <- read.delim( "data/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",header = T,as.is = T)
pancan_cnv[1:5,1:5]
rownames(pancan_cnv) <- pancan_cnv$Sample
colnames(pancan_cnv) <- gsub(colnames(pancan_cnv),pattern = "[.]",replacement = "-")
pancan_cnv <- pancan_cnv[,colnames(pancan_cnv) %in% pancan_mut$sample]
save(pancan_cnv,file="data/pancan_cnv.RData")

pancan_mut$quantile_expression <- rep(NA,nrow(pancan_mut))
pancan_mut$wt_expression_median <- rep(NA,nrow(pancan_mut))
pancan_mut$wt_expression_mean <- rep(NA,nrow(pancan_mut))
pancan_mut$n_wt <- rep(NA,nrow(pancan_mut))
  
sum(!pancan_mut$sample %in%colnames(pancan_cnv))
# [1] 64528
sum(!pancan_mut$gene %in%rownames(pancan_cnv))
# [1] 7642
for(i in 1:nrow(pancan_mut)) # calculate the mean/median wt expression and the quantile
{ 
  expression_wt <- get.wt.expression(pancan_mut$gene[i], pancan_mut$cancer_abbr[i],pancan_mut$nmfc[i])
  pancan_mut$n_wt[i] <- length(expression_wt)
  if(length(expression_wt)>=5)
  {
    pancan_mut$wt_expression_mean[i] <- mean(as.numeric(expression_wt),na.rm	=T)
    pancan_mut$wt_expression_median[i] <- median(as.numeric(expression_wt),na.rm	=T)
    pancan_mut$quantile_expression[i] <- sum(expression_wt<=pancan_mut$expression[i])/length(expression_wt)
  }
}
pancan_mut_expression <- pancan_mut
any(is.na(pancan_mut$wt_expression_median))
save(pancan_mut_expression,file = "data/results/pancan_mut_classfied_with_quantileexpression.RData")


### pancan_mut2: samples without RNA-seq data-----
# ~Classify the mutations into NMD and non-NMD----------
#------------------------------------------------------
# combine the mutation files of five cancer types
nrow(pancan_mut2)
# [1] 516118
# exclude mutations outside coding sequence
# pancan_mut2<- pancan_mut2[!(pancan_mut2$effect %in% c("Splice_Site","3'UTR","5'UTR","IGR","Intron","RNA","targeted_region")),] 

# [1] 1344589

# filter the mutations without annotations in cds db
has.anno <- which(pancan_mut2$geneId %in% cds37$geneid)
# pancan_mut2 <- pancan_mut2[has.anno,]; rm(has.anno)
length(has.anno)
nrow(pancan_mut2) 
# [1] 1343298
save(pancan_mut2,file="data/results/pancan_mut2_before_assess")
# see if the gene has correct CDS annotation: start with start codon, end to stop codon
pancan_mut2$build  <- 37

# find dominate transcript for each gene in each cancer
assess <- data.frame(gene=as.numeric(unique(pancan_mut2$geneId[has.anno]))) # extract all genes

assess$assess <- NA # 17280
# see if the gene has correct CDS annotation: start with start codon, end to stop codon
for(i in 1:nrow(assess)) 
  assess$assess[i] <-  assess.annotation(assess$gene[i])

table(assess$assess)
sum(pancan_mut2$geneId %in% assess$gene[assess$assess==0])
sum(!pancan_mut2$geneId %in% assess$gene[assess$assess==0])

pancan_mut2$classified <- T
pancan_mut2$classified[-has.anno] <- F
pancan_mut2$classified[!pancan_mut2$geneId %in% assess$gene[assess$assess==0]] <- F
table(pancan_mut2$classified)

# Classify the mutated sequences into NMD T or F----------
nmd.rl <- matrix(ncol=8, nrow = nrow(pancan_mut2))
colnames(nmd.rl) <- names(classify.nmd(pancan_mut[1,]))
cancer_type$Cancer[order(cancer_type$Number_total_mutations,decreasing = T)]
sample <- which(pancan_mut2$cancer %in% cancer_type$Cancer[order(cancer_type$Number_total_mutations,decreasing = T)][1:24])
length(sample)
# [1] 473819
pancan_mut2 <- pancan_mut2[sample,]
table(pancan_mut2$classified)
# FALSE   TRUE 
# 181045 292774
length(unique(pancan_mut2$sample[pancan_mut2$cancer_abbr=="UCEC"]))
# [1] 248
table(pancan_mut2$cancer_abbr)
# ACC   BLCA   BRCA   CESC   COAD   DLBC   ESCA    GBM   HNSC   KIRC   KIRP    LGG   LIHC   LUAD   LUSC     OV 
# 5534  17036  12144  10318  25379   1586   8058  15367  14656   4010   2139   2613   6761  50630   6132  14333 
# PAAD   PRAD   READ   SARC   SKCM   STAD   THCA   UCEC 
# 1740   1564  11388   7434  47924  20646   3704 182723
for(i in 1:nrow(pancan_mut2))
{
  if(pancan_mut2$classified[i])
  {
    if( (pancan_mut2$geneId[i] %in% tx37$GENEID & pancan_mut2$build[i]==37))
      nmd.rl[i,] <- classify.nmd(pancan_mut2[i,])
    else nmd.rl[,"note"] <- "can't find geneID in txdb"
  }
  
}

pancan_mut2 <- cbind(pancan_mut2,nmd.rl[1:nrow(pancan_mut2),-3])
save(pancan_mut2,file="panCan17_results/pancan_mut2_classified.RData")

sum(is.na(pancan_mut2$mut_nmd))
# [1] 222834
table(pancan_mut2$note[!is.na(pancan_mut2$mut_nmd)])[table(pancan_mut2$note[!is.na(pancan_mut2$mut_nmd)])>0]

colnames(pancan_mut2)[which(colnames(pancan_mut2)=="mut_nmd")] <- "nmd.detail2"
pancan_mut2$nmd.detail2[grep("splice site mutated",pancan_mut2$note)] <- NA
sum(is.na(pancan_mut2$nmd.detail2))
# [1] 232686
table(pancan_mut2$nmd.detail2)

# FALSE   TRUE 
# 226951  14182 
pancan_mut2$nmd.detail<- NA
pancan_mut2$nmd.detail[which(pancan_mut2$nmd.detail2==T)] <- "NMD"
pancan_mut2$nmd.detail[pancan_mut2$nmd.detail2==F] <- "nonNMD"
pancan_mut2$nmd.detail[(grepl("cannot find ORF",pancan_mut2$note) | is.na(pancan_mut2$PTC.Stop))&!is.na(pancan_mut2$nmd.detail2)] <- "nonstop"

table(pancan_mut2[,c("effect","nmd.detail")])
which(pancan_mut2$effect=="Silent" & pancan_mut2$nmd.detail2==T)

grep("reference unmatched",pancan_mut2$note)
# integer(0)
which(pancan_mut2$effect=="Missense_Mutation" & pancan_mut2$nmd.detail2==T) # some missense SNP are NMD
grep("non-stop",pancan_mut2$note) # should be free of wierd wildtype CDS

pancan_mut$classified <- T
pancan_mut2$expression <- NA

pancan_mut_all <- rbind(pancan_mut[,match(colnames(pancan_mut2),colnames(pancan_mut))],pancan_mut2)
head(pancan_mut_all)


pancan_mut_all[pancan_mut_all$effect=="Frame_Shift_Del","effect"] <- "FSD"
pancan_mut_all$effect[pancan_mut_all$effect=="Frame_Shift_Ins"] <- "FSI"
pancan_mut_all$effect[pancan_mut_all$effect=="In_Frame_Del"] <- "IFD"
pancan_mut_all$effect[pancan_mut_all$effect=="In_Frame_Ins"] <- "IFI"
pancan_mut_all$effect[pancan_mut_all$effect=="Missense_Mutation"] <- "MM"
pancan_mut_all$effect[pancan_mut_all$effect=="Nonsense_Mutation"] <- "NM"
sum(grepl("single exon",pancan_mut_all$note)&pancan_mut_all$nmd.detail2==T,na.rm = T)
# [1] 0

save(pancan_mut_all,file="panCan17_results/pancan_mut_all_classified.RData")

pancan_mut_all$quantile_expression <- pancan_mut_expression$quantile_expression[match(paste(pancan_mut_all$gene,pancan_mut_all$sample),
                                                                                      paste(pancan_mut_expression$gene,pancan_mut_expression$sample))]

pancan_mut_all$new.annotation <- NA
pancan_mut_all$new.annotation[pancan_mut_all$nmd.detail2==T] <- "NMD-elicit"
pancan_mut_all$new.annotation[pancan_mut_all$nmd.detail2==F&pancan_mut_all$have.ptc==T] <- "NMD-escape"
pancan_mut_all$new.annotation[pancan_mut_all$have.ptc==F&pancan_mut_all$effect%in%c("3'UTR","5'UTR","Intron","RNA","Silent")] <- "Silent"
pancan_mut_all$new.annotation[pancan_mut_all$have.ptc==F&pancan_mut_all$effect%in%c("FSI","FSD")]  <- "non-PTC frameshift"
pancan_mut_all$new.annotation[pancan_mut_all$have.ptc==F&pancan_mut_all$effect%in%c("IFI","IFD")]  <- "non-PTC Inframe"
pancan_mut_all$new.annotation[pancan_mut_all$have.ptc==F&pancan_mut_all$effect%in%c("MM","NM","Nonstop_Mutation",
                                                                                    "Splice_Site","Translation_Start_Site")]  <- "non-PTC nonsilent"

pancan_mut_expression$new.annotation <- NA
pancan_mut_expression$new.annotation[pancan_mut_expression$nmd.detail2==T] <- "NMD-elicit"
pancan_mut_expression$new.annotation[pancan_mut_expression$nmd.detail2==F&pancan_mut_expression$have.ptc==T] <- "NMD-escape"
pancan_mut_expression$new.annotation[pancan_mut_expression$have.ptc==F&pancan_mut_expression$effect%in%c("3'UTR","5'UTR","Intron","RNA","Silent")] <- "Silent"
pancan_mut_expression$new.annotation[pancan_mut_expression$have.ptc==F&pancan_mut_expression$effect%in%c("FSI","FSD")]  <- "non-PTC frameshift"
pancan_mut_expression$new.annotation[pancan_mut_expression$have.ptc==F&pancan_mut_expression$effect%in%c("IFI","IFD")]  <- "non-PTC Inframe"
pancan_mut_expression$new.annotation[pancan_mut_expression$have.ptc==F&pancan_mut_expression$effect%in%c("MM","NM","Nonstop_Mutation",
                                                                                                         "Splice_Site","Translation_Start_Site")]  <- "non-PTC nonsilent"


# add cnv data------
# CNV data
pancan_cnv <- read.delim( "data/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",header = T,as.is = T)
pancan_cnv[1:5,1:5]
rownames(pancan_cnv) <- pancan_cnv$Sample
colnames(pancan_cnv) <- gsub(colnames(pancan_cnv),pattern = "[.]",replacement = "-")
pancan_mut$cnv <- mapply(x= pancan_mut$gene, y=pancan_mut$sample,FUN =get.cnv)

get.cnv <- function(x,y) # get the CNV value of gene x in sample y
{
  if(x %in% rownames(pancan_cnv) & y %in% colnames(pancan_cnv))
  {
    return(pancan_cnv[rownames(pancan_cnv)==x,colnames(pancan_cnv)==y])
  }else{
    return(NA)
  }
}   
pancan_mut_all$cnv <- NA

for(i in 1:nrow(pancan_mut_all))
{
    pancan_mut_all$cnv[i] <- get.cnv(pancan_mut_all$gene[i],pancan_mut_all$sample[i])
}

