# read the data

pancan_mut <- read.delim( "data/PANCAN_mutation",header = T,as.is = T)
pancan_seq <- read.delim( "data/PANCAN_HiSeqV2",header = T,as.is = T)
pancan_clinic <- read.delim( "data/PANCAN_clinicalMatrix",header = T,as.is = T)
cancer_type <- read.csv("data/cancer_list.csv",header = T,as.is = T)

pancan_mut$cancer <- pancan_clinic$X_primary_disease[match(pancan_mut$sample,pancan_clinic$sampleID)]
pancan_mut$cancer_abbr <- cancer_type$Abbr[match(pancan_mut$cancer,cancer_type$Cancer)]
cancer_type$Num_all_mut_original <- table(pancan_mut$cancer_abbr)[match(cancer_type$Abbr,names(table(pancan_mut$cancer_abbr)))]

# preprocess pancan_seq
rownames(pancan_seq) <- pancan_seq$Sample
pancan_seq <- pancan_seq[,-1]
colnames(pancan_seq) <- gsub(colnames(pancan_seq),pattern = "[.]",replacement = "-")
pancan_seq[1:5,1:5]
dim(pancan_seq)
# [1] 20530 10459

# add gene entrez id to pancan_mut
library(org.Hs.eg.db) # map gene name to id by org.Hs.eg.db
e2s = toTable(org.Hs.egSYMBOL)
pancan_mut$geneId <- e2s$gene_id[match(pancan_mut$gene,e2s$symbol)]
rm(e2s)

sum(is.na(pancan_mut$geneId))
# [1] 92707
pancan_mut <- pancan_mut[!is.na(pancan_mut$geneId),]
dim(pancan_mut) #[1] 1895455      13
View(data.frame(table(pancan_mut$cancer )[order(table(pancan_mut$cancer ),decreasing = T)]))

mut_tb <-table(pancan_mut[,c("sample","gene")])
sum(mut_tb)
save(mut_tb,file="data/pan_cancer_mut_anno_table.RData")
rm(mut_tb)

sample.gene <- apply(pancan_mut[,1:5],1,paste,collapse=" ")
deplicated.mut.sample.gene <- duplicated(sample.gene)
sum(deplicated.mut.sample.gene) # [1] 89933
pancan_mut <- pancan_mut[-which(sample.gene %in% sample.gene[deplicated.mut.sample.gene]),]
dim(pancan_mut)
# [1] 1733639      13
rm(sample.gene)

# filter by presence in RNA-seq data------
pancan_mut <- pancan_mut[(pancan_mut$gene %in% rownames(pancan_seq))&(pancan_mut$sample %in% colnames(pancan_seq)),]
# [1] 1400222  14
length(unique(pancan_mut$gene))# [1] 17718
length(unique(pancan_mut$sample)) #[1] 7708

# annotate the multi-mutation in one gene from the same sample
deplicated.mut.sample.gene <- paste(pancan_mut$sample,pancan_mut$gene)[duplicated(paste(pancan_mut$sample,pancan_mut$gene))]
pancan_mut$duplicated <- F
pancan_mut$duplicated[which(paste(pancan_mut$sample,pancan_mut$gene) %in% deplicated.mut.sample.gene)] <- T
rm(deplicated.mut.sample.gene)

table(pancan_mut$duplicated)
# FALSE    TRUE 
# 1130496  269726

cancer_type$Num_all_mut_fist_filtered <- table(pancan_mut$cancer_abbr)[match(cancer_type$Abbr,names(table(pancan_mut$cancer_abbr)))]

# full pancan mutations -----
pancan_mut_all <- read.delim("data/PANCAN_mutation",header = T,as.is = T)
pancan_mut_all$cancer <- pancan_clinic$X_primary_disease[match(pancan_mut_all$sample,pancan_clinic$sampleID)]
pancan_mut_all$cancer_abbr <- cancer_type$Abbr[match(pancan_mut_all$cancer,cancer_type$Cancer)]

pancan_mut_all$geneId <- e2s$gene_id[match(pancan_mut_all$gene,e2s$symbol)]

sample.gene <- apply(pancan_mut_all[,1:5],1,paste,collapse=" ")
deplicated.mut.sample.gene <- duplicated(sample.gene)
sum(deplicated.mut.sample.gene) # [1] 92809
pancan_mut_all <- pancan_mut_all[-which(sample.gene %in% sample.gene[deplicated.mut.sample.gene]),]
dim(pancan_mut_all)
# [1] 1821109      15
rm(sample.gene)

#~ annotate the multi-mutation in one gene from the same sample-----
deplicated.mut.sample.gene <- paste(pancan_mut_all$sample,pancan_mut_all$gene)[duplicated(paste(pancan_mut_all$sample,pancan_mut_all$gene))]
pancan_mut_all$duplicated <- F
pancan_mut_all$duplicated[which(paste(pancan_mut_all$sample,pancan_mut_all$gene) %in% deplicated.mut.sample.gene)] <- T
rm(deplicated.mut.sample.gene)

pancan_mut2 <- pancan_mut_all[!(apply(pancan_mut_all[,1:5],1,paste,collapse=" ") %in% apply(pancan_mut[,1:5],1,paste,collapse=" ")),]
