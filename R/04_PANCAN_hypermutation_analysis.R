#----------------------------------
# study the hypermutations and NMD
#----------------------------------

## table of samples ##------
sample_df <- data.frame(unique(pancan_mut_all[,c("sample","cancer_abbr")]))
colnames(sample_df)[2]<-"cancer"

# the number of mutations per sample
tb <- table(pancan_mut_all$sample)
sample_df$mut <- tb[match(sample_df$sample,names(tb))]
# the number of NMD mutations per sample
tb <- table(pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T])
sample_df$NMD_mut <- tb[match(sample_df$sample,names(tb))]
# the number of NMD-escape mutations per sample
tb <- table(pancan_mut_all$sample[pancan_mut_all$new.annotation=="NMD-escape"])
sample_df$NMD_mut <- tb[match(sample_df$sample,names(tb))]
# the number of NMD-elicit mutations in TSGs per sample
tb <- table(pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T&grepl("TSG",pancan_mut_all$role_in_cancer)])
sample_df$TSG_NMD_mut <- tb[match(sample_df$sample,names(tb))]

# the number of mutations in TSGs per sample
tb <- table(pancan_mut_all$sample[grepl("TSG",pancan_mut_all$role_in_cancer)])
sample_df$TSG_mut <- tb[match(sample_df$sample,names(tb))]

sample_df[is.na(sample_df)] <- 0
sample_df$TCGA.identifier <-  substr(sample_df$sample,start = 1,stop = 12)

head(sample_df)
sample_df$hypermutated <- NA
sample_df$hypermutated[sample_df$mut<100] <- F
sample_df$hypermutated[sample_df$mut>=1000] <- T
hm_cancer <- c("BLCA","COAD","LUAD","SKCM","STAD","UCEC")
hypermutated_cancer <- c("BLCA", "COAD", "LUAD" ,"SKCM" ,"STAD", "UCEC") #c("COAD","LUAD","STAD","READ")
# ACC BLCA BRCA CESC COAD DLBC ESCA  GBM HNSC KIRC  LGG LIHC LUAD LUSC PAAD READ SARC SKCM STAD UCEC 
# 5   40   19   12   50    3    3    3   18    9    3   14   78   12    1    5    6  144   65   32 
# BLCA COAD LUAD SKCM STAD UCEC 

table(sample_df$hypermutated)
length(is.na(sample_df$hypermutated))
# FALSE  TRUE 
# 4456   346
ggplot(sample_df[sample_df$cancer%in%hm_cancer,],aes(x=(as.numeric(mut))))+
    geom_histogram(binwidth=250)+facet_wrap(~cancer,ncol=2)+
    geom_vline(xintercept = 1000,lty=2)
# define hypermutaed samples------
# hypermutated samples
hm_sample <- sample_df$sample[which((sample_df$hypermutated==T)&sample_df$cancer %in% hm_cancer) ]
length(hm_sample) # 267
# non-hypermutated samples
nonhm_sample <- sample_df$sample[which((sample_df$hypermutated==F)&sample_df$cancer %in% hm_cancer) ]
length(nonhm_sample) # 602
pancan_mut_all$hypermutated <- NA
pancan_mut_all$hypermutated[pancan_mut_all$sample%in%hm_sample] <- T
pancan_mut_all$hypermutated[pancan_mut_all$sample%in%nonhm_sample] <- F

# add the microsatellite stability data------
# read the mss annotation table-----
MS_table <- read.csv(file = "data/MS_table_10.1038_nm.4191.csv",as.is = T)

sample_df <- cbind(sample_df,MS_table[match(sample_df$TCGA.identifier,MS_table$Sample.Name),c(5,7)])

pancan_mut_all$TCGA.identifier <-  substr(pancan_mut_all$sample,start = 1,stop = 12)
pancan_mut_all$MSS <- MS_table$MOSAIC.classification[match(pancan_mut_all$TCGA.identifier,MS_table$Sample.Name)]
pancan_mut_all$MSUnstable.Loci <- MS_table$Proportion.Unstable.Loci[match(pancan_mut_all$TCGA.identifier,MS_table$Sample.Name)]
pancan_mut_all$MSS_by_pcr <- MS_table$MSI.PCR.Result[match(pancan_mut_all$TCGA.identifier,MS_table$Sample.Name)]
table(pancan_mut_all[pancan_mut_all$cancer_abbr%in%hm_cancer&pancan_mut_all$gene%in%gene_dnadamage,c("MSS","new.annotation")])

# define HM&MSI samples------
hm_mss_sample <- na.omit(sample_df$sample[(sample_df$hypermutated==T)&sample_df$MOSAIC.classification=="MSS"&sample_df$cancer %in%hm_cancer ])
msi_sample <- na.omit(sample_df$sample[sample_df$MOSAIC.classification=="MSI-H"&sample_df$cancer %in% hm_cancer ])
nonhm_mss_sample <- na.omit(sample_df$sample[(sample_df$hypermutated==F)&sample_df$MOSAIC.classification=="MSS"&sample_df$cancer %in% hm_cancer ])

pancan_mut_all$sample_hm_mss <- NA
pancan_mut_all$sample_hm_mss[pancan_mut_all$sample%in%hm_mss_sample] <- "hm_mss"
pancan_mut_all$sample_hm_mss[pancan_mut_all$sample%in%msi_sample] <- "hm_msi"
pancan_mut_all$sample_hm_mss[pancan_mut_all$sample%in%nonhm_mss_sample] <- "nonhm_mss"

# 1 gene enriched NMD in hypermutated samples----
# there are enriched of multiple mutations in one samples
hypermutated_sample_gene <- data.frame(gene=  unique(pancan_mut_all$gene[(pancan_mut_all$sample%in%hm_sample)|(pancan_mut_all$sample%in%nonhm_sample)]))

# number of nmd mutations-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)])
hypermutated_sample_gene$nmd <-tb[match(hypermutated_sample_gene$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% nonhm_sample)])
hypermutated_sample_gene$nmd_bg <-tb[match(hypermutated_sample_gene$gene,names(tb))]

hypermutated_sample_gene$nmd_norm <-hypermutated_sample_gene$nmd/sum(hypermutated_sample_gene$nmd)*100000
hypermutated_sample_gene$nmd_bg_norm <-hypermutated_sample_gene$nmd_bg/sum(hypermutated_sample_gene$nmd_bg)*100000

# number of silent mutation-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="Silent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% hm_sample)])
hypermutated_sample_gene$silent <-tb[match(hypermutated_sample_gene$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio =="Silent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% nonhm_sample)])
hypermutated_sample_gene$silent_bg <-tb[match(hypermutated_sample_gene$gene,names(tb))]

hypermutated_sample_gene$silent_norm <-hypermutated_sample_gene$silent/sum(hypermutated_sample_gene$silent)*100000
hypermutated_sample_gene$silent_bg_norm <-hypermutated_sample_gene$silent_bg/sum(hypermutated_sample_gene$silent_bg)*100000


# arrange the table-----
hypermutated_sample_gene[is.na(hypermutated_sample_gene)] <- 0

hypermutated_sample_gene$cancer_gene <- F
hypermutated_sample_gene$cancer_gene[hypermutated_sample_gene$gene %in% cancergene$Gene.Symbol] <- T
hypermutated_sample_gene$TSG <-F
hypermutated_sample_gene$TSG[hypermutated_sample_gene$gene %in% cancergene$Gene.Symbol[cancergene$Role.in.Cancer!="oncogene"]]<-T
table( hypermutated_sample_gene$TSG)
# FALSE  TRUE 
# 22674   496

# the z score for each gene for each of six cancers-----
hypermutated_sample_gene_z_score <- data.frame(gene=hypermutated_sample_gene$gene,
                                               BLCA=rep(0,nrow(hypermutated_sample_gene)),
                                               COAD=rep(0,nrow(hypermutated_sample_gene)),
                                               LUAD=rep(0,nrow(hypermutated_sample_gene)),
                                               SKCM=rep(0,nrow(hypermutated_sample_gene)),
                                               STAD=rep(0,nrow(hypermutated_sample_gene)),
                                               UCEC=rep(0,nrow(hypermutated_sample_gene)))
for(i in 1:nrow(hypermutated_sample_gene))
{
    hypermutated_sample_gene_z_score[i,2:7] <- U_pan$Z_score[match(paste(hypermutated_sample_gene_z_score$gene[i],hm_cancer),paste(U_pan$gene,U_pan$cancer))]
}

# the mean z score for each gene for six cancers-----
for(i in 1:nrow(hypermutated_sample_gene))
{
    hypermutated_sample_gene$z_score[i] <- mean(U_pan$Z_score[U_pan$gene==hypermutated_sample_gene$gene[i]&
                                                                  U_pan$cancer%in%hm_cancer],na.rm=T)
}
hypermutated_sample_gene$z_score_range <- NA
hypermutated_sample_gene$z_score_range[hypermutated_sample_gene$z_score<=0.3] <- "NMD sensitive"
hypermutated_sample_gene$z_score_range[hypermutated_sample_gene$z_score>=0.4] <- "NMD insensitive"
hypermutated_sample_gene$z_score_range[hypermutated_sample_gene$z_score<0.4&hypermutated_sample_gene$z_score>0.3] <- "intermediate"
hypermutated_sample_gene$z_score_range[is.nan(hypermutated_sample_gene$z_score)] <- "no z-score"

tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene$nmd_COAD <-tb[match(hypermutated_sample_gene$gene,names(tb))]

tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene$nmd_STAD <-tb[match(hypermutated_sample_gene$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="UCEC"])
hypermutated_sample_gene$nmd_UCEC <-tb[match(hypermutated_sample_gene$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="LUAD"])
hypermutated_sample_gene$nmd_LUAD <-tb[match(hypermutated_sample_gene$gene,names(tb))]

tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="SKCM"])
hypermutated_sample_gene$nmd_SKCM <-tb[match(hypermutated_sample_gene$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="BLCA"])
hypermutated_sample_gene$nmd_BLCA <-tb[match(hypermutated_sample_gene$gene,names(tb))]

# compare NMD mutation and silent mutation-----
hypermutated_sample_gene$p.value_ns<- NA
hypermutated_sample_gene$ratio_ns <- NA
for(i in 1:nrow(hypermutated_sample_gene)) #compare  nmd/silence between hyper/non-hyper
{
    m <- matrix(c(hypermutated_sample_gene$nmd_norm[i],hypermutated_sample_gene$nmd_bg_norm[i],
                  hypermutated_sample_gene$silent_norm[i],hypermutated_sample_gene$silent_bg_norm[i]),nrow = 2,byrow = T,
                dimnames =list(c("hypermutated", "normal"),
                               c("NMD", "silent")))
    if(!any(is.na(m)))
    {
        hypermutated_sample_gene$p.value_ns[i] <- fisher.test(m)$p.value
        hypermutated_sample_gene$ratio_ns[i] <- fisher.test(m)$estimate
    }
}
hypermutated_sample_gene$ratio_ns_noninf <- (hypermutated_sample_gene$nmd_norm+1)/(1+hypermutated_sample_gene$nmd_bg_norm)*
    (hypermutated_sample_gene$silent_bg_norm+1)/(1+hypermutated_sample_gene$silent_norm)

hypermutated_sample_gene$p.adj <- p.adjust(hypermutated_sample_gene$p.value_ns,method = "BH")
index <- (hypermutated_sample_gene$ratio_ns_noninf>1&
                      hypermutated_sample_gene$p.value_ns<0.005&hypermutated_sample_gene$z_score<=0.3)
write.csv(hypermutated_sample_gene,file = "data/table/genelist_sixCancers.csv",row.names = F)               

# connect gene lists with go------
hm_sample_gene_go <- select(org.Hs.eg.db,
                            keys = as.character(hypermutated_sample_gene$gene),
                            columns=c("ENTREZID","SYMBOL","GENENAME","GO"),
                            keytype="SYMBOL")
hm_sample_gene_go <- hm_sample_gene_go[,colnames(hm_sample_gene_go)!="EVIDENCE"]
hm_sample_gene_go <- unique(hm_sample_gene_go)
hm_sample_gene_go$TERM  <- select(GO.db,
                                  keys = as.character(hm_sample_gene_go$GO),
                                  columns=c("GOID","TERM"),
                                  keytype="GOID")[,2]
hm_sample_gene_go <- cbind(hypermutated_sample_gene[match(hm_sample_gene_go$SYMBOL,hypermutated_sample_gene$gene),],hm_sample_gene_go)
View(hm_sample_gene_go)

write.csv(hm_sample_gene_go,file = "PANCAN_results/supplementary tables/20161129hm_sample_gene_go_six Cancer.csv",row.names = F)
# six cancers: GOEA------
# this one works
gene.universe<- rep(0,length(unique(hypermutated_sample_gene$gene)))
names(gene.universe) <- unique(hypermutated_sample_gene$gene)
index <- which(hypermutated_sample_gene$ratio_ns_noninf>1&
                   hypermutated_sample_gene$p.adj<0.05&
                   hypermutated_sample_gene$z_score<=0.3)
gene.universe[names(gene.universe)%in%hypermutated_sample_gene$gene[index]] <-1
table(gene.universe)
#BP-------
GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata, classicFisher = resultFisher,weight=resultweight,topNodes = 1000 ,numChar = 120)
write.csv(tb_hypermutated_goea,file = "data/table/NMDtoSilent_BP_sixCancers2.csv",row.names = F)               

# MF------
GOdata <- new("topGOdata", ontology ="MF", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata, classicFisher = resultFisher,weight=resultweight,topNodes = 1000 ,numChar = 120)
write.csv(tb_hypermutated_goea,file = "data/table/NMDtoSilent_MF_sixCancers.csv",row.names = F)               

## nmd to nonPTC nonsilent----
gene.universe<- rep(0,length(unique(hypermutated_sample_gene$gene[hypermutated_sample_gene$z_score<=0.3])))
names(gene.universe) <- unique(hypermutated_sample_gene$gene[hypermutated_sample_gene$z_score<=0.3])
index <- which(hypermutated_sample_gene$ratio_no_scaled>1&
                       hypermutated_sample_gene$p.adj<0.05&
                   hypermutated_sample_gene$z_score<=0.3)
gene.universe[names(gene.universe)%in%hypermutated_sample_gene$gene[index]] <-1
table(gene.universe)
#--------
GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea<- GenTable(GOdata, weight=resultweight,classicFisher = resultFisher,topNodes = 1000 ,numChar = 120)
write.csv(tb_hypermutated_goea,
          file = "data/table/NMDtononPTCnonSilent_BP_sixCancers.csv",row.names = F)               

GOdata <- new("topGOdata", ontology ="MF", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea_no <- GenTable(GOdata, weight=resultweight,classicFisher = resultFisher,topNodes = 1000 ,numChar = 120)
View(tb_hypermutated_goea[tb_hypermutated_goea$classicFisher<0.05&tb_hypermutated_goea$weight<0.8,])
write.csv(tb_hypermutated_goea,
          file = "data/table/NMDtononPTCnonSilent_MF_sixCancers.csv",row.names = F)               

# ~~~repeat the analysis for COAD----
# 1 gene enriched NMD in hypermutated samples----
# there are enriched of multiple mutations in one samples
hypermutated_sample_gene2 <- data.frame(gene=  unique(pancan_mut_all$gene[((pancan_mut_all$sample%in%hm_sample)|
                                                                              (pancan_mut_all$sample%in%nonhm_sample))&pancan_mut_all$cancer_abbr=="COAD"]))

# number of nmd mutations-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene2$nmd <-tb[match(hypermutated_sample_gene2$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% nonhm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene2$nmd_bg <-tb[match(hypermutated_sample_gene2$gene,names(tb))]

hypermutated_sample_gene2$nmd_norm <-hypermutated_sample_gene2$nmd/sum(hypermutated_sample_gene2$nmd,na.rm = T)*10000
hypermutated_sample_gene2$nmd_bg_norm <-hypermutated_sample_gene2$nmd_bg/sum(hypermutated_sample_gene2$nmd_bg,na.rm = T)*10000

# number of silent mutation-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="Silent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene2$silent <-tb[match(hypermutated_sample_gene2$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio =="Silent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% nonhm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene2$silent_bg <-tb[match(hypermutated_sample_gene2$gene,names(tb))]

hypermutated_sample_gene2$silent_norm <-hypermutated_sample_gene2$silent/sum(hypermutated_sample_gene2$silent,na.rm = T)*10000
hypermutated_sample_gene2$silent_bg_norm <-hypermutated_sample_gene2$silent_bg/sum(hypermutated_sample_gene2$silent_bg,na.rm = T)*10000

# number of non-PTC nonsilent SN mutation-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="non-PTC nonsilent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene2$nonPTCnonsilent <-tb[match(hypermutated_sample_gene2$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="non-PTC nonsilent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% nonhm_sample)&pancan_mut_all$cancer_abbr=="COAD"])
hypermutated_sample_gene2$nonPTCnonsilent_bg <-tb[match(hypermutated_sample_gene2$gene,names(tb))]

hypermutated_sample_gene2$nonPTCnonsilent_norm <-hypermutated_sample_gene2$nonPTCnonsilent/sum(hypermutated_sample_gene2$nonPTCnonsilent,na.rm = T)*10000
hypermutated_sample_gene2$nonPTCnonsilent_bg_norm <-hypermutated_sample_gene2$nonPTCnonsilent_bg/sum(hypermutated_sample_gene2$nonPTCnonsilent_bg,na.rm = T)*10000

head(hypermutated_sample_gene2)

# arrange the table-----
hypermutated_sample_gene2[is.na(hypermutated_sample_gene2)] <- 0

hypermutated_sample_gene2$cancer_gene <- F
hypermutated_sample_gene2$cancer_gene[hypermutated_sample_gene2$gene %in% cancergene$Gene.Symbol] <- T
hypermutated_sample_gene2$TSG <-F
hypermutated_sample_gene2$TSG[hypermutated_sample_gene2$gene %in% cancergene$Gene.Symbol[cancergene$Role.in.Cancer!="oncogene"]]<-T
table( hypermutated_sample_gene2$TSG)

# the mean z score for each gene for six cancers-----
hypermutated_sample_gene2$z_score<- NA
for(i in 1:nrow(hypermutated_sample_gene2))
{
    if(any(U_pan$gene==hypermutated_sample_gene2$gene[i]&
           U_pan$cancer=="COAD"))
    {
        hypermutated_sample_gene2$z_score[i] <- U_pan$Z_score[U_pan$gene==hypermutated_sample_gene2$gene[i]&
                                                                  U_pan$cancer=="COAD"]
    }
}
hypermutated_sample_gene2$z_score_range <- NA
hypermutated_sample_gene2$z_score_range[hypermutated_sample_gene2$z_score<=0.3] <- "NMD sensitive"
hypermutated_sample_gene2$z_score_range[hypermutated_sample_gene2$z_score>=0.4] <- "NMD insensitive"
hypermutated_sample_gene2$z_score_range[hypermutated_sample_gene2$z_score<0.4&hypermutated_sample_gene2$z_score>0.3] <- "intermediate"
hypermutated_sample_gene2$z_score_range[is.nan(hypermutated_sample_gene2$z_score)] <- "no z-score"

hypermutated_sample_gene2$wt_mean <- NA
for(i in 1:nrow(hypermutated_sample_gene2))
{
    if(any(U_pan$gene==hypermutated_sample_gene2$gene[i]&
           U_pan$cancer=="COAD"))
    hypermutated_sample_gene2$wt_mean[i] <- U_pan$bg_expression_mean[U_pan$gene==hypermutated_sample_gene2$gene[i]&
                                                                             U_pan$cancer=="COAD"]
}

write.csv(hypermutated_sample_gene2,file = "Pdata/table/genelist_COAD.csv",row.names = F)               

# compare NMD mutation and silent mutation-----
hypermutated_sample_gene2$p.value_ns<- NA
hypermutated_sample_gene2$ratio_ns <- NA
for(i in 1:nrow(hypermutated_sample_gene2)) #compare  nmd/silence between hyper/non-hyper
{
    m <- matrix(c(hypermutated_sample_gene2$nmd_norm[i],hypermutated_sample_gene2$nmd_bg_norm[i],
                  hypermutated_sample_gene2$silent_norm[i],hypermutated_sample_gene2$silent_bg_norm[i]),nrow = 2,byrow = T,
                dimnames =list(c("hypermutated", "normal"),
                               c("NMD", "silent")))
    if(!any(is.na(m)))
    {
        hypermutated_sample_gene2$p.value_ns[i] <- fisher.test(m)$p.value
        hypermutated_sample_gene2$ratio_ns[i] <- fisher.test(m)$estimate
    }
}
hypermutated_sample_gene2$ratio_ns_noninf <- (hypermutated_sample_gene2$nmd_norm+1)/(1+hypermutated_sample_gene2$nmd_bg_norm)*
    (hypermutated_sample_gene2$silent_bg_norm+1)/(1+hypermutated_sample_gene2$silent_norm)

# compare NMD mutation and nonPTCnonsilent SN mutation-----
hypermutated_sample_gene2$p.value_no<- NA
hypermutated_sample_gene2$ratio_no <- NA
for(i in 1:nrow(hypermutated_sample_gene2)) #compare  nmd/silence between hyper/non-hyper
{
    m <- matrix(c(hypermutated_sample_gene2$nmd_norm[i],hypermutated_sample_gene2$nmd_bg_norm[i],
                  hypermutated_sample_gene2$nonPTCnonsilent_norm[i],hypermutated_sample_gene2$nonPTCnonsilent_bg_norm[i]),nrow = 2,byrow = T,
                dimnames =list(c("hypermutated", "normal"),
                               c("NMD", "silent")))
    if(!any(is.na(m)))
    {
        hypermutated_sample_gene2$p.value_no[i] <- fisher.test(m)$p.value
        hypermutated_sample_gene2$ratio_no[i] <- fisher.test(m)$estimate
    }
}
hypermutated_sample_gene2$ratio_no_scaled <- (hypermutated_sample_gene2$nmd_norm+1)/(hypermutated_sample_gene2$nmd_bg_norm+1)*
    (hypermutated_sample_gene2$nonPTCnonsilent_bg_norm+1)/
    (hypermutated_sample_gene2$nonPTCnonsilent_norm+1)

View(hypermutated_sample_gene2[(hypermutated_sample_gene2$p.value_ns<0.05|hypermutated_sample_gene2$p.value_no<0.05)&
                                  hypermutated_sample_gene2$z_score_range=="NMD sensitive",])
View(hypermutated_sample_gene2)
index <- (hypermutated_sample_gene2$ratio_ns_noninf>1&
              hypermutated_sample_gene2$p.value_ns<0.005&hypermutated_sample_gene2$z_score<=0.3)

# add entrez id to improve the precision----
symbol2id <- select(org.Hs.eg.db,keys =as.character(hypermutated_sample_gene2$gene),
                    columns = c("ENTREZID","GENENAME","ALIAS","SYMBOL"),keytype ="ALIAS" )
length(unique(symbol2id$ALIAS[which(duplicated(symbol2id$ALIAS))]))
alas <- unique(symbol2id$ALIAS[symbol2id$ALIAS == symbol2id$SYMBOL & symbol2id$ALIAS%in% (symbol2id$ALIAS[which(duplicated(symbol2id$ALIAS))])])

symbol2id2 <- symbol2id[-which(symbol2id$ALIAS%in% alas&symbol2id$ALIAS!=symbol2id$SYMBOL),]
alas <- unique(symbol2id2$ALIAS[which(duplicated(symbol2id2$ALIAS))])
View(symbol2id2[symbol2id2$ALIAS %in% alas,])
hypermutated_sample_gene2<- cbind(hypermutated_sample_gene2[match(symbol2id2$ALIAS,hypermutated_sample_gene2$gene),],symbol2id2[,c("ENTREZID","SYMBOL")])
write.csv(hypermutated_sample_gene2,file = "data/table/genelist_COAD.csv",row.names = F)               

# connect gene lists with go------
hm_sample_gene_go <- select(org.Hs.eg.db,
                                keys = as.character(hypermutated_sample_gene2$ENTREZID[index]),
                            columns=c("ENTREZID","SYMBOL","GENENAME","GO"),
                            keytype="ENTREZID")
hm_sample_gene_go <- hm_sample_gene_go[,colnames(hm_sample_gene_go)!="EVIDENCE"]
hm_sample_gene_go <- unique(hm_sample_gene_go)
hm_sample_gene_go$TERM  <- select(GO.db,
                                  keys = as.character(hm_sample_gene_go$GO),
                                  columns=c("GOID","TERM"),
                                  keytype="GOID")[,2]
hm_sample_gene_go <- cbind(hypermutated_sample_gene2[match(hm_sample_gene_go$SYMBOL,hypermutated_sample_gene2$gene),],hm_sample_gene_go)
View(hm_sample_gene_go)

# GOEA------
# this one works
gene.universe<- rep(0,length(unique(hypermutated_sample_gene2$gene[hypermutated_sample_gene2$z_score<=0.3])))
names(gene.universe) <- unique(hypermutated_sample_gene2$gene[hypermutated_sample_gene2$z_score<=0.3])
index <- which(hypermutated_sample_gene2$ratio_ns_noninf>1&
              hypermutated_sample_gene2$p.value_ns<0.05&
                  (hypermutated_sample_gene2$z_score<=0.3))
gene.universe[names(gene.universe)%in%hypermutated_sample_gene2$gene[index]] <-1
table(gene.universe)
# gene.universe
# 0     1 
# 16526    78 


GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                topNodes = 1000 ,numChar = 120)
View(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,])
write.csv(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,],file = "PANCAN_results/hypermutation/20161116GO/NMDtoSilent_BP_COAD.csv",row.names = F)               

GOdata <- new("topGOdata", ontology ="MF", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test",annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)
View(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,])
write.csv(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,],file = "PANCAN_results/hypermutation/20161116GO/NMDtoSilent_MF_COAD.csv",row.names = F)               

#nmd to others--------
gene.universe<- rep(0,length(unique(hypermutated_sample_gene2$gene[hypermutated_sample_gene2$z_score<=0.3])))
names(gene.universe) <-  unique(hypermutated_sample_gene2$gene[hypermutated_sample_gene2$z_score<=0.3])
index <- which(hypermutated_sample_gene2$ratio_no_scaled>1&
                   hypermutated_sample_gene2$p.value_no<0.05& # try 0.001 or 0.05
                   (hypermutated_sample_gene2$z_score<=0.3))

gene.universe[names(gene.universe)%in%hypermutated_sample_gene2$gene[index]] <-1
table(gene.universe)
#-BP
GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)

write.csv(tb_hypermutated_goea,file = "PANCAN_results/hypermutation/20161116GO/NMDtoOtherMut_BP_COAD.csv",row.names = F)               
#-MF
GOdata <- new("topGOdata", ontology ="MF", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)

write.csv(tb_hypermutated_goea,file = "data/table/NMDtoOtherMut_MF_COAD.csv",row.names = F)               



##~~STAD------
hypermutated_sample_gene3 <- data.frame(gene=  unique(pancan_mut_all$gene[((pancan_mut_all$sample%in%hm_sample)|
                                                                               (pancan_mut_all$sample%in%nonhm_sample))&pancan_mut_all$cancer_abbr=="STAD"]))

hypermutated_sample_gene3 <- hypermutated_sample_gene3[!is.na(hypermutated_sample_gene3$gene),]
# number of nmd mutations-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene3$nmd <-tb[match(hypermutated_sample_gene3$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio=="NMD-elicit")&
                                   (pancan_mut_all$nmd.detail2==T)&
                                   (pancan_mut_all$sample %in% nonhm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene3$nmd_bg <-tb[match(hypermutated_sample_gene3$gene,names(tb))]

hypermutated_sample_gene3$nmd_norm <-hypermutated_sample_gene3$nmd/sum(hypermutated_sample_gene3$nmd,na.rm = T)*10000
hypermutated_sample_gene3$nmd_bg_norm <-hypermutated_sample_gene3$nmd_bg/sum(hypermutated_sample_gene3$nmd_bg,na.rm = T)*10000

# number of silent mutation-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="Silent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene3$silent <-tb[match(hypermutated_sample_gene3$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotatio =="Silent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% nonhm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene3$silent_bg <-tb[match(hypermutated_sample_gene3$gene,names(tb))]

hypermutated_sample_gene3$silent_norm <-hypermutated_sample_gene3$silent/sum(hypermutated_sample_gene3$silent,na.rm = T)*10000
hypermutated_sample_gene3$silent_bg_norm <-hypermutated_sample_gene3$silent_bg/sum(hypermutated_sample_gene3$silent_bg,na.rm = T)*10000

# number of non-PTC nonsilent SN mutation-----
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="non-PTC nonsilent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% hm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene3$nonPTCnonsilent <-tb[match(hypermutated_sample_gene3$gene,names(tb))]
tb<- table(pancan_mut_all$gene[(pancan_mut_all$new.annotation =="non-PTC nonsilent")&
                                   (pancan_mut_all$nmd.detail2==F)&
                                   (pancan_mut_all$sample %in% nonhm_sample)&pancan_mut_all$cancer_abbr=="STAD"])
hypermutated_sample_gene3$nonPTCnonsilent_bg <-tb[match(hypermutated_sample_gene3$gene,names(tb))]

hypermutated_sample_gene3$nonPTCnonsilent_norm <-hypermutated_sample_gene3$nonPTCnonsilent/sum(hypermutated_sample_gene3$nonPTCnonsilent,na.rm = T)*10000
hypermutated_sample_gene3$nonPTCnonsilent_bg_norm <-hypermutated_sample_gene3$nonPTCnonsilent_bg/sum(hypermutated_sample_gene3$nonPTCnonsilent_bg,na.rm = T)*10000

head(hypermutated_sample_gene3)

# arrange the table-----
hypermutated_sample_gene3[is.na(hypermutated_sample_gene3)] <- 0

hypermutated_sample_gene3$cancer_gene <- F
hypermutated_sample_gene3$cancer_gene[hypermutated_sample_gene3$gene %in% cancergene$Gene.Symbol] <- T
hypermutated_sample_gene3$TSG <-F
hypermutated_sample_gene3$TSG[hypermutated_sample_gene3$gene %in% cancergene$Gene.Symbol[cancergene$Role.in.Cancer!="oncogene"]]<-T
table( hypermutated_sample_gene3$TSG)

# the mean z score for each gene for six cancers-----
hypermutated_sample_gene3$z_score<- NA
for(i in 1:nrow(hypermutated_sample_gene3))
{
    if(any(U_pan$gene==hypermutated_sample_gene3$gene[i]&
           U_pan$cancer=="STAD"))
    {
        hypermutated_sample_gene3$z_score[i] <- U_pan$Z_score[U_pan$gene==hypermutated_sample_gene3$gene[i]&
                                                                  U_pan$cancer=="STAD"]
    }
}
hypermutated_sample_gene3$z_score_range <- NA
hypermutated_sample_gene3$z_score_range[hypermutated_sample_gene3$z_score<=0.3] <- "NMD sensitive"
hypermutated_sample_gene3$z_score_range[hypermutated_sample_gene3$z_score>=0.4] <- "NMD insensitive"
hypermutated_sample_gene3$z_score_range[hypermutated_sample_gene3$z_score<0.4&hypermutated_sample_gene3$z_score>0.3] <- "intermediate"
hypermutated_sample_gene3$z_score_range[is.nan(hypermutated_sample_gene3$z_score)] <- "no z-score"

hypermutated_sample_gene3$wt_mean <- NA
for(i in 1:nrow(hypermutated_sample_gene3))
{
    if(any(U_pan$gene==hypermutated_sample_gene3$gene[i]&
           U_pan$cancer=="STAD"))
        hypermutated_sample_gene3$wt_mean[i] <- U_pan$bg_expression_mean[U_pan$gene==hypermutated_sample_gene3$gene[i]&
                                                                             U_pan$cancer=="STAD"]
}
# compare NMD mutation and silent mutation-----
hypermutated_sample_gene3$p.value_ns<- NA
hypermutated_sample_gene3$ratio_ns <- NA
for(i in 1:nrow(hypermutated_sample_gene3)) #compare  nmd/silence between hyper/non-hyper
{
    m <- matrix(c(hypermutated_sample_gene3$nmd_norm[i],hypermutated_sample_gene3$nmd_bg_norm[i],
                  hypermutated_sample_gene3$silent_norm[i],hypermutated_sample_gene3$silent_bg_norm[i]),nrow = 2,byrow = T,
                dimnames =list(c("hypermutated", "normal"),
                               c("NMD", "silent")))
    if(!any(is.na(m)))
    {
        hypermutated_sample_gene3$p.value_ns[i] <- fisher.test(m)$p.value
        hypermutated_sample_gene3$ratio_ns[i] <- fisher.test(m)$estimate
    }
}
hypermutated_sample_gene3$ratio_ns_noninf <- (hypermutated_sample_gene3$nmd_norm+1)/(1+hypermutated_sample_gene3$nmd_bg_norm)*
    (hypermutated_sample_gene3$silent_bg_norm+1)/(1+hypermutated_sample_gene3$silent_norm)

# compare NMD mutation and nonPTCnonsilent SN mutation-----
hypermutated_sample_gene3$p.value_no<- NA
hypermutated_sample_gene3$ratio_no <- NA
for(i in 1:nrow(hypermutated_sample_gene3)) #compare  nmd/silence between hyper/non-hyper
{
    m <- matrix(c(hypermutated_sample_gene3$nmd_norm[i],hypermutated_sample_gene3$nmd_bg_norm[i],
                  hypermutated_sample_gene3$nonPTCnonsilent_norm[i],hypermutated_sample_gene3$nonPTCnonsilent_bg_norm[i]),nrow = 2,byrow = T,
                dimnames =list(c("hypermutated", "normal"),
                               c("NMD", "silent")))
    if(!any(is.na(m)))
    {
        hypermutated_sample_gene3$p.value_no[i] <- fisher.test(m)$p.value
        hypermutated_sample_gene3$ratio_no[i] <- fisher.test(m)$estimate
    }
}
hypermutated_sample_gene3$ratio_no_scaled <- (hypermutated_sample_gene3$nmd_norm+1)/(hypermutated_sample_gene3$nmd_bg_norm+1)*
    (hypermutated_sample_gene3$nonPTCnonsilent_bg_norm+1)/
    (hypermutated_sample_gene3$nonPTCnonsilent_norm+1)

View(hypermutated_sample_gene3[(hypermutated_sample_gene3$p.value_ns<0.05|hypermutated_sample_gene3$p.value_no<0.05)&
                                   hypermutated_sample_gene3$z_score_range=="NMD sensitive",])
View(hypermutated_sample_gene3)
index <- (hypermutated_sample_gene3$ratio_ns_noninf>1&
              hypermutated_sample_gene3$p.value_ns<0.005&hypermutated_sample_gene3$z_score<=0.3)


# connect gene lists with go------
hm_sample_gene_go <- select(org.Hs.eg.db,
                            keys = as.character(hypermutated_sample_gene3$ENTREZID[index]),
                            columns=c("ENTREZID","SYMBOL","GENENAME","GO"),
                            keytype="ENTREZID")
hm_sample_gene_go <- hm_sample_gene_go[,colnames(hm_sample_gene_go)!="EVIDENCE"]
hm_sample_gene_go <- unique(hm_sample_gene_go)
hm_sample_gene_go$TERM  <- select(GO.db,
                                  keys = as.character(hm_sample_gene_go$GO),
                                  columns=c("GOID","TERM"),
                                  keytype="GOID")[,2]
hm_sample_gene_go <- cbind(hypermutated_sample_gene3[match(hm_sample_gene_go$SYMBOL,hypermutated_sample_gene3$gene),],hm_sample_gene_go)
View(hm_sample_gene_go)

write.csv(hm_sample_gene_go,file = "data/table/hm_sample_gene_go.csv",row.names = F)
# GOEA------

gene.universe<- rep(0,length(unique(hypermutated_sample_gene3$gene[hypermutated_sample_gene3$z_score<=0.3])))
names(gene.universe) <- unique(hypermutated_sample_gene3$gene[hypermutated_sample_gene3$z_score<=0.3])
index <- which(hypermutated_sample_gene3$ratio_ns_noninf>1&
                   hypermutated_sample_gene3$p.value_ns<0.05&
                   (hypermutated_sample_gene3$z_score<=0.3))
gene.universe[names(gene.universe)%in%hypermutated_sample_gene3$gene[index]] <-1
table(gene.universe)

GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)
View(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,])
write.csv(tb_hypermutated_goea,file = "data/table/NMDtoSilent_BP_STAD.csv",row.names = F)               


GOdata <- new("topGOdata", ontology ="MF", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)
View(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,])
write.csv(tb_hypermutated_goea,file = "data/table/NMDtoSilent_MF_STAD.csv",row.names = F)               


#nmd to others--------
gene.universe<- rep(0,length(unique(hypermutated_sample_gene3$gene[hypermutated_sample_gene3$z_score<=0.3])))
names(gene.universe) <- unique(hypermutated_sample_gene3$gene[hypermutated_sample_gene3$z_score<=0.3])
index <- which(hypermutated_sample_gene3$ratio_no_scaled>1&
                   hypermutated_sample_gene3$p.value_no<0.05&
                   (hypermutated_sample_gene3$z_score<=0.3))

gene.universe[names(gene.universe)%in%hypermutated_sample_gene3$gene[index]] <-1
table(gene.universe)
#--------
GOdata <- new("topGOdata", ontology ="BP", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)
write.csv(tb_hypermutated_goea,file = "data/table/NMDtononPTCnonSilent_BP_STAD.csv",row.names = F)               

GOdata <- new("topGOdata", ontology ="MF", allGenes = gene.universe,geneSel = function(x) x==1, 
              description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "alias")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
tb_hypermutated_goea <- GenTable(GOdata,classicFisher = resultFisher,
                                 weight=resultweight,
                                 topNodes = 1000 ,numChar = 120)
View(tb_hypermutated_goea[tb_hypermutated_goea$Significant>=2,])
write.csv(tb_hypermutated_goea,file = "data/table/NMDtononPTCnonSilent_MF_STAD.csv",row.names = F)               


write.csv(hypermutated_sample_gene3,file = "data/table/genelist_STAD.csv",row.names = F)   


# Validate correlation between NMD and RNA-seq for LARP4B EIF5B-------
# for STAD in TCGA 

sample_stad <- sample_df[sample_df$cancer=="STAD",]
dim(sample_stad)


# PTEN-----
sample_stad$nmd_PTEN <- F
sample_stad$ptc_PTEN <- F
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="PTEN"&pancan_mut_all$nmd.detail2==T]
sample_stad$nmd_PTEN[sample_stad$sample%in%tb] <- T
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="PTEN"&pancan_mut_all$have.ptc==T]
sample_stad$ptc_PTEN[sample_stad$sample%in%tb] <- T

sample_stad$mut_PTEN <- F
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="PTEN"&pancan_mut_all$new.annotation!="Silent"]
sample_stad$mut_PTEN[sample_stad$sample%in%tb] <- T


sample_stad$expression_PTEN<- NA
tb <- pancan_seq[rownames(pancan_seq)=="PTEN",na.omit(match(sample_stad$sample,colnames(pancan_seq)))]
sample_stad$expression_PTEN[match(names(tb),sample_stad$sample)] <- as.numeric(tb)



# LARP4B
sample_stad$nmd_LARP4B <- F
sample_stad$ptc_LARP4B <- F
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="LARP4B"&pancan_mut_all$nmd.detail2==T]
sample_stad$nmd_LARP4B[sample_stad$sample%in%tb] <- T
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="LARP4B"&pancan_mut_all$have.ptc==T]
sample_stad$ptc_LARP4B[sample_stad$sample%in%tb] <- T

sample_stad$mut_LARP4B <- F
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="LARP4B"&pancan_mut_all$new.annotation!="Silent"]
sample_stad$mut_LARP4B[sample_stad$sample%in%tb] <- T

sample_stad$expression_LARP4B<- NA
tb <- pancan_seq[rownames(pancan_seq)=="LARP4B",na.omit(match(sample_stad$sample,colnames(pancan_seq)))]
sample_stad$expression_LARP4B[match(names(tb),sample_stad$sample)] <- as.numeric(tb)

wilcox.test(sample_stad$expression_LARP4B[sample_stad$nmd_LARP4B==T],  sample_stad$expression_LARP4B[sample_stad$nmd_LARP4B==F])
wilcox.test(sample_stad$expression_PTEN[sample_stad$nmd_PTEN==T],  sample_stad$expression_PTEN[sample_stad$nmd_PTEN==F])

# EIF5B-----
sample_stad$nmd_EIF5B <- F
sample_stad$ptc_EIF5B <- F
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="EIF5B"&pancan_mut_all$nmd.detail2==T]
sample_stad$nmd_EIF5B[sample_stad$sample%in%tb] <- T
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="EIF5B"&pancan_mut_all$have.ptc==T]
sample_stad$ptc_EIF5B[sample_stad$sample%in%tb] <- T
sample_stad$mut_EIF5B <- F
tb <- pancan_mut_all$sample[pancan_mut_all$gene=="EIF5B"&pancan_mut_all$new.annotation!="Silent"]
sample_stad$mut_EIF5B[sample_stad$sample%in%tb] <- T

sample_stad$cnv_EIF5B <- NA
tb <- pancan_cnv[rownames(pancan_cnv)=="EIF5B",na.omit(match(sample_stad$sample,colnames(pancan_cnv)))]
sample_stad$cnv_EIF5B[match(names(tb),sample_stad$sample)] <- as.numeric(tb)


sample_stad$expression_EIF5B<- NA
tb <- pancan_seq[rownames(pancan_seq)=="EIF5B",na.omit(match(sample_stad$sample,colnames(pancan_seq)))]
sample_stad$expression_EIF5B[match(names(tb),sample_stad$sample)] <- as.numeric(tb)

# plots----
ggplot(sample_stad,aes(y = as.numeric(mut),x=nmd_TTN)) + geom_boxplot()
    labs(x="Mutation in EIF5B",y="Expression level of EIF5B")+theme(legend.position="none")

    ggplot(sample_stad,aes(y = expression_EIF5B,x=tp53)) + geom_boxplot()
    labs(x="Mutation in EIF5B",y="Expression level of EIF5B")+theme(legend.position="none")
    

ggplot(sample_stad,aes(y = expression_EIF5B,x=paste(nmd_EIF5B,ptc_EIF5B,mut_EIF5B),fill=paste(nmd_EIF5B,ptc_EIF5B,mut_EIF5B))) + 
    geom_boxplot()+
    scale_x_discrete(breaks=c("TRUE TRUE TRUE","FALSE TRUE TRUE", "FALSE FALSE TRUE","FALSE FALSE FALSE" ),
                         labels=c( "NMD-elicit", "NMD-escape","nonsilent\nnon-PTC","no nonsilent"))+
    labs(x="Mutation in EIF5B",y="Expression level of EIF5B")+theme(legend.position="none")

ggplot(sample_stad,aes(y = expression_LARP4B,x=paste(nmd_LARP4B,ptc_LARP4B,mut_LARP4B),fill=paste(nmd_LARP4B,ptc_LARP4B,mut_LARP4B))) + 
    geom_boxplot()+
    scale_x_discrete(breaks=c("TRUE TRUE TRUE","FALSE TRUE TRUE", "FALSE FALSE TRUE","FALSE FALSE FALSE" ),
                     labels=c( "NMD-elicit", "NMD-escape","nonsilent\nnon-PTC","no nonsilent"))+
    labs(x="Mutation in LARP4B",y="Expression level of LARP4B")+theme(legend.position="none")

ggplot(sample_stad,aes(y = expression_PTEN,x=paste(nmd_PTEN,ptc_PTEN,mut_PTEN),fill=paste(nmd_PTEN,ptc_PTEN,mut_PTEN))) + 
    geom_boxplot()+
    scale_x_discrete(breaks=c("TRUE TRUE TRUE","FALSE TRUE TRUE", "FALSE FALSE TRUE","FALSE FALSE FALSE" ),
                     labels=c( "NMD-elicit", "NMD-escape","nonsilent\nnon-PTC","no nonsilent"))+
    labs(x="Mutation in PTEN",y="Expression level of PTEN")+theme(legend.position="none")


ggplot(sample_stad,aes(y = as.numeric(mut),x=paste(ptc_EIF5B,ptc_LARP4B)) )+ geom_boxplot()+
    scale_x_discrete(breaks=c("TRUE TRUE","FALSE TRUE", "TRUE FALSE","FALSE FALSE" ),
                     labels=c( "both", "LARP4B only","EIF5B only","neither"))+
    labs(x="NMD-elicit mutation in LARP4B and EIF5B",y="Mutation frequency")

    scale_x_discrete(breaks=c("TRUE TRUE TRUE","FALSE TRUE TRUE", "FALSE FALSE TRUE","FALSE FALSE FALSE" ),
                     labels=c( "NMD-elicit", "NMD-escape","nonsilent\nnon-PTC","no nonsilent"))+
    labs(x="Mutation in LARP4B",y="Expression level of LARP4B")+theme(legend.position="none")



ggplot(sample_stad,
       aes(y=as.numeric(mut),x = expression_EIF5B,col=paste(ptc_EIF5B,mut_EIF5B),shape=MOSAIC.classification)) + 
    geom_point(alpha=0.5)+
    geom_hline(yintercept = 1000,lty=2) + 
    scale_color_discrete(name="mutation in EIF5B",
                        breaks=c("FALSE FALSE", "FALSE TRUE", "TRUE TRUE"),
                        labels=c("no nonsilent", "nonsilent non-PTC", "PTC"))+
    labs(y="Number of mutations in sample",x="Expression level of EIF5B")

ggplot(sample_stad,aes(y=as.numeric(mut),x = expression_LARP4B,col=paste(ptc_LARP4B,mut_LARP4B),shape=MOSAIC.classification)) + 
    geom_point(alpha=0.5)+
    geom_hline(yintercept = 1000,lty=2)+ 
    scale_color_discrete(name="mutation in LARP4B",
                         breaks=c("FALSE FALSE", "FALSE TRUE", "TRUE TRUE"),
                         labels=c("no nonsilent", "nonsilent non-PTC", "PTC"))+
    labs(y="Number of mutations in sample",x="Expression level of LARP4B")

ggplot(data=sample_stad,aes(y=as.numeric(mut),x = expression_PTEN,col=paste(ptc_PTEN,mut_PTEN),shape=MOSAIC.classification)) + 
    geom_point(alpha=0.5)+
    geom_hline(yintercept = 1000,lty=2)+ 
    scale_color_discrete(name="mutation in PTEN",
                         breaks=c("FALSE FALSE", "FALSE TRUE", "TRUE TRUE"),
                         labels=c("no nonsilent", "nonsilent non-PTC", "PTC"))+
    labs(y="Number of mutations in sample",x="Expression level of PTEN")
    


ggplot(sample_stad,aes(y=as.numeric(mut),x = expression_EIF5B,col=nmd_EIF5B)) + geom_point(alpha=0.5)
 
ggplot(sample_stad,aes(size=as.numeric(mut),x = expression_LARP4B,y=expression_EIF5B,col=paste(ptc_EIF5B,ptc_LARP4B))) + geom_point(alpha=0.5)+
    theme(legend.position="bottom")+geom_vline(xintercept = median(sample_stad$expression_LARP4B,na.rm = T))+
    geom_hline(yintercept = median(sample_stad$expression_EIF5B,na.rm = T))

# --------------------------------------------------------------#
# ---------------------*large heatmap plot----------------------
# --------------------------------------------------------------#

# GO to gene------
my_keys <- keys(org.Hs.eg.db)
gene2go <-select(org.Hs.eg.db,
                 keys = as.character(my_keys),
                 columns=c("ENTREZID","SYMBOL","GENENAME","GO"),
                 keytype="ENTREZID")
dim(gene2go)

interested_go <- c("GO:0007010","GO:0016568","GO:0007050","GO:0043087","GO:0044822","GO:0006281","GO:0042981")
interested_go <- c("GO:0016568","GO:0044822","GO:0006281")

interested_go_term <- select(GO.db, keys = as.character(interested_go ),columns=c("GOID", "TERM"),keytype="GOID")[,2]
go_mut_tb <-  data.frame(sample =  unique(pancan_mut_all$sample[pancan_mut_all$cancer_abbr%in%hm_cancer&
                                                                    pancan_mut_all$sample%in%c(hm_sample,nonhm_sample)]))

for(i in 1:length(interested_go))
{
    go_mut_tb <- cbind(go_mut_tb,rep(F,nrow(go_mut_tb)))
    mmr_gene <- gene2go$SYMBOL[which(gene2go$GO==interested_go[i])] #genes in the GO
    go_mut_tb[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T&pancan_mut_all$gene%in%mmr_gene],i+1] <-T
}
colnames(go_mut_tb) <- c("sample",interested_go_term )
head(go_mut_tb)
# get the dataframe for plotting-----
go_mut_tb2 <- data.frame(sample=rep(go_mut_tb$sample,length(interested_go)),
                         GO = rep(interested_go_term,each=nrow(go_mut_tb)),
                         NMD_mutated =as.vector(as.matrix(go_mut_tb[,-1])))
go_mut_tb2 <- rbind(go_mut_tb2 ,data.frame(sample=go_mut_tb$sample,
                                           GO = rep("microsatellite stability",length(go_mut_tb$sample)),
                                           NMD_mutated =pancan_mut_all$MSS[match(go_mut_tb$sample,pancan_mut_all$sample)]))
go_mut_tb2$cancer <- pancan_mut_all$cancer_abbr[match(go_mut_tb2$sample,pancan_mut_all$sample)]

# add the mutation frequency
go_mut_tb2 <- rbind(go_mut_tb2 ,data.frame(sample=go_mut_tb$sample,
                                           GO = rep("log10 mutation freq",length(go_mut_tb$sample)),
                                           NMD_mutated =as.character(round(log10(sample_df$mut))[match(go_mut_tb$sample,sample_df$sample)]),
                                           cancer = sample_df$cancer[match(go_mut_tb$sample,sample_df$sample)]))


# add the mutation of TP53, PTEN, POLE, POLD1------
# interested_genes <- c("HDAC4","BRD3","SVIL","LRP1","MSH2","MSH6")
interested_genes <- c("LARP4B","EIF5B","YTHDC1",
                      "SMARCAD1","CHD1","CHD8","HDAC4","BRD3","TTF1",
                      "MSH2", "MSH6","ASCC3","PTEN")
for(i in 1:length(interested_genes ))
{
    go_mut_tb2_add <- data.frame(sample=go_mut_tb$sample,
                                 GO = rep(interested_genes[i],length(go_mut_tb$sample)),
                                 NMD_mutated =rep(NA,length(go_mut_tb$sample)),
                                 cancer = sample_df$cancer[match(go_mut_tb$sample,sample_df$sample)])
    nonsilent <- c("FSD","FSI","IFD","IFI","MM","NM","Splice_Site","Translation_Start_Site")
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$nmd.detail2==F&(!pancan_mut_all$effect%in%c("FSD","FSI"))&
                                                                               pancan_mut_all$effect%in% nonsilent &pancan_mut_all$gene==interested_genes[i]]] <- "other nonsilent"
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$nmd.detail2==F&pancan_mut_all$effect%in%c("FSD","FSI")&pancan_mut_all$gene==interested_genes[i]]] <- "frameshift"
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$new.annotation=="NMD-elicit"&pancan_mut_all$gene==interested_genes[i]]] <- "NMD-elicit"
    go_mut_tb2 <- rbind(go_mut_tb2 ,
                        go_mut_tb2_add)
    
}


go_mut_tb2_add <- data.frame(sample=go_mut_tb$sample,
                             GO = rep("hypermutation",length(go_mut_tb$sample)),
                             NMD_mutated =rep(NA,length(go_mut_tb$sample)),
                             cancer = sample_df$cancer[match(go_mut_tb$sample,sample_df$sample)])
go_mut_tb2_add$NMD_mutated[go_mut_tb2_add$sample%in%nonhm_sample]<-"non-hypermutated"
go_mut_tb2_add$NMD_mutated[go_mut_tb2_add$sample%in%hm_sample]<-"hypermutated"
go_mut_tb2 <- rbind(go_mut_tb2 ,go_mut_tb2_add)
go_mut_tb2$NMD_mutated <- factor(go_mut_tb2$NMD_mutated,levels=c("0","1","2","3","4","MSI-H","MSS","FALSE","TRUE","other nonsilent",
                                                                 "NMD-elicit","frameshift","hypermutated","non-hypermutated"))

go_mut_tb2$sample <- factor(go_mut_tb2$sample,levels = unique(go_mut_tb2$sample[go_mut_tb2$GO=="log10 mutation freq"][order(sample_df$mut[match(go_mut_tb2$sample[go_mut_tb2$GO=="log10 mutation freq"],sample_df$sample)],decreasing = T)]))
go_mut_tb2$sample <- factor(go_mut_tb2$sample,levels = unique(go_mut_tb2$sample[go_mut_tb2$GO=="microsatellite stability"][order(paste(go_mut_tb2$NMD_mutated[go_mut_tb2$GO=="microsatellite stability"],(10000000-sample_df$mut[match(go_mut_tb2$sample[go_mut_tb2$GO=="microsatellite stability"],sample_df$sample)])),decreasing = F)]))


go_mut_tb2$GO <- factor(go_mut_tb2$GO,levels=c(interested_go_term,interested_genes,
                                               "microsatellite stability","log10 mutation freq","hypermutation"))


# new more dense firgue----
go_mut_tb3 <- unique(go_mut_tb2[go_mut_tb2$sample%in%c(hm_sample,nonhm_sample)&go_mut_tb2$cancer%in%hm_cancer,c("sample","cancer")])
go_mut_tb3 <- data.frame(sample=go_mut_tb3$sample,
                         GO=rep("cancer",nrow(go_mut_tb3)),
                         NMD_mutated=go_mut_tb3$cancer)
go_mut_tb3 <- rbind(go_mut_tb2[go_mut_tb2$sample%in%c(hypermutated_sample,nonhypermutated_sample),-4],go_mut_tb3)
go_mut_tb3 <- go_mut_tb3[go_mut_tb3$sample%in%sample_df$sample[!is.na(sample_df$MOSAIC.classification)],]
head(go_mut_tb3)
dim(go_mut_tb3 )

samplesforplot=unique(go_mut_tb3$sample)
samplesforplot <- samplesforplot[order(paste(sample_df$cancer[match(samplesforplot,sample_df$sample)],
                                             sample_df$MOSAIC.classification[match(samplesforplot,sample_df$sample)],
                                             (100000000-sample_df$mut[match(samplesforplot,sample_df$sample)])))]

go_mut_tb3$sample <- factor(go_mut_tb3$sample,
                            levels = samplesforplot)


# choose color
levels(go_mut_tb3$NMD_mutated)
# [1] "0"                "1"                "2"                "3"                "4"                "MSI-H"           
# [7] "MSS"              "FALSE"            "TRUE"             "other nonsilent"  "frameshift"       "NMD-elicit"      
# [13] "hypermutated"     "non-hypermutated" "BLCA"             "COAD"             "LUAD"             "SKCM"            
# [19] "STAD"             "UCEC"           

color_heatmap <- rep("",20)
color_heatmap[1:5] <- brewer.pal(name = "Blues",n=9)[c(1,3,5,7,9)]# colors for expression levels
color_heatmap[6:7] <- brewer.pal(name = "Set2",n=8)[c(2,5)]# colors for expression levels
color_heatmap[8:9] <- brewer.pal(name = "Greys",n=9)[c(3,9)]# colors for expression levels
color_heatmap[c(10,11,12)] <- brewer.pal(name = "RdYlGn",n=9)[c(7,1,4)]
color_heatmap[c(13:14)] <- brewer.pal(name = "RdYlGn",n=9)[c(8,3)]
color_heatmap[c(15:20)] <- brewer.pal(name = "Accent",n=8)[1:6]

# axis.text.x=element_blank(),
ggplot(go_mut_tb3[go_mut_tb3$sample%in% samplesforplot,],
       aes(x=sample,y=GO,fill=NMD_mutated))+geom_tile(col="white")+
    theme(axis.line=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),axis.text.x=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    theme(legend.position="bottom")+
    labs(title = "")+scale_fill_manual(values=color_heatmap)+
    guides(fill = guide_legend(nrow = 3,title =""))

pancan_mut_all$new.annotation <- as.character(pancan_mut_all$new.annotation)
pancan_mut_all$new.annotation[is.na(pancan_mut_all$new.annotation)]<- "unclassified"
go_mut_sample <- data.frame(sample=levels(go_mut_tb3$sample))
go_mut_sample <- as.data.frame(table(pancan_mut_all[pancan_mut_all$sample %in% go_mut_sample$sample,
                      c("new.annotation","sample")]))
head(go_mut_sample )


go_mut_sample$sample <- factor(go_mut_sample$sample,levels=samplesforplot)


ggplot(go_mut_sample[go_mut_sample$sample%in% samplesforplot,],aes(x=sample,y=Freq,fill=new.annotation))+
    geom_bar(stat="identity",position = "stack")+labs(x="",fill="")+
    theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),legend.position="top")



ggplot(go_mut_tb3[go_mut_tb3$sample %in% sample_df$sample[sample_df$cancer=="STAD"],],aes(x=sample,y=GO,fill=NMD_mutated))+geom_tile(col="white")+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(angle=30,size=4),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    theme(legend.position="bottom")+
    labs(title = "")+scale_fill_manual(values=color_heatmap)+
    guides(fill = guide_legend(nrow = 5,title =""))

go_mut_tb2[which(go_mut_tb2$NMD_mutated==4&go_mut_tb2$cancer=="STAD"),]

write.csv(pancan_mut_all[pancan_mut_all$gene=="LARP4B"&pancan_mut_all$cancer_abbr=="STAD",],"PANCAN_results/STAD/LARP4Bmut_STAD.csv")

pancan_mut_all$n_mut_in_sample <- NA
pancan_mut_all$n_mut_in_sample <- sample_df$mut[match(pancan_mut_all$sample,sample_df$sample)]

write.csv(pancan_mut_all[pancan_mut_all$gene=="LARP4B"&pancan_mut_all$cancer_abbr=="STAD",],"PANCAN_results/STAD/LARP4Bmut_STAD.csv",row.names = F)
write.csv(pancan_mut_all[pancan_mut_all$gene=="EIF5B"&pancan_mut_all$cancer_abbr=="STAD",],"PANCAN_results/STAD/EIF5Bmut_STAD.csv",row.names = F)
write.csv(pancan_mut_all[pancan_mut_all$gene=="PTEN"&pancan_mut_all$cancer_abbr=="STAD",],"PANCAN_results/STAD/PTENmut_STAD.csv",row.names = F)


write.csv(pancan_mut_all[pancan_mut_all$gene=="LARP4B",],"PANCAN_results/STAD/LARP4Bmut_allcancers.csv",row.names = F)
write.csv(pancan_mut_all[pancan_mut_all$gene=="EIF5B",],"PANCAN_results/STAD/EIF5Bmut_allcancers.csv",row.names = F)


# nmf

pancan_mut_all$nmf_cluster <- nmf_clust$cluster[match(pancan_mut_all$sample,nmf_clust$SampleName)]

# --------------------------------------------------------------#
# ---------------------*large heatmap plot for STAD----------------------
# --------------------------------------------------------------#

# GO to gene------
my_keys <- keys(org.Hs.eg.db)
gene2go <-select(org.Hs.eg.db,
                 keys = as.character(my_keys),
                 columns=c("ENTREZID","SYMBOL","GENENAME","GO"),
                 keytype="ENTREZID")
dim(gene2go)
interested_go <- c("GO:0016568","GO:0044822","GO:0006281")

interested_go_term <- select(GO.db, keys = as.character(interested_go ),columns=c("GOID", "TERM"),keytype="GOID")[,2]
go_mut_tb <-  data.frame(sample =  unique(pancan_mut_all$sample[pancan_mut_all$cancer_abbr=="STAD"]))

for(i in 1:length(interested_go))
{
    go_mut_tb <- cbind(go_mut_tb,rep(F,nrow(go_mut_tb)))
    mmr_gene <- gene2go$SYMBOL[which(gene2go$GO==interested_go[i])] #genes in the GO
    go_mut_tb[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$nmd.detail2==T&pancan_mut_all$gene%in%mmr_gene],i+1] <-T
}
colnames(go_mut_tb) <- c("sample",interested_go_term )
head(go_mut_tb)
# get the dataframe for plotting-----
go_mut_tb2 <- data.frame(sample=rep(go_mut_tb$sample,length(interested_go)),
                         GO = rep(interested_go_term,each=nrow(go_mut_tb)),
                         NMD_mutated =as.vector(as.matrix(go_mut_tb[,-1])))
go_mut_tb2 <- rbind(go_mut_tb2 ,data.frame(sample=go_mut_tb$sample,
                                           GO = rep("microsatellite stability",length(go_mut_tb$sample)),
                                           NMD_mutated =pancan_mut_all$MSS[match(go_mut_tb$sample,pancan_mut_all$sample)]))
go_mut_tb2$cancer <- pancan_mut_all$cancer_abbr[match(go_mut_tb2$sample,pancan_mut_all$sample)]

# add the mutation frequency
go_mut_tb2 <- rbind(go_mut_tb2 ,data.frame(sample=go_mut_tb$sample,
                                           GO = rep("log10 mutation freq",length(go_mut_tb$sample)),
                                           NMD_mutated =as.character(round(log10(sample_df$mut))[match(go_mut_tb$sample,sample_df$sample)]),
                                           cancer = sample_df$cancer[match(go_mut_tb$sample,sample_df$sample)]))


# add the mutation of TP53, PTEN, POLE, POLD1------
interested_genes <- c("LARP4B","EIF5B","PTEN")
for(i in 1:length(interested_genes ))
{
    go_mut_tb2_add <- data.frame(sample=go_mut_tb$sample,
                                 GO = rep(interested_genes[i],length(go_mut_tb$sample)),
                                 NMD_mutated =rep(NA,length(go_mut_tb$sample)),
                                 cancer = sample_df$cancer[match(go_mut_tb$sample,sample_df$sample)])
    nonsilent <- c("FSD","FSI","IFD","IFI","MM","NM","Splice_Site","Translation_Start_Site")
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$nmd.detail2==F&(!pancan_mut_all$effect%in%c("FSD","FSI"))&
                                                                               pancan_mut_all$effect%in% nonsilent &pancan_mut_all$gene==interested_genes[i]]] <- "other nonsilent"
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$nmd.detail2==F&pancan_mut_all$effect%in%c("FSD","FSI")&pancan_mut_all$gene==interested_genes[i]]] <- "frameshift"
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% pancan_mut_all$sample[pancan_mut_all$new.annotation=="NMD-elicit"&pancan_mut_all$gene==interested_genes[i]]] <- "NMD-elicit"
    
    go_mut_tb2 <- rbind(go_mut_tb2 ,
                        go_mut_tb2_add)
    
}


table(go_mut_tb2$NMD_mutated)

# rearrannge the position of x/y axis-----

# table(go_mut_tb2$GO)

go_mut_tb2$NMD_mutated <- factor(as.character(go_mut_tb2$NMD_mutated),
                                 levels=c("0","1","2","3","4","MSI-H","MSS","FALSE","TRUE","other nonsilent",
                                          "NMD-elicit","frameshift"))

go_mut_tb2$sample <- factor(go_mut_tb2$sample,levels = unique(go_mut_tb2$sample[go_mut_tb2$GO=="log10 mutation freq"][order(sample_df$mut[match(go_mut_tb2$sample[go_mut_tb2$GO=="log10 mutation freq"],sample_df$sample)],decreasing = T)]))
go_mut_tb2$sample <- factor(go_mut_tb2$sample,levels = unique(go_mut_tb2$sample[go_mut_tb2$GO=="microsatellite stability"][order(paste(go_mut_tb2$NMD_mutated[go_mut_tb2$GO=="microsatellite stability"],(10000000-sample_df$mut[match(go_mut_tb2$sample[go_mut_tb2$GO=="microsatellite stability"],sample_df$sample)])),decreasing = F)]))


go_mut_tb2$GO <- factor(go_mut_tb2$GO,levels=c(interested_go_term,interested_genes,
                                               "microsatellite stability","log10 mutation freq",
                                               "hypermutation"))


# new more dense firgue----
go_mut_tb3 <- go_mut_tb2
head(go_mut_tb3)
dim(go_mut_tb3 )

samplesforplot=unique(go_mut_tb3$sample)
samplesforplot <- samplesforplot[order(paste(sample_df$cancer[match(samplesforplot,sample_df$sample)],
                                             (100000000-sample_df$mut[match(samplesforplot,sample_df$sample)]),
                                             sample_df$MOSAIC.classification[match(samplesforplot,sample_df$sample)]))]

go_mut_tb3$sample <- factor(go_mut_tb3$sample,
                            levels = samplesforplot)


# choose color
levels(go_mut_tb3$NMD_mutated)
tiff("figures/Figure3b_new.tiff", width =9, height = 2.7, units = 'in', res = 500)

ggplot(go_mut_tb3[go_mut_tb3$sample%in% samplesforplot,],
       aes(x=sample,y=GO,fill=NMD_mutated))+geom_tile(col="white")+
    theme(axis.line=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),axis.text.x=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    theme(legend.position="none")+
    labs(title = "")+scale_fill_manual(values=color_heatmap[1:12] )+
    guides(fill = guide_legend(title ="",byrow = F))
dev.off()

pancan_mut_all$effect2 <- NA
pancan_mut_all$effect2[pancan_mut_all$effect%in%c("FSD","FSI")] <- "FS"
pancan_mut_all$effect2[pancan_mut_all$effect%in%c("IFD","IFI")] <- "IF"
pancan_mut_all$effect2[pancan_mut_all$effect%in%c("MM","NM")] <- "SN"
pancan_mut_all$effect2[pancan_mut_all$effect%in%c("RNA","Silent")] <- "silent"
pancan_mut_all$effect2[is.na(pancan_mut_all$effect2)] <- "others"

go_mut_sample <- data.frame(sample=levels(go_mut_tb3$sample))
go_mut_sample <- as.data.frame(table(pancan_mut_all[pancan_mut_all$sample %in% go_mut_sample$sample,
                                                    c("effect2","sample")]))
head(go_mut_sample )


go_mut_sample$sample <- factor(go_mut_sample$sample,levels=samplesforplot)

tiff("figures/FigureS32.tiff", width =9, height = 3.5, units = 'in', res = 500)

ggplot(go_mut_sample[go_mut_sample$sample%in% samplesforplot,],aes(x=sample,y=Freq,fill=effect2))+
    geom_bar(stat="identity",position = "stack")+labs(x="",fill="")+
    theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),legend.position="top")

dev.off()
