# validate the findings in a new gastric cancer data------

#hg19 == GRCh37
mutation_list<- read.csv(file = "data/gastric_validation_data/ng.2983-S2_mutations.csv",
                         header = T, stringsAsFactors=FALSE)
gcv_clinic <- read.csv(file = "data/gastric_validation_data/ng.2983-S2_clinic.csv",
                       header = T,as.is=T)


gcv_mut <- data.frame(sample=mutation_list$Sample.code,
                            chr=mutation_list$Chromosome,
                            start=mutation_list$Start,
                            end=mutation_list$End,
                            reference = mutation_list$Reference.Base,
                            alt = mutation_list$Mutant.allele,
                            gene=mutation_list$gene,
                            effect = mutation_list$Mutation.Type,
                            zygosity=mutation_list$Mutant.zygosity,
                            Amino_Acid_Change = mutation_list$cds.coordinate,
                            cancer_abbr = rep("STAD_v" , nrow(mutation_list)),
                            geneId = NA )
library(org.Hs.eg.db) # map gene name to id by org.Hs.eg.db
e2s = toTable(org.Hs.egSYMBOL)
gcv_mut$geneId <- e2s$gene_id[match(gcv_mut$gene,e2s$symbol)]
sum(is.na(gcv_mut$geneId))
rm(e2s)


i <- sapply(gcv_mut, is.factor)
gcv_mut[i] <- lapply(gcv_mut[i], as.character)
rm(mutation_list)
dim(gcv_mut)

gcv_mut$classifed <- T
gcv_mut$classifed[-which(gcv_mut$geneId%in% tx37$GENEID)] <- F
assess <- data.frame(gene=unique(as.numeric(as.character(gcv_mut$geneId))) ) # extract all genes
head(assess )
assess$assess <- NA # 17280
assess <- assess[!is.na(assess$gene),]
for(i in 1:nrow(assess)) # see if the gene has correct CDS annotation: start with start codon, end to stop codon
{
    if(assess$gene[i] %in% tx37$GENEID)
        assess$assess[i] <-  assess.annotation(assess$gene[i])
}
    
table(assess$assess)

gcv_mut$classifed[-which(gcv_mut$geneId %in% assess$gene[assess$assess==0])] <- F
table(gcv_mut$classifed)

nmd.class.results <- matrix(ncol=8, nrow = nrow(gcv_mut))
colnames(nmd.class.results) <- names(classify.nmd(gcv_mut[1,]))

for(i in 1:nrow(gcv_mut))
{
    if(gcv_mut$classifed[i] == T) # it cana be classified
    {
        if( (gcv_mut$geneId[i] %in% tx37$GENEID))
            nmd.class.results[i,] <- classify.nmd(gcv_mut[i,])
        else nmd.class.results[,"note"] <- "can't find geneID in txdb"
    }

}
gcv_mut <- cbind(gcv_mut,nmd.class.results[,-3])
gcv_mut$msi <- NA
gcv_mut$msi [gcv_clinic$Microsatellite.instability.statusb[match(gcv_mut$sample,gcv_clinic$Sample.code)]==0] <- "MSS"
    gcv_mut$msi[gcv_clinic$Microsatellite.instability.statusb[match(gcv_mut$sample,gcv_clinic$Sample.code)]==1] <- "MSI"
table(gcv_mut$msi)

tb <- table(mutation_list$Sample.code)
gcv_clinic$n_mut <- as.numeric(tb[match(gcv_clinic$Sample.code,names(tb))])
# (gcv_mut$gene == "LARP4B")

gcv_mut$hypermutated <- F
gcv_mut$hypermutated[gcv_clinic$n_mut[match(gcv_mut$sample ,gcv_clinic$Sample.code)]>=1000] <- T

gcv_mut$n_mut <- gcv_clinic$n_mut[match(gcv_mut$sample ,gcv_clinic$Sample.code)]

colnames(gcv_mut)[14] <- "nmd.detail2"


gcv_mut$new.annotation <- NA
gcv_mut$new.annotation[gcv_mut$nmd.detail2==T] <- "NMD-elicit"
gcv_mut$new.annotation[gcv_mut$nmd.detail2==F&gcv_mut$have.ptc==T] <- "NMD-escape"
gcv_mut$new.annotation[gcv_mut$have.ptc==F&gcv_mut$effect%in%c("S")] <- "non-PTC SNP"
gcv_mut$new.annotation[gcv_mut$have.ptc==F&gcv_mut$effect%in%c("I","D")]  <- "non-PTC indels"
table(gcv_mut$new.annotation)

write.table(gcv_mut,file = "data/table/Gastric_validation_corhort_NMD_classifyresults.txt",sep="\t",row.names = F)

# heatmap----


my_keys <- keys(org.Hs.eg.db)
gene2go <-select(org.Hs.eg.db,
                 keys = as.character(my_keys),
                 columns=c("ENTREZID","SYMBOL","GENENAME","GO"),
                 keytype="ENTREZID")
dim(gene2go)
interested_go <- c("GO:0016568","GO:0044822","GO:0006281")

interested_go_term <- select(GO.db, keys = as.character(interested_go ),columns=c("GOID", "TERM"),keytype="GOID")[,2]
go_mut_tb <-  data.frame(sample =  unique(gcv_mut$sample))

for(i in 1:length(interested_go))
{
    go_mut_tb <- cbind(go_mut_tb,rep(F,nrow(go_mut_tb)))
    mmr_gene <- gene2go$SYMBOL[which(gene2go$GO==interested_go[i])] #genes in the GO
    go_mut_tb[go_mut_tb$sample %in% gcv_mut$sample[gcv_mut$nmd.detail2==T&gcv_mut$gene%in%mmr_gene],i+1] <-T
}
colnames(go_mut_tb) <- c("sample",interested_go_term )
head(go_mut_tb)
# get the dataframe for plotting-----
go_mut_tb2 <- data.frame(sample=rep(go_mut_tb$sample,length(interested_go)),
                         GO = rep(interested_go_term,each=nrow(go_mut_tb)),
                         NMD_mutated =as.vector(as.matrix(go_mut_tb[,-1])))
go_mut_tb2 <- rbind(go_mut_tb2 ,data.frame(sample=go_mut_tb$sample,
                                           GO = rep("microsatellite stability",length(go_mut_tb$sample)),
                                           NMD_mutated =gcv_mut$msi[match(go_mut_tb$sample,gcv_mut$sample)]))
go_mut_tb2$cancer <- gcv_mut$cancer_abbr[match(go_mut_tb2$sample,gcv_mut$sample)]

# add the mutation frequency
go_mut_tb2 <- rbind(go_mut_tb2 ,data.frame(sample=go_mut_tb$sample,
                                           GO = rep("log10 mutation freq",length(go_mut_tb$sample)),
                                           NMD_mutated =as.character(round(log10(gcv_clinic$n_mut))[match(go_mut_tb$sample,gcv_clinic$Sample.code)]),
                                           cancer = rep("STAD_v",length(go_mut_tb$sample))))

gcv_mut$effect2 <- NA


# add the mutation of TP53, PTEN, POLE, POLD1------
interested_genes <- c("LARP4B","EIF5B","PTEN")
for(i in 1:length(interested_genes ))
{
    go_mut_tb2_add <- data.frame(sample=go_mut_tb$sample,
                                 GO = rep(interested_genes[i],length(go_mut_tb$sample)),
                                 NMD_mutated =rep(NA,length(go_mut_tb$sample)),
                                 cancer = rep("STAD_v",length(go_mut_tb$sample)))
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% gcv_mut$sample[gcv_mut$have.ptc==F&gcv_mut$gene==interested_genes[i]]] <- "nonsilent"
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% gcv_mut$sample[gcv_mut$nmd.detail2==F&
                                                                        gcv_mut$effect%in%c("D","I")&
                                                                        (!as.character(nchar(gcv_mut$reference))%in%as.character((1:100)*3))&
                                                                        (!as.character(nchar(gcv_mut$alt))%in%as.character((1:100)*3))&
                                                                        gcv_mut$gene==interested_genes[i]]] <- "frameshift"
    go_mut_tb2_add$NMD_mutated[go_mut_tb$sample %in% gcv_mut$sample[gcv_mut$nmd.detail2==T&gcv_mut$gene==interested_genes[i]]] <- "NMD-elicit"
    

    go_mut_tb2 <- rbind(go_mut_tb2 ,
                        go_mut_tb2_add)
    
}



table(go_mut_tb2$NMD_mutated)

# rearrannge the position of x/y axis-----

# table(go_mut_tb2$GO)

go_mut_tb2$NMD_mutated <- factor(go_mut_tb2$NMD_mutated,levels=c("1","2","3","MSI","MSS","FALSE","TRUE","nonsilent","NMD-elicit","NMD-escape"))

# go_mut_tb2$sample <- factor(go_mut_tb2$sample,levels = unique(go_mut_tb2$sample[go_mut_tb2$GO=="log10 mutation freq"][order(gcv_clinic$mut[match(go_mut_tb2$sample[go_mut_tb2$GO=="log10 mutation freq"],gcv_clinic$sample)],decreasing = T)]))
go_mut_tb2$sample <- factor(go_mut_tb2$sample,levels = unique(go_mut_tb2$sample[go_mut_tb2$GO=="microsatellite stability"][order(paste(go_mut_tb2$NMD_mutated[go_mut_tb2$GO=="microsatellite stability"],(10000000-gcv_clinic$mut[match(go_mut_tb2$sample[go_mut_tb2$GO=="microsatellite stability"],gcv_clinic$sample)])),decreasing = F)]))


go_mut_tb2$GO <- factor(go_mut_tb2$GO,levels=c(interested_go_term,interested_genes,
                                               "microsatellite stability","log10 mutation freq","hypermutation"))


# new more dense firgue----
go_mut_tb3 <- go_mut_tb2
samplesforplot=unique(go_mut_tb3$sample)
samplesforplot <- samplesforplot[order(paste((100000000-gcv_clinic$n_mut[match(samplesforplot,gcv_clinic$Sample.code)]),
                                             gcv_clinic$Microsatellite.instability.statusb[match(samplesforplot,gcv_clinic$Sample.code)]))]

go_mut_tb3$sample <- factor(go_mut_tb3$sample,
                            levels = samplesforplot)


# choose color
levels(go_mut_tb3$NMD_mutated)
# [1] "1"          "2"          "3"          "MSI"        "MSS"        "FALSE"     
# [7] "TRUE"       "nonsilent"  "NMD-elicit" "NMD-escape"     


tiff("figures/Figure.tiff", width =9, height = 3.2, units = 'in', res = 500)
ggplot(go_mut_tb3,aes(x=sample,y=GO,fill=NMD_mutated))+geom_tile(col="white")+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    theme(legend.position="bottom")+
    labs(title = "")+scale_fill_manual(values=color_heatmap[c(2:4,6:9,11:12)])+
    guides(fill = guide_legend(nrow = 5,title =""))
dev.off()

# the barplot----

go_mut_sample <- data.frame(sample=levels(go_mut_tb3$sample))
go_mut_sample <- as.data.frame(table(gcv_mut[gcv_mut$sample %in% go_mut_sample$sample,
                                                    c("effect","sample")]))
head(go_mut_sample )


go_mut_sample$sample <- factor(go_mut_sample$sample,levels=samplesforplot)

tiff("figures/FigureS4-2.tiff", width =9, height = 3.5, units = 'in', res = 500)

ggplot(go_mut_sample[go_mut_sample$sample%in% samplesforplot,],aes(x=sample,y=Freq,fill=effect))+
    geom_bar(stat="identity",position = "stack")+labs(x="",fill="")+scale_fill_discrete(labels=c("Del","Ins","SN"))+
    theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),legend.position="top")

dev.off()

