#Preprocess-------- -------------------
preprocess.hiseq <- function(seq)
{
  rownames(seq) <- seq[,1]
  seq <- seq[,-1]
  colnames(seq) <- gsub(colnames(seq),pattern = "\\.",replacement = "-")
  return(seq)
}

#--------------------------------------------------------------
#detect the location of Open Reading Frame---------------------
#http://www.talkstats.com/showthread.php/28335-Finding-start-and-stop-codons-in-a-DNA-sequence
library("Biostrings")

# The START and STOP codons are stored in the variables start_codons and stop_codons
start_codons <- c("ATA","ATG")
stop_codons <- c("TGA","TAA","TAG")

# We use the function matchPattern to find the START codons in the sequence and in the inverted sequence
# Their positions are stored in start_pos and revstart_pos and are sorted

# find the position of *start* codon in +strand
findStartPos <- function(start_codons=c("ATA","ATG"), seq)
{
  start_pos <- c()
  #revstart_pos <- c()
  for (codon in start_codons)
  {
    matches <- matchPattern(codon,seq)
    start_pos <- c(start_pos, start(matches))
    #revmatches <- matchPattern(reverseComplement(DNAString(codon)),seq)
    #revstart_pos <- c(revstart_pos, start(revmatches))
  }
  start_pos <- sort(start_pos)
  #revstart_pos <- sort(revstart_pos, decreasing = TRUE)
  return(start_pos)#return(list(start_pos,revstart_pos))
}

# find the position of *stop* codon in +strand
findStopPos <- function(stop_codons = c("TGA","TAA","TAG"),seq)
{
  stop_pos <- c()
  #revstop_pos <- c()
  for (codon in stop_codons)
  {
    matches <- matchPattern(codon,seq)
    stop_pos <- c(stop_pos, start(matches))
    #revmatches <- matchPattern(reverseComplement(DNAString(codon)),seq)
    #revstop_pos <- c(revstop_pos, start(revmatches))
  }
  stop_pos <- sort(stop_pos)
  #revstop_pos <- sort(revstop_pos, decreasing = TRUE)
  return(stop_pos)#return(list(stop_pos,revstop_pos))
}

# find the ORF in '+' strand
findORF <- function(start_pos, k, stop_pos)
{  # k : Minimum size of Open Reading Frame

  stop_pointers<- c(0,0,0) # frame 1,2,3
  count <- 0
  ORF <- list()
  for (current_start in start_pos) {
    frame <- (current_start%%3) + 1 # decide which frame by position of start codon
    stop_pointer <- stop_pointers[frame] # 0 at first
    if (stop_pointer <= length(stop_pos) && # there are other stop positions
        (stop_pointer == 0 || stop_pos[stop_pointer] < current_start)) # current stop position is 0 or < start position
    {
      stop_pointer <- stop_pointer + 1 # pick up a larger stop postition
      while ( (stop_pointer <= length(stop_pos))
              && ((stop_pos[stop_pointer] <= current_start)
                  || (((stop_pos[stop_pointer]%%3) + 1) != frame)) ) # not on the same frame with start position
      {
        stop_pointer <- stop_pointer + 1 # try next stop position
      }
      stop_pointers[frame] <- stop_pointer # get the proper stop position

      if (stop_pointer <= length(stop_pos))
      {
        if ((stop_pos[stop_pointer] + 2 - current_start + 1) > k )
        {
          count <- count + 1
          new_ORF <- c(n=count,Frame=frame,Start=current_start,Stop=stop_pos[stop_pointer],Length=stop_pos[stop_pointer]+2-current_start+1)
          # Sequence=subseq(seq,current_start,stop_pos[stop_pointer]+2)
          ORF[[length(ORF)+1]] <- new_ORF

        }
      }
    }
  }
  ORF <- t(as.data.frame(ORF, as.is=T))
  rownames(ORF) <- NULL
  return(ORF)
}

# find the ORF in '-' strand
findrevORF <- function(revstart_pos, k, revstop_pos)
{
  revstop_pointers <- c(0,0,0)
  count <- 0
  ORF <- list()
  for (current_revstart in revstart_pos)
  {
    current_revstart <- current_revstart + 2 # start position on the reverse complement sequence
    frame <- (current_revstart%%3) + 1
    revstop_pointer <- revstop_pointers[frame]
    if (revstop_pointer <= length(revstop_pos) && (revstop_pointer == 0 || revstop_pos[revstop_pointer] > current_revstart))
    {
      revstop_pointer <- revstop_pointer + 1
      while ( (revstop_pointer <= length(revstop_pos))
              && ((revstop_pos[revstop_pointer] + 2 >= current_revstart)
                  || ((((revstop_pos[revstop_pointer] + 2)%%3) + 1) != frame)) )
      {
        revstop_pointer <- revstop_pointer + 1
      }
      revstop_pointers[frame] <- revstop_pointer

      if (revstop_pointer <= length(revstop_pos))
      {
        if ((current_revstart - revstop_pos[revstop_pointer])+1 > k)
        {
          count <- count + 1
          new_ORF <- c(n=count,Frame=-frame,Start=current_revstart,Stop=revstop_pos[revstop_pointer]+2,Length=current_revstart-revstop_pos[revstop_pointer])
          ORF[[length(ORF)+1]] <- new_ORF
        }
      }
    }
  }
  ORF <- t(as.data.frame(ORF, as.is=T))
  rownames(ORF) <- NULL
  return(ORF)
}

#-----------------------------------------
# NMD classification functions------------
# retrieve the CDS sequence of a transcript & position of exon-exon junctions
get_tx_seq <- function(tx_pos, ncbi_build)
{
  pos<- numeric() # position on chromosome
  exon_jun <- numeric() # junction positions on chromosome
  if(ncbi_build==36)
  {
    genome <- genome36
  }
  if(ncbi_build==37)
  {
    genome <- genome37
  }
  tx_pos <- tx_pos[!is.na(tx_pos$CDSSTART),] # some genes have no CDS
  if(nrow(tx_pos)>0)
  {
    for(j in 1:nrow(tx_pos))
    {
      # seq <- paste(seq,as.character(getSeq(genome, tx_pos$EXONCHROM[j], tx_pos$CDSSTART[j], tx_pos$CDSEND[j])),sep = "")
      start = tx_pos$CDSSTART[j] #
      end = tx_pos$CDSEND[j]
      start= as.numeric(start)
      end=as.numeric(end)
      pos<- c(pos,start:end)
      #exon_jun <- c(exon_jun,1,rep(0,end-start-1),1)
    }

    #seq <- as.character(seq)
    # +/- strand on which the transcript is
    strand <- unique(tx_pos$EXONSTRAND)
    if(length(strand) > 1)
      print(paste(tx_pos$TXID ,"Wrong strand!"))
    if(strand == "+")
    {
      pos <- sort(pos,decreasing = F)
      exon_jun <- rep(0,length(pos))
      exon_jun[(pos[-1]-pos[-length(pos)])>1] <- 1
      exon_jun[which(exon_jun==1)+1] <- 1
    }

    if(strand == "-")
    {
      pos <- sort(pos,decreasing = T)
      exon_jun <- rep(0,length(pos))
      exon_jun[(pos[-length(pos)]-pos[-1])>1] <- 1
      exon_jun[which(exon_jun==1)+1] <-1
    }

    return(list(wt_pos = pos, exonJunction = exon_jun))


    #return(list(seq = seq,  wt_pos = pos, exonJunction = exon_jun,cds = c(tx_pos$CDSSTART[1],tx_pos$CDSEND[nrow(tx_pos)])))
  }
  # no cds will return NULL
}

# Get the mutated exon exon junction
get_mut_exjun <- function(wt_exjun,mut_start,mut_end, mut_nt)
{
  note <-""
  if(any(wt_exjun[mut_start:mut_end]==1)) # mutation on the exon-intron junction
  {
    note <- "mutated exon intron junction"
    mut_exjun <-  c(wt_exjun[1:(mut_start-1)],rep(-1,length(mut_nt)), wt_exjun[(mut_end+1):length(wt_exjun)])
  }else if(mut_start==1)
  {
    mut_exjun <- c(rep(-1,length(mut_nt)), wt_exjun[(mut_end+1):length(wt_exjun)])
    note <- "mutated start"
  }else if(mut_end==length(wt_exjun))
  {
    mut_exjun <- c(wt_exjun[1:(mut_start-1)],rep(-1,length(mut_nt)))
    note <- "mutated end"
  }else if(all(mut_start:mut_end %in% 1:length(wt_exjun))& all(wt_exjun[mut_start:mut_end]!=1)) # mutation is not at the start or the end
  {
    mut_exjun <- c(wt_exjun[1:(mut_start-1)],rep(-1,length(mut_nt)), wt_exjun[(mut_end+1):length(wt_exjun)])
  }  else
  {
    mut_exjun <-""
    note <- "unknown mutation"
  }
  return(list(mut_exjun,note))
}

# Get the mutated sequence
get_mut_seq_exjun <- function(wt_exjun,wt_seq,mut_start,mut_end, mut_nt,ref_nt,note,txid)
{
    mut_seq = NA
    mut_exjun = NA
  if(mut_nt == "-") # deletion
  {
    if(ref_nt== substr(wt_seq,mut_start,mut_end))
    {
      mut_seq <- paste(c(substr(wt_seq,1,mut_start-1),substr(wt_seq,mut_end+1,nchar(wt_seq))),collapse= "")
    }
    else note<- "reference unmatched"

    if(mut_start>1 & mut_end < length(wt_exjun)) # mutation is not at the start or the end
      mut_exjun <- c(wt_exjun[1:(mut_start-1)], wt_exjun[(mut_end+1):length(wt_exjun)])
    else if(mut_start==1)
      mut_exjun <- c(wt_exjun[(mut_end+1):length(wt_exjun)])
    else if(mut_end==length(wt_exjun))
      mut_exjun <- c(wt_exjun[1:(mut_start-1)])
  }
  if(ref_nt == "-") # insertion
  {
    mut_seq <- paste(c(substr(wt_seq,1,mut_start),mut_nt,substr(wt_seq,mut_end,nchar(wt_seq))),collapse= "")
    if(all(mut_start:mut_end %in% 1:length(wt_exjun)))
      mut_exjun <- get_mut_exjun(wt_exjun = wt_exjun,mut_nt = mut_nt,mut_start = mut_start,mut_end = mut_end)
    mut_exjun <-mut_exjun[[1]]
  }

  if(ref_nt != "-" & mut_nt != "-") # deletion and insertion in same site
  {
    if(ref_nt== substr(wt_seq,mut_start,mut_end))
    {
        mut_seq <- paste(c(substr(wt_seq,1,mut_start-1),mut_nt,substr(wt_seq,mut_end+1,nchar(wt_seq))),collapse= "")
        mut_exjun <- get_mut_exjun(wt_exjun = wt_exjun,mut_nt = mut_nt,mut_start = mut_start,mut_end = mut_end)
        if(mut_exjun[[2]]!="") note <- paste(note,mut_exjun[[2]],sep = " ")
        mut_exjun <-mut_exjun[[1]]
    }

    else
    {
        note <- "reference unmatched"
    }


  }
  return(list(seq = mut_seq,exjun = mut_exjun,note= note))
}

# get TPM of genes for all mutations
get.expression.vector <- function(mut)
{
  mut$expression <- rep(NA, nrow(mut))

  for(i in 1:nrow(mut))
  {
    mut$expression[i] <-  as.numeric(as.character(meta_tpm[rownames(meta_tpm)==mut$Entrez_Gene_Id[i],
                                                           colnames(meta_tpm)==mut$sample[i] ] ))
  }
  return(mut)
}

# get gene annotation
get.anno <- function(mut_entry)
{
  gene_id <-unlist(mut_entry["geneId"])

  # get the transcripts of the genes

    tx <- tx37[tx37$GENEID==gene_id,]
    cds <- cds37[cds37$geneid==gene_id,]

  return(list(tx=tx,cds = cds))
}


get.anno.by.geneID <- function(gene_id)
{
  # get the transcripts of the genes

  tx <- tx37[tx37$GENEID==gene_id,]
  cds <- cds37[cds37$geneid==gene_id,]

  return(list(tx=tx,cds = cds))
}
# get gene annotation from old build ref
get.anno36 <- function(mut_entry)
{
  gene_id <-unlist(mut_entry[13])

  # get the transcripts of the genes
  tx <- tx36[tx36$GENEID==gene_id,]
  cds <- cds36[cds36$geneid==gene_id,]

  return(list(tx=tx,cds = cds))
}

# note: remove splice site mutation first
# get the longest cds
get.longest.cds <- function(mut_entry,ref) # combine
{
  gene_id <-as.character (mut_entry["geneId"])

  #get gene annotations
  if(ref==37)
  {
    anno <- get.anno(mut_entry)
    tx <- anno$tx
    cds <- anno$cds
  }
  if(ref==36)
  {
    anno <- get.anno36(mut_entry)
    tx <- anno$tx
    cds <- anno$cds
  }

  cds <- cds[which.max(cds$width),]
  tx <- tx[tx$TXID == cds$txid,]

  return(list(gene = as.character(cds$geneid), txid = as.character(cds$txid) , seq= as.character(cds$seq), tx = tx))
}

get.longest.cds.by.geneID <- function(gene_id) # combine
{

  anno <- get.anno.by.geneID(gene_id)
  tx <- anno$tx
  cds <- anno$cds
  cds <- cds[which.max(cds$width),]
  tx <- tx[tx$TXID == cds$txid,]

  return(list(gene = as.character(cds$geneid), txid = as.character(cds$txid) , seq= as.character(cds$seq), tx = tx))
}

# get the mutated sequence
get.mut.seq <- function(mut_entry)
{
  cds <- get.longest.cds(mut_entry = mut_entry,ref=37)
  ref <- 37
  tx <- cds$tx
  txid <- cds$txid
  wt_seq <- cds$seq
  n.exon <- nrow(tx[!is.na(tx$CDSSTART),])
  wt_pos <- get_tx_seq(tx,ncbi_build = ref)
  wt_exjun <- wt_pos$exonJunction
  wt_pos <- wt_pos$wt_pos
  strand <- as.character (unique(tx$EXONSTRAND))
  note<- ""
  mut_seq <- NA
  mut_exjun <- NA
    if(length(cds$gene) >0)
  {


    if(strand == "+") # gene on the + strand
    {
      ref_nt <- as.character(mut_entry[5])
      mut_nt <- as.character(mut_entry[6])
      mut_start <- as.numeric(mut_entry[3])
      mut_end <- as.numeric(mut_entry[4])
    }
    if( strand  == "-")
    {
      ref_nt <- as.character(mut_entry[5])
      mut_nt <- as.character(mut_entry[6])
      mut_start <- as.numeric(mut_entry[3])
      mut_end <-  as.numeric(mut_entry[4])
      ref_nt <- as.character(reverseComplement(DNAString(ref_nt)))
      mut_nt <- as.character(reverseComplement(DNAString(mut_nt)))

    }
    if(!any(mut_start:mut_end %in% min(wt_pos):max(wt_pos))) # mutation out of the transcript
    {


      note <- "mutation outside"

      mut_seq <- wt_seq
      mut_exjun <- wt_exjun

    } else if( all(mut_start:mut_end %in% min(wt_pos):max(wt_pos)) & (!any(mut_start:mut_end %in% wt_pos))) # mutation position in intron
    {
      note <- "mutation in intron"
      mut_seq <- wt_seq
      mut_exjun <- wt_exjun
    }else if( all(mut_start:mut_end %in% wt_pos)) # mutation in CDS
    {
      mut_start <- which(wt_pos==mut_start)
      mut_end <- which(wt_pos==mut_end)
      note<-""
    } else # part of the mutation in exon
    {
      note  <-  "mutation partially in cds"
      mut_seq <- NA
      mut_exjun <- NA
    }


    if(note=="" ) # check if mutation is on splice site
    {
      mut_start <- c(mut_start,mut_end)
      mut_end <- max(mut_start)
      mut_start <- min(mut_start)
      # use function
      mut_exjun<- get_mut_seq_exjun(wt_exjun = wt_exjun,wt_seq=wt_seq,mut_nt = mut_nt,mut_start = mut_start,
                                    mut_end = mut_end, ref_nt = ref_nt,txid = txid, note =note)
      note <- mut_exjun$note
      mut_seq <- mut_exjun$seq
      mut_exjun <- mut_exjun$exjun
    }else{
      mut_seq <- NA
      mut_exjun <- NA
    }

    return(list(mut_seq= mut_seq, mut_exjun = mut_exjun, note= note,wt_seq=wt_seq,n.exon=n.exon))
  }

}

# find where is the last exon-exon junction or decide if it is a single-exon gene
last.exexjun.pos <- function(mut_exjun)
{
  note <- ""
  if((sum(mut_exjun==1)) %in% ((1:1000)*2))
  {
    last_exon_exon_junction<- which(mut_exjun==1)[length( which(mut_exjun==1))-1]
    note <- ""
  }else if(sum(mut_exjun==1) %in% ((1:1000)*2+1))
  {
    last_exon_exon_junction<- which(mut_exjun==1)[length( which(mut_exjun==1))-1]
    note <- "splice site mutated"
  }else if(sum(mut_exjun==1)==1)
  {
    last_exon_exon_junction<- which(mut_exjun==1)
    note<- "splice site mut"
  }else if(sum(mut_exjun==1)==0)
  {
    last_exon_exon_junction<-length(mut_exjun)
    note <- "single exon"
  }
  return(list(last_exon_exon_junction=last_exon_exon_junction,note=note))
}

# classify
mut.info <- function(mut_entry)
{
  mut_seq <- get.mut.seq(mut_entry) # get the mutated seq, wt seq, pos of mut exon-exon junctions and number of exons
  note <- mut_seq$note
  wt_seq <- mut_seq$wt_seq
  n.exon <- mut_seq$n.exon
  if(length(note)>1)
    note <- paste(note,collapse = " ")
  note <- gsub("  ", "", note, fixed = TRUE)
  mut_exjun <- mut_seq$mut_exjun
  mut_seq <- mut_seq$mut_seq
  if(!is.na(mut_seq))
  {
    last.exexjun <- last.exexjun.pos(mut_exjun)
    if(last.exexjun$note!="")  note <- paste(note,last.exexjun$note,sep=",")
    last.exexjun <- last.exexjun$last_exon_exon_junction
  }else
  {
    last.exexjun = NA
  }
  return(c(mut_seq=mut_seq,last.exexjun=last.exexjun,note=note,wt_seq = wt_seq, n.exon = n.exon))
}

# classify the mutations into NMD or non-NMD, and return the position of premature stop codon in CDS
classify.nmd <- function(mut_entry)
{
  mut_info <- mut.info(mut_entry )
  mut_seq <-as.character( mut_info[1] )
  wt_seq <- as.character(mut_info[4])
  n.exon <- as.numeric(mut_info[5])
  last_exon_exon_junction <- as.numeric(mut_info[2])
  nmd.detail <- NA
  nmd.detail2 <- NA
  note1<-as.character(mut_info[3])
  note<-""
  pos.ptc <- NA
  have.ptc <- F
  length <- NA
  if(!is.na(mut_seq))
  {
    length <- nchar(mut_seq)
    if(nchar(mut_seq)>0) # it has mutated sequences
    {

      start_pos = findStartPos(seq=mut_seq) # find the start position
      stop_pos = findStopPos(seq=mut_seq)# find the stop position in the frame

      if(length(start_pos)==0) # cannot find start codon
      {
        note <- "no start codon"
      }else if(length(stop_pos)==0) # cannot find start codon
      {
        note <- "no stop codon"
      }else{  # find start and stop codon
        start_pos = min(start_pos) # choose the cloest as the start position

        ORF <- findORF(start_pos=start_pos,k = 0,stop_pos=stop_pos)

        if( nrow(ORF)==0)
        {
          nmd.detail  = F
          note = "cannot find ORF"
          pos.ptc <- Inf # no stop codon
        }else # find ORF
        {
          ORF <-  ORF[which.min( ORF[,"Stop"]),]
          pos.ptc <- ORF["Stop"]
          if(pos.ptc > nchar(mut_seq)-3) # normal coding sequence no PTC
          {
            have.ptc <- F
            note = "normal sc"
            nmd.detail <- F
          }else{
            if(pos.ptc <= nchar(mut_seq)-3) have.ptc <- T
            if(n.exon<1.1)
            {
              nmd.detail <- F
              note = "single exon"
            }else if(pos.ptc < (last_exon_exon_junction -50) & pos.ptc > 200){
              nmd.detail <- T
            }else if(pos.ptc >= (last_exon_exon_junction -50)& pos.ptc <(nchar(mut_seq)-2))
            {
              nmd.detail <- F
              note = "PTC too close to leej"
            }else if(pos.ptc < 200)
            {
              nmd.detail <- F
              note = "PTC too close to start codon"
            }else if(is.infinite(pos.ptc)){
              nmd.detail <- F
              note = "nonstop"
            }else{
              nmd.detail <- F
              note = "unclear"
            }
          } # have PTC

        }

      }

    }
    # find the new start and stop codons in the sequence
    if(nchar(wt_seq)>0) # it has wt sequences
    {
      start_pos = findStartPos(seq=wt_seq) # find the start position
      stop_pos = findStopPos(seq=wt_seq)# find the stop position in the frame

      if(length(start_pos)==0) # cannot find start codon
      {
        note <- paste(note,"no start codon in wt")
      } else if(length(stop_pos)==0) # cannot find start codon
      {
        note <-paste(note, "no stop codon in wt" )
      }else # find start and stop codon
      {
        start_pos = min(start_pos) # choose the cloest as the start position
        if(start_pos!=1) note <-paste(note, "wt cds not start with start codon" )
        ORF <- findORF(start_pos=start_pos,k = 0,stop_pos=stop_pos)

        if(nrow(ORF)==0)
        {
          nmd.detail2  = F
          note = paste(note, "can't find ORF in wt")
        }else{ # find ORF in WT
          ORF <-  ORF[which.min( ORF[,"Start"]),]

          if(ORF["Stop"] ==nchar(wt_seq)){
            nmd.detail2 <- T
            note <-paste(note, "PTC in wt" )
          }else nmd.detail2 <- F

        }
      }

    }
  }
  note <- paste(note1,note)
  return(c( mut_nmd = nmd.detail,note=note, wt_nmd = nmd.detail2, PTC = pos.ptc,have.ptc=have.ptc,mutseq_length=length,last_exon_exon_junction=last_exon_exon_junction,n.exon=n.exon))
}

# if the number if whole number?
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

get.build <- function(mut_entry)
{
  # which ref build it use 36 or 37
  ref.no = NA
  cds <- get.longest.cds(mut_entry,ref = 36)
  tx <- cds$tx
  txid <- cds$txid
  wt_seq <- cds$seq
  n.exon <- nrow(tx[!is.na(tx$CDSSTART),])
  wt_pos <- get_tx_seq(tx,ncbi_build = 36)
  wt_exjun <- wt_pos$exonJunction
  wt_pos <- wt_pos$wt_pos
  strand <- as.character (unique(tx$EXONSTRAND))
  note<- ""
  mut_start <- as.numeric(mut_entry[3])
  mut_end <- as.numeric(mut_entry[4])

  if(!is.null(wt_pos))
  {
    if(any(mut_start:mut_end %in% wt_pos))
    {
      ref.no = 36
    }
  }
  cds <- get.longest.cds(mut_entry,ref = 37)
  tx <- cds$tx
  txid <- cds$txid
  n.exon <- nrow(tx[!is.na(tx$CDSSTART),])
  wt_pos <- get_tx_seq(tx,ncbi_build = 37)
  wt_pos <- wt_pos$wt_pos
  strand <- as.character (unique(tx$EXONSTRAND))

  if(strand == "+") # gene on the + strand
  {
    mut_start <- as.numeric(mut_entry[3])
    mut_end <- as.numeric(mut_entry[4])
  }
  if( strand  == "-")
  {
    mut_start <- as.numeric(mut_entry[3])
    mut_end <-  as.numeric(mut_entry[4])
  }
  if(any(mut_start:mut_end %in% wt_pos)&is.na(ref.no))
  {
    ref.no = 37
  }else if(any(mut_start:mut_end %in% wt_pos))
  {
    ref.no = paste(ref.no ,"37")
  }
  return(ref.no)
}
# see if the gene has correct CDS annotation
assess.annotation <- function(gene)
{
  cds <- get.longest.cds.by.geneID(gene_id = gene)
  # whether it is right coding sequence
  if(nchar(cds$seq)>0){
    if(is.wholenumber(nchar(cds$seq)/3))
    {
      wt_seq <- cds$seq
      start_pos = findStartPos(seq=wt_seq) # find the start position
      stop_pos = findStopPos(seq=wt_seq)# find the stop position in the frame
      if(length(start_pos)==0) # no start codon found
      {
        index=1 #"no start pos"
      }
      else
      {
        if(min(start_pos)!=1) # not start wiht start codon
        {
          index=2 # "not start with start codon"
        }
        else
        { # do start with start codon
          if(length(stop_pos)==0) # have no stop codon
          {
            index=3 #"have no stop codon"
          }
          else
          { # do have stop codon
            ORF <- findORF(start_pos=min(start_pos),k = 0,stop_pos=stop_pos)
            if(nrow(ORF)==0)
            {
              index=4 # "no ORF found"
            }
            else
            {
              if(min(ORF[,4])+2 < nchar(wt_seq))
              {
                index=5 # "PTC in WT"
              }
              else
              {
                index=0 # "right annotation"
              }

            }
          }
        }
      }


    }else
    {
      index=6#return("length not 3 multiple")
    }

  }else index=7#return("no cds")

  return(index)
}

#-------------------------------------------------------
# Calculate the features of mutations ------------------
#-------------------------------------------------------
# We will calculate:
# downstreamexons
# donordist
# utr3diststart
# X50bp
# stopdist
# relCDSpos
# utr3distend startdist
# nexons utr3size acceptordist trnscaffected utr5distend af
# utr5size onecodingexon relcDNApos utr5diststart EncExp
# GerpN mamPhyloP lofflag priPhyloP minDistTSS GerpS verPhyloP


feature.mut <- function(mut_entry) # get features of mutations based on the mut pos, gene id and NCBI build
{
  features <- c("X50bp","startdist","stopdist","utr3size","utr5size",
                "trnscaffected","ntrnsc","nexons","donordist","acceptordist","onecodingexon",
                "relcDNApos","relCDSPos","downstreamexons")
  mut_feature <- rep(NA,length(features))
  names(mut_feature) <- features

  ptc.pos <- mut_entry["pos_ptc"]

  gene_id <-as.numeric( mut_entry["Entrez_Gene_Id"]) # gene id
  mut_start <- as.numeric(mut_entry["Start_Position"]) # mut start position
  mut_end <- as.numeric(mut_entry["End_Position"]) # mut end position
  #get gene annotations
  anno <- get.anno(mut_entry[c(2,4)])
  tx <- anno$tx
  cds <- anno$cds

  if(mut_entry["nmd.detail2"]==T) mut_feature["X50bp"] <- 1
  if(mut_entry["nmd.detail2"]==F) mut_feature["X50bp"] <- 0
  mut_feature["ntrnsc"] <- length(unique(tx$TXID))-1 # number of alternative isoforms
  mut_feature["trnscaffected"] <-length(unique(tx$TXID[ mut_start > tx$TXSTART & mut_start < tx$TXEND ])) # how many isoforms are affected

  longest.cds <- get.longest.cds(mut_entry = mut_entry)
  ref <- as.character(mut_entry[4])
  tx <- longest.cds$tx
  txid <- longest.cds$txid
  wt_seq <- longest.cds$seq
  n.exon <- nrow(tx[!is.na(tx$CDSSTART),])
  mut_feature["nexons"] <- n.exon # number of exons
  if(n.exon == 1) mut_feature["onecodingexon"] <- 1
  if(n.exon > 1) mut_feature["onecodingexon"] <- 0
  wt_pos <- get_tx_seq(tx,ncbi_build = ref)
  wt_exjun <- wt_pos$exonJunction
  wt_pos <- wt_pos$wt_pos
  strand <- as.character (unique(tx$EXONSTRAND))
  note<- ""
  mut_feature["startdist"] <- min(mut_start,mut_end) - min(wt_pos)# distance to start codon
  mut_feature["stopdist"] <-  max(wt_pos) - max(mut_start,mut_end)# distance to stop codon
  if(grepl(36,mut_entry["Ncbi_Build"]))
  {
    mut_feature["utr3size"] <- hg18_txlens$utr3_len[hg18_txlens$tx_id==txid]
    mut_feature["utr5size"] <- hg18_txlens$utr5_len[hg18_txlens$tx_id==txid]
  }
  if(grepl(37,mut_entry["Ncbi_Build"]))
  {
    mut_feature["utr3size"] <- hg19_txlens$utr3_len[hg19_txlens$tx_id==txid]
    mut_feature["utr5size"] <- hg19_txlens$utr5_len[hg19_txlens$tx_id==txid]
  }

  mut_feature["relcDNApos"] <- (mut_start-min(tx$TXSTART[1],tx$TXEND[1]))/abs(tx$TXSTART[1]-tx$TXEND[1])
  mut_feature["relCDSPos"] <- match(mut_end,wt_pos)/length(wt_pos)
  if(n.exon >1)
  {
    donors <- which(wt_exjun==1)[(1:(length(which(wt_exjun==1))/2))*2-1]
    donors <- wt_pos[donors]
    acceptors <- which(wt_exjun==1)[(1:(length(which(wt_exjun==1))/2))*2]
    acceptors <- wt_pos[acceptors]

    mut_feature["donordist"] <- min(c(abs(donors- mut_start), abs(donors-mut_end)))
    mut_feature["acceptordist"] <-  min(c(abs(acceptors-mut_start), abs(acceptors-mut_end)))
    mut_feature["downstreamexons"]<- sum(acceptors >= max(mut_start,mut_end))
  }else{
    mut_feature["downstreamexons"]<- 0
  }




  return(mut_feature)
}

feature.ptc <- function(mut_entry) # get features of mutations based on the mut pos, gene id and NCBI build
{
  ptc_feature <- rep(NA,length(features))
  names(ptc_feature) <- features
  if(mut_entry["Variant_Classification"] %in% c("FSD","FSI","NM"))
  {
    ptc.pos <- mut_entry["pos_ptc"]

    gene_id <-as.numeric( mut_entry["Entrez_Gene_Id"]) # gene id
    mut_start <- as.numeric(mut_entry["Start_Position"]) # mut start position
    mut_end <- as.numeric(mut_entry["End_Position"]) # mut end position
    #get gene annotations
    anno <- get.anno(mut_entry[c(2,4)])
    tx <- anno$tx
    cds <- anno$cds

    if(mut_entry["nmd.detail2"]==T) ptc_feature["X50bp"] <- 1
    if(mut_entry["nmd.detail2"]==F) ptc_feature["X50bp"] <- 0
    ptc_feature["ntrnsc"] <- length(unique(tx$TXID))-1 # number of alternative isoforms
    ptc_feature["trnscaffected"] <-length(unique(tx$TXID[ mut_start > tx$TXSTART & mut_start < tx$TXEND ])) # how many isoforms are affected

    longest.cds <- get.longest.cds(mut_entry = mut_entry)
    ref <- as.character(mut_entry[4])
    tx <- longest.cds$tx
    txid <- longest.cds$txid
    wt_seq <- longest.cds$seq
    n.exon <- nrow(tx[!is.na(tx$CDSSTART),])
    ptc_feature["nexons"] <- n.exon # number of exons
    if(n.exon == 1) ptc_feature["onecodingexon"] <- 1
    if(n.exon > 1) ptc_feature["onecodingexon"] <- 0
    wt_pos <- get_tx_seq(tx,ncbi_build = ref)
    wt_exjun <- wt_pos$exonJunction
    wt_pos <- wt_pos$wt_pos
    strand <- as.character (unique(tx$EXONSTRAND))
    note<- ""
    ptc_feature["startdist"] <- min(mut_start,mut_end) - min(wt_pos)# distance to start codon
    ptc_feature["stopdist"] <-  max(wt_pos) - max(mut_start,mut_end)# distance to stop codon
    if(grepl(36,mut_entry["Ncbi_Build"]))
    {
      ptc_feature["utr3size"] <- hg18_txlens$utr3_len[hg18_txlens$tx_id==txid]
      ptc_feature["utr5size"] <- hg18_txlens$utr5_len[hg18_txlens$tx_id==txid]
    }
    if(grepl(37,mut_entry["Ncbi_Build"]))
    {
      ptc_feature["utr3size"] <- hg19_txlens$utr3_len[hg19_txlens$tx_id==txid]
      ptc_feature["utr5size"] <- hg19_txlens$utr5_len[hg19_txlens$tx_id==txid]
    }

    ptc_feature["relcDNApos"] <- (mut_start-min(tx$TXSTART[1],tx$TXEND[1]))/abs(tx$TXSTART[1]-tx$TXEND[1])
    ptc_feature["relCDSPos"] <- match(mut_end,wt_pos)/length(wt_pos)
    if(n.exon >1)
    {
      donors <- which(wt_exjun==1)[(1:(length(which(wt_exjun==1))/2))*2-1]
      donors <- wt_pos[donors]
      acceptors <- which(wt_exjun==1)[(1:(length(which(wt_exjun==1))/2))*2]
      acceptors <- wt_pos[acceptors]

      ptc_feature["donordist"] <- min(c(abs(donors- mut_start), abs(donors-mut_end)))
      ptc_feature["acceptordist"] <-  min(c(abs(acceptors-mut_start), abs(acceptors-mut_end)))
      ptc_feature["downstreamexons"]<- sum(acceptors >= max(mut_start,mut_end))
    }else{
      ptc_feature["downstreamexons"]<- 0
    }
  }

  return(ptc_feature)
}

get.wt.expression <- function(gene,cancer,nmf)
{
  if(!is.na(nmf))
  {
    mut_sample <- unique(pancan_mut$sample[pancan_mut$cancer_abbr==cancer&pancan_mut$nmfc==nmf])
  }else{
    mut_sample <- unique(pancan_mut$sample[pancan_mut$cancer_abbr==cancer])
  }
  wt_sample <- rownames(mut_tb[rownames(mut_tb) %in% mut_sample,])[mut_tb[rownames(mut_tb) %in% mut_sample,colnames(mut_tb)==gene]==0]
  if(any(rownames(pancan_cnv)==gene) & sum(wt_sample %in% colnames(pancan_cnv))>10)
    wt_sample <- wt_sample[wt_sample %in% colnames(pancan_cnv)[pancan_cnv[rownames(pancan_cnv)==gene,] ==0]]

  expression_wt <- pancan_seq[rownames(pancan_seq)==gene, colnames(pancan_seq) %in% wt_sample]
  return(expression_wt)
}
