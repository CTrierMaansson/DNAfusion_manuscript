x #name of BAM file
mate_discovery_cigar <- function(x){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(bamsignals)
  library(ggseqlogo)
  library(ggplot2)
  `%ni%` <- Negate(`%in%`)
  start_time <- Sys.time()
  what <- c("mpos", "pos", "seq","cigar")
  which <- GRanges(seqnames="chr2", IRanges(start = 42169461, end = 42330874)) #EML4 gene location
  param <- ScanBamParam(which = which, what = what)
  bam <- scanBam(x, param = param)
  reads <- data.frame(sequences = bam$`chr2:42169461-42330874`$seq,
                      mate = bam$`chr2:42169461-42330874`$mpos,
                      position = bam$`chr2:42169461-42330874`$pos,
                      cigar = bam$`chr2:42169461-42330874`$cigar)
  reads <- reads %>% filter(mate < 29921586 & mate > 29192774) #ALKpostion
  if (length(reads$mate)<2){ #More than one read in EML4 with mate in ALK is required 
    res <- "No ALK was detected"
    return(res)
  }
  else{
    res <- reads
    clip_reads <- reads %>% filter(cigar != "96M") #Reads which align perfectly are excluded because they cannot discover breakpoint
    clip_reads <- clip_reads %>% filter(!grepl("D",cigar)) #Reads with artificial indels are removed
    clip_reads <- clip_reads %>% filter(!grepl("I",cigar))
  }
  fun <- function(str) sub("\\M.*", "",str) #Index of matching reads
  index1 <- sapply(clip_reads$cigar,FUN = fun)
  fun1 <- function(ind1){ #Removal of reads short clipped in the beginning of read, meaning the first bases in the read is matching ALK and not EML4
    if (length(ind1)>1){
      return(NA) #These reads have their S bases align to ALK and their remaining M bases matching EML4 AFTER the breakpoint
    }
    else{
      return(ind1)
    }
  }
  fun <- function(ind){
    splits <- strsplit(ind,split = "S")
    splits <- lapply(splits, as.numeric)
    splits_ind <- lapply(splits,FUN = fun1)
    return(splits_ind)
  }
  index2 <- sapply(index1,FUN = fun)
  clip_reads$indeces <- index2 #This number represents the bp number in the read where the break is
  clip_reads <- clip_reads %>% filter(!is.na(indeces))
  EML4_fun <- function(inp){ #Finding the sequence of EML4 based on the index of the breakpoint and the bases before the breakpoint
    return(substring(inp$sequences,(inp$indeces-19),inp$indeces))
  }
  ALK_fun <- function(inp){ #Finding the sequence of ALK based on the index of the breakpoint and the bases following the breakpoint
    return(substring(inp$sequences,(inp$indeces+1),(inp$indeces+20)))
  }
  char_df <- clip_reads %>% dplyr::select(sequences,indeces)
  EML4_seq <- apply(char_df, FUN = EML4_fun, MARGIN = 1) #Seq at 3' EML4 break
  EML4_tab <- table(EML4_seq) 
  ALK_seq <- apply(char_df, FUN = ALK_fun, MARGIN = 1) #Seq at 5' ALK break
  ALK_tab <- table(ALK_seq)
  clip_reads$indeces <- as.numeric(clip_reads$indeces)
  break_pos <- clip_reads$position + clip_reads$indeces-1 #bp number of last EML4 base
  break_pos_tab <- table(break_pos)
  uniqe_seq <- unique(clip_reads$sequences) #Identifying unique reads
  stop_pos <- as.numeric(names(which.max(break_pos_tab))) #Finding the position with the most hits -> most likely breakpoint
  depth <- bamCoverage(x,GRanges(
    seqnames = "chr2",IRanges(start=(stop_pos),end=stop_pos+1)),
    mapqual=0,verbose=F) #Read depth at breakpoint
  ALK20 <- ALK_seq[nchar(ALK_seq)>19] #Identified first sequence in ALK, only sequences with length 20 are analyzed
  EML420 <- EML4_seq[nchar(EML4_seq)>19] #Identified Last sequence in EML4, only sequences with length 20 are analyzed
  EML4_gg <- ggseqlogo(EML420, method = "prob")+ #Consensus plot EML4
    ggtitle("Last EML4 sequence")
  ALK_gg <- ggseqlogo(ALK20, method = "prob")+ #Consensus plot ALK
    ggtitle("First ALK sequence")
  end_time <- Sys.time()
  times <- end_time - start_time
  return(list(clipped_reads = clip_reads,
              last_EML4 = EML4_tab,
              breakpoint = break_pos_tab,
              first_ALK = ALK_tab,
              unique_reads = uniqe_seq,
              read_depth = max(depth[1]),
              runtime = times,
              consensus = gridExtra::grid.arrange(EML4_gg, ALK_gg)))
}