
library(Rsamtools)
library(bamsignals)
x #name on .txt file containing each gene, name of chromosome, start location, end location and strand
gr_ALK <- function(x){
  library(GenomicRanges)
  y <- read.table(x, header = T)
  y$chromStart <- as.numeric(y$chromStart) 
  y$chromEnd <- as.numeric(y$chromEnd) 
  granges <- GRanges(y$chrom, 
                     IRanges(start = y$chromStart, end = y$chromEnd),
                     strand = y$strand)
  mcols(granges)$SYMBOL <- y$name
  mcols(granges)$length <- y$length
  return(granges)
}
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
grs <- gr_ALK("AVENIO ALK partners.txt")
grs

x #Granges object returned by gr_ALK()
y #name of BAM-file
ALK_partner_peaks <- function(x,y){
  library(GenomicRanges)
  library(Rsamtools)
  library(bamsignals)
  peaks <- c()
  for (i in 1:length(x)){
    peaks[i] <- max(bamCoverage(y, x, mapqual = 0, verbose = F)[i]) 
  }
  z <- data.frame(gene = x$SYMBOL, peak_depth = peaks)
  z <- z[order(-z$peak_depth),]
  return(z)
}
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/ALK patients")
ALK_positive_peaks <- ALK_partner_peaks(grs, "ALK_positive.bam")
ALK_positive_peaks
ctDNA_negative <- ALK_partner_peaks(grs, "ctDNA_negative.bam")
ctDNA_negative
ctDNA_positive <- ALK_partner_peaks(grs, "ctDNA_positive.bam")
ctDNA_positive
max(bamProfile("ALK_positive.bam", grs, verbose = F)[19])
?bamProfile

df <- read.table("AVENIO ALK partners.txt", header = T)
colnames(df) <- c("chr", "start", "end", "name", "strand", "length")
df$start <- as.numeric(df$start)
df$end <- as.numeric(df$end)
df[1,]
df
library(ctDNAtools)
library(GenomicAlignments)
library(BiocParallel)
library(dplyr)
ALK_positive <- bam("ALK_positive.bam", "ALK_positive.bam.bai")
se_input <-summarizeOverlaps(features = grs, reads = ALK_positive, ignore.strand = T, mode = "Union")
assays(se_input)$counts
summarize_fragment_size(bam = "ALK_positive.bam", regions = df[1,])
get_fragment_size(bam = "ALK_positive.bam", targets = df[5,],
                  mapqFilter = 0, ignore_trimmed = F, different_strands = F,
                  simple_cigar = F)$size
bin_fragment_size(bam = "ALK_positive.bam", targets = df[1,])
?get_fragment_size()
which.max(analyze_fragmentation(bam = "ALK_positive.bam", targets = df[1,])$WPS)
analyze_fragmentation(bam = "ALK_positive.bam", targets = df[1,])[25954,]
hist(analyze_fragmentation(bam = "ALK_positive.bam", targets = df[1,])$WPS)
res_EML <- analyze_fragmentation(bam = "ALK_positive.bam", targets = df[1,])
res_NPM <- analyze_fragmentation(bam = "ALK_positive.bam", targets = df[2,])
library(dplyr)
ress <- res_EML %>% filter(WPS < -4)
ress
hist(res_NPM$WPS)
min(res_EML$WPS)
df <- df %>% filter(name != "MSH6")
wpss <- c()
for(i in 1:length(df$chr)){
  res <- analyze_fragmentation(bam = "ALK_positive.bam", targets = df[i,])
  wpss[i] <- min(res$WPS)
  print(i)
}

df
wpss
which.min(wpss)
df[5,]
res <- analyze_fragmentation(bam = "ALK_positive.bam", targets = df[5,])
hist(res$WPS) 
ress <- res %>% filter(WPS < -4)
res