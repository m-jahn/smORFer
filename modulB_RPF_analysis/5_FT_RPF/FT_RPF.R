# load library
library(tidyverse)
library(seqinr)
library(Biostrings)

# get parameters
get_param <- commandArgs(trailingOnly = TRUE)

# read RPF counts per nucleotide
rpf <- read.table(get_param[1], sep = "\t", header = F, stringsAsFactors = F)
# get genes
u_genes <- unique(rpf[,1:6])

# cutoff for selection of candidates: (FT at period 3) / (FT between 1.5 and 3) 
cutoff <- 0

## function for fourier transform input counts
# counts: input count
# it will only take a vector which is of length 3 of dividable by 3, others are skipped
# output(in default mode): 4 values (FT at 1.5, FT at 3, mean of FT between 1.5 and 3, mean > 3) 
# normalze: per default it normalizes the transformed values by length (which we have to do)
# full_output: if TRUE will output list with all data, e.g. to plot 
get_FT_signal = function(counts, normalize = T, full_output = F){
  # check if length diviable by 3
  if( length(counts)%%3 == 0){
    # get frequency (x-axis) for FT transform
    # this is a bit complicated to explain and understand (we can try later)
    freq <- length(counts)/(0:(length(counts)-1))
    # perform fast FT
    ft <- abs(fft(counts))
    # normalize by length
    if(normalize){
      ft <- ft / length(counts)
    }
    # get identity of period 3 and 1.5 together with mean of inter-regions
    idx3 <- which(freq == 3)
    idx15 <- which(freq == 1.5)
    res <- list(
      ft15 = ft[idx15],
      ft3 = ft[idx3],
      mean_ftlower3 = mean(ft[(idx3+1):(idx15-1)]),
      mean_fthigher3 = mean(ft[(idx15+1):(length(ft))])
    )
    # return
    if(!full_output){
      res
    } else {
      list(
        res = res,
        ft = ft,
        freq = freq
      )
    }
    
    # else skip with message
  } else {
    cat('skipped\n')
    NULL
  }
}


### calculate FT for RPFs
verified <- rpf %>%
  as_tibble %>%
  group_by(V1, V2, V3, V4, V5, V6) %>%
  summarize(ft_result = list(get_FT_signal(V8))) %>%
  tidyr::unnest_wider(ft_result) %>%
  mutate(FT_ratio = ft3/mean_ftlower3) %>%
  filter(FT_ratio > cutoff) %>%
  select(-ft15, -ft3, -mean_ftlower3, -mean_fthigher3)

# write candidates
write.table(verified, get_param[2], sep = "\t", col.names = F, row.names = F, quote = F)
