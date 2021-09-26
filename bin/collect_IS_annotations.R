#!/usr/bin/env Rscript
library(optparse)
library(dplyr)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Tab-delimited file of ITR positions", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output prefix of IS annotations", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("One argument must be supplied for -i/--input (input file).n", call.=FALSE)
}

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("One argument must be supplied for -o/--out (output file).n", call.=FALSE)
}

### Read Palidis data
itrs <- read.delim(opt$input)
swap <- itrs$position1 - itrs$position2 > 0
position1 <- itrs$position1[swap]
read1 <- itrs$read1[swap]
position2 <- itrs$position2[swap]
read2 <- itrs$read2[swap]
itrs$position1[swap] <- position2
itrs$read1[swap] <- read2
itrs$position2[swap] <- position1
itrs$read2[swap] <- read1

### Get exact positions and estimates, filter two estimates
itrs$exact1 <- grepl('_LCoord_', itrs$read1)
itrs$exact2 <- grepl('_LCoord_', itrs$read2)
itrs <- itrs[!(!itrs$exact1 & !itrs$exact2),]

### Get lengths of ITRs
itrs$itr1_length <- sapply(strsplit(gsub(".*_LCoord_", "", itrs$read1), "_RCoord_"), function(.x) as.numeric(.x[2]) + 1 - as.numeric(.x[1]))
itrs$itr2_length <- sapply(strsplit(gsub(".*_LCoord_", "", itrs$read2), "_RCoord_"), function(.x) as.numeric(.x[2]) + 1 - as.numeric(.x[1]))

### Filter lengths of ITRs above 50
if (length(which(itrs$itr1_length > 50)) > 0) {
  itrs <- itrs[-which(itrs$itr1_length > 50),]
}
if (length(which(itrs$itr2_length > 50)) > 0) {
  itrs <- itrs[-which(itrs$itr2_length > 50),]
}

# Get read names
reads_with_itr_cluster <- itrs %>% select(itr_cluster, read1, read2)

# Get start and end positions and summary
itrs_summary <- itrs %>%
  select(sample_id, contig, position1, position2, exact1, exact2, itr1_length, itr2_length, itr_cluster) %>%
  unique() %>%
  mutate(itr1_start_position = position1, 
         itr1_end_position = position1 + itr1_length - 1, 
         itr2_start_position = position2, 
         itr2_end_position = position2 + itr2_length - 1) %>%
  mutate(itr1_start_position = replace(itr1_start_position, is.na(itr1_end_position), NA),
         itr2_start_position = replace(itr2_start_position, is.na(itr2_end_position), NA)) %>%
  select(sample_id, contig, itr1_start_position, itr1_end_position, itr2_start_position, itr2_end_position, itr_cluster) 

### Write output
write.table(itrs_summary, paste0(opt$output, "_insertion_sequence_annotations.tab"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(reads_with_itr_cluster, paste0(opt$output, "_reads_itr_clusters.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
