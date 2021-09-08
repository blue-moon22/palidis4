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

### Function to group ITRs into sub clusters based on their position
get_sub_clusters <- function(df) {
  count <- 0
  for (i in 1:length(df$subcluster)) {
    if (df$range[i] > 50) {
      count = count + 1
    }
    df$subcluster[i] <- count
  }
  return(df)
}

### Get ITR positions 
itrs_summary <- itrs %>%
  group_by(sample_id, contig, itr_cluster, exact1, exact2) %>%
  arrange(position1) %>%
  # Get difference in start positions of ISs
  mutate(range = position1 - lag(position1, default = position1[1])) %>% 
  mutate(subcluster = NA) %>%
  group_by(sample_id, contig, itr_cluster, exact1, exact2) %>%
  arrange(position1) %>%
  # Get sub clusters based on start position range
  group_modify(~get_sub_clusters(.x)) %>%
  mutate(subcluster1 = subcluster) %>%
  filter(!is.na(subcluster1)) %>%
  arrange(position2) %>%
  # Get difference in end positions of ISs
  mutate(range =  position2 - lag(position2, default = position2[1])) %>%
  mutate(subcluster = NA) %>%
  arrange(position2) %>%
  # Get sub clusters based on end position range
  group_modify(~get_sub_clusters(.x)) %>%
  rename(subcluster2 = subcluster) %>%
  filter(!is.na(subcluster2)) %>%
  ungroup() %>%
  group_by(sample_id, contig, itr_cluster, subcluster1, subcluster2, exact1, exact2) %>%
  # Get the end positions of ITR1 and ITR2
  summarise(itr1_start_position = min(position1), 
            itr1_end_position = min(position1) + itr1_length - 1, 
            itr2_start_position = max(position2), 
            itr2_end_position = max(position2) + itr2_length - 1) %>%
  # Do the same as above but not ignore grouping by ITR cluster
  group_by(sample_id, contig, exact1, exact2) %>%
  arrange(itr1_start_position) %>%
  mutate(range = itr1_start_position - lag(itr1_start_position, default = itr1_start_position[1])) %>%
  mutate(subcluster = NA) %>%
  group_by(sample_id, contig, exact1, exact2) %>%
  arrange(itr1_start_position) %>%
  group_modify(~get_sub_clusters(.x)) %>%
  mutate(subcluster1 = subcluster) %>%
  filter(!is.na(subcluster1)) %>%
  arrange(itr2_end_position) %>%
  mutate(range =  itr2_start_position - lag(itr2_start_position, default = itr2_start_position[1])) %>%
  mutate(subcluster = NA) %>%
  arrange(itr2_start_position) %>%
  group_modify(~get_sub_clusters(.x)) %>%
  mutate(subcluster2 = subcluster) %>%
  filter(!is.na(subcluster2)) %>%
  ungroup() %>%
  # Join read names
  group_by(sample_id, contig, subcluster1, subcluster2, exact1, exact2) %>%
  summarise(itr1_start_position = min(itr1_start_position),
            itr1_end_position = max(itr1_end_position),
            itr2_start_position = min(itr2_start_position),
            itr2_end_position = max(itr2_end_position),
            itr_clusters = paste(unique(itr_cluster), collapse = ';'),
            read1 = paste(unique(read1), collapse = ';'),
            read2 = paste(unique(read2), collapse = ';'),.groups = 'rowwise') %>%
  ungroup() %>%
  # Clean output
  mutate(itr1_start_position = replace(itr1_start_position, is.na(itr1_end_position), NA),
         itr2_start_position = replace(itr2_start_position, is.na(itr2_end_position), NA)) %>%
  select(sample_id, contig, itr1_start_position, itr1_end_position, itr2_start_position, itr2_end_position, itr_clusters, read1, read2)

# Get read names
reads_with_itr_cluster <- itrs[itrs$itr_cluster %in% unlist(strsplit(itrs_summary$itr_clusters, ";")),] %>% select(itr_cluster, read1, read2)
ir_read_names <- c(reads_with_itr_cluster$read1, reads_with_itr_cluster$read2)

### Write output
write.table(itrs_summary, paste0(opt$output, "_insertion_sequence_annotations.tab"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(ir_read_names, paste0(opt$output, "_read_names.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
write.table(reads_with_itr_cluster, paste0(opt$output, "_reads_itr_clusters.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
