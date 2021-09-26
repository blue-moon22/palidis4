#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(igraph)
library(stringr)
library(tidyr)

option_list = list(
  make_option(c("-c", "--clusters"), type="character", default=NULL, 
              help="Tab-delimited file of CD-HIT clusters", metavar="character"),
  make_option(c("-t", "--txt"), type="character", default=NULL, 
              help="Suffix of txt files with reads and non-catalog ITR clusters e.g. _reads_itr_clusters.txt", metavar="character"),
  make_option(c("-a", "--annot"), type="character", default=NULL, 
              help="Suffix of insertion sequence annotations e.g. _insertion_sequence_annotations.tab", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output prefix of IS annotations catalog", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$clusters)){
  print_help(opt_parser)
  stop("One argument must be supplied for -c/--clusters (input file).n", call.=FALSE)
}

if (is.null(opt$txt)){
  print_help(opt_parser)
  stop("One argument must be supplied for -t/--txt (pattern).n", call.=FALSE)
}

if (is.null(opt$annot)){
  print_help(opt_parser)
  stop("One argument must be supplied for -a/--annot (pattern).n", call.=FALSE)
}

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("One argument must be supplied for -o/--output (output prefix).n", call.=FALSE)
}

# Read tab file of cd-hit clusters
itr_clusters <- read.table(opt$clusters, header = TRUE)

# Read read + old itr cluster text files
itr_pairs <- data.frame()
txt_files <- list.files(pattern = opt$txt)
for (i in 1:length(txt_files)) {
  tmp <- read.table(txt_files[i], header = TRUE)
  tmp$sample_id <- gsub(opt$txt, "", txt_files[i])
  itr_pairs <- rbind(itr_pairs, tmp)
}

# Assign new clusters
new_itr_clusters <- left_join(itr_pairs, itr_clusters, by = c("read1"="seq")) %>%
  left_join(itr_clusters, by = c("read2"="seq")) %>%
  filter(!(is.na(cd_hit_cluster.x) & is.na(cd_hit_cluster.y))) %>%
  group_by(cd_hit_cluster.x, cd_hit_cluster.y) %>%
  mutate(cd_hit_cluster = paste(unique(c(cd_hit_cluster.x, cd_hit_cluster.y)), collapse = ';')) %>%
  mutate(cd_hit_cluster = gsub("NA;", "", cd_hit_cluster),
         cd_hit_cluster = gsub(";NA", "", cd_hit_cluster)) %>%
  ungroup()
  

# Re-assign multiple clusters
new_itr_clusters2 <- data.frame()
for(i in 1:length(new_itr_clusters$cd_hit_cluster)){
  itr_clusters <- str_split(new_itr_clusters$cd_hit_cluster[i], ";")[[1]]
  for(j in 1:length(itr_clusters)){
    new_itr_clusters2 <- rbind(new_itr_clusters2, c(itr_clusters[1], itr_clusters[j]))
  }
}
clusters <- clusters(graph.data.frame(new_itr_clusters2))
itr_clusters <- with(clusters, 
                      data.frame(
                        id = names(membership), 
                        group = membership, 
                        number_in_group = csize[membership]
                      )
) %>% arrange(group)

new_itr_pairs <- new_itr_clusters %>%
  separate_rows(cd_hit_cluster, sep = ";") %>%
  left_join(itr_clusters, by = c("cd_hit_cluster"="id")) %>%
  mutate(itr_cluster = as.character(itr_cluster)) %>%
  select(sample_id, itr_cluster, group) %>%
  rename(itr_cluster_catalog = group) %>%
  unique()

# Open all tab files and re-assign ITR cluster
tab_files <- list.files(pattern = opt$annot)

annot <- data.frame()
for (i in 1:length(tab_files)){
  tmp <- read.table(tab_files[i], header = TRUE)  
  annot <- rbind(annot, tmp)
}

# Add new itr clusters
new_annots <- annot %>%
  inner_join(new_itr_pairs, by = c("sample_id", "itr_cluster")) %>%
  select(sample_id, contig, itr1_start_position, itr1_end_position, itr2_start_position, itr2_end_position, itr_cluster_catalog) %>%
  unique()

write.table(new_annots, paste0(opt$output, "_insertion_sequence_annotations_catalog.tab"), row.names = FALSE, quote = FALSE, sep = "\t")