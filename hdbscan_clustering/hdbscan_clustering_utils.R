## Installation
## Install the following package to visualize data. 
# BiocManager::install('ape')
# BiocManager::install('argparse')
# BiocManager::install('dplyr')
# BiocManager::install('ggplot2')
# BiocManager::install('ggtree')
# BiocManager::install('gtools')
# BiocManager::install('RColorBrewer')
# BiocManager::install('devtools')
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")

## Note: here are the package version that I used to plot figure
## but I think the script can work with almost other versions. 
## For ComplexHeatmap, I recommend to use the latest version 
# [1] "ape: 5.6.2"
# [1] "argparse: 2.1.6"
# [1] "ComplexHeatmap: 2.8.0"
# [1] "dplyr: 1.0.9"
# [1] "ggplot2: 3.3.6"
# [1] "ggtree: 3.0.4"
# [1] "gtools: 3.9.3"
# [1] "RColorBrewer: 1.1.3"


suppressPackageStartupMessages({
  library(dplyr)
  library(signals)
})  
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)

## Reference: please cite the signals package if you use these scripts
## https://github.com/shahcompbio/signals
## Noted: original version from Tyler Funnel 
## Modified scripts by Hoa Tran
## Hoa: adding some scripts to run this function but main functions are from signal packgage, Marc William
## Note: comment main() function in this file








##   Note: copy number file name with row name is chr regions, and column name is cell id
## filtered_CNV_fn file has this format: 
#                         AT11391-A98166B-R45-C08
# 1_2000001_2500000                       1
# 1_3000001_3500000                       1
# 1_4000001_4500000                       1

## Note: usually we have 4375 filtered genomic bins, but sometimes, you may exclude some noisy bins from a reference list, 
## and have ~3950 bins for analysis 

get_filtered_data <- function(filtered_CNV_fn){
  cnv <- data.table::fread(filtered_CNV_fn) %>% as.data.frame()
  rownames(cnv) <- cnv$V1 # chr regions are row names
  cnv$V1 <- NULL
  print(dim(cnv))
  ##   Note: format of cell id is: sampleId-libraryId-rowId-colId, ex: "AT11391-A98166B-R45-C08"
  filtered_cells <- colnames(cnv)
  
  ##   Note: format of chr region is: chr_start_end, ex: '1_2000001_2500000'
  filtered_chr_regions <- rownames(cnv)
  return(list(copynumber=cnv,
              filtered_cells=filtered_cells, 
              filtered_chr_regions=filtered_chr_regions))
}



## pctcells: percentage of minimum cells to consider as one clone
# datatag <- 'A138848B'
# reads_fn <- lib_fn
# pctcells = 0.05
# min_dist_grain_level = 0.05
# remove_outliers=TRUE
# probJump=0.99
# probAVG=0.975

# datatag <- 'A138850A'
# reads_fn <- lib_fn
# pctcells = 0.05
# min_dist_grain_level = 0.01
# remove_outliers=FALSE # if not, removing small interesting population here
# probJump=0.99
# probAVG=0.995

# datatag <- 'A138979A' # noisy clusters, need repeat hdbscan
# reads_fn <- lib_fn
# pctcells = 0.05
# min_dist_grain_level = 0.1
# remove_outliers=TRUE # if not, removing small interesting population here
# probJump=0.99
# probAVG=0.99

datatag <- 'A144171A-AT23025' # noisy clusters, need repeat hdbscan
# datatag <- 'A144171A-AT23026'
reads_fn <- lib_fn
pctcells = 0.05
min_dist_grain_level = 0.1
remove_outliers=TRUE # if not, removing small interesting population here
probJump=0.995
probAVG=0.995

datatag <- 'A144168A'
reads_fn <- lib_fn
pctcells = 0.05
min_dist_grain_level = 0.1
remove_outliers=TRUE # if not, removing small interesting population here
probJump=0.995
probAVG=0.995
hdbscran_viz <- function(filtered_CNV_fn, reads_fn,
                         save_dir, datatag, grouping_file, 
                         pctcells = 0.05, min_dist_grain_level = 0.05, 
                         remove_outliers=TRUE, probJump=0.995, probAVG=0.995){
  
  
  res <- get_filtered_data(filtered_CNV_fn)
  copynumber <- res$copynumber
  if(!dir.exists(save_dir)){
    dir.create(save_dir)
  }
  cells_used <- colnames(copynumber)
  
  if(remove_outliers){
    jump_rank <- compute_jump_cells(copynumber)
    # # Filter cells by jump/stuff
    bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = probJump) 
                                              | avgCNA > quantile(x=avgCNA, probs = probAVG)) %>% dplyr::select(cell_id)
    # bad_cells <- jump_rank %>% dplyr::filter( njumps > quantile(x=njumps, probs = prob)) %>% dplyr::select(cell_id)
    
    bad_cells <- bad_cells$cell_id
    print(paste0("Number of bad cells (outliers jump or large CNA values): ",length(bad_cells)))
    
    
    cells_used <- cells_used[!cells_used %in% bad_cells]
    print(paste0('Number of filtered cells after removing outlier cells is: ', length(cells_used)))
    copynumber <- copynumber[,cells_used]
    # dim(copynumber)
  }
  
  
  ## To Do: removing bad cells here
  
  ## Input hmmcopy reads file
  reads_df <- data.table::fread(reads_fn) %>% as.data.frame()
  # reads_df <- data.table::fread(paste0(download_dir,library_id,'/hmmcopy/',library_id,'_reads.csv.gz')) %>% as.data.frame()
  # print(dim(reads_df))
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% cells_used)
  
  reads_df$chr_desc <- paste0(reads_df$chr,'_',reads_df$start,'_',reads_df$end)
  reads_df <- reads_df %>%
    dplyr::filter(chr_desc %in% res$filtered_chr_regions)
  
  # print(dim(reads_df))
  ncells <- length(unique(reads_df$cell_id))
  
  
  ##Note:   here I use field = "copy", copy number values as the input for cell clustering,
  ## instead of copy number state at 'state' column
  ## you can use different columns, and comparing the results of clustering
  min_cell_per_cluster <- 10
  
  # min_dist_grain_level <- 0.1
  ## min_dist: minimum distance between embedded points, smaller values more clusters, default 0.1
  clusters <- signals::umap_clustering(reads_df,
                                        minPts = max(round(pctcells * ncells), min_cell_per_cluster),
                                        field = "copy",
                                        min_dist = min_dist_grain_level)
  
  # clusters <- signals::umap_clustering(reads_df, 
  #                                      minPts = min(round(pctcells * ncells), 30), 
  #                                      field = "copy")
  
  tree <- clusters$tree
  clones <- clusters$clustering
  
  print(summary(as.factor(clones$clone_id)))
  res_clones <- unique(clones$clone_id)
  # print(res_clones)
  if('0' %in% res_clones & !'None' %in% res_clones){
    clones$clone_id <- ifelse(clones$clone_id=='0','None',clones$clone_id)
  }

  
  output_fn <- paste0(save_dir, datatag,'_hdbscan_cell_cn_tree_heatmap.png')
  png(output_fn, height = 2*700, width=2*1200, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  dev.off()
  
  
  # data.table::fwrite(clones, paste0(save_dir,datatag,'_',added_time,'_cell_clones.csv'))
  data.table::fwrite(clones, paste0(save_dir,datatag,'_cell_clones.csv'))
  
  # added_time <- gsub(':','',format(Sys.time(), "%Y%b%d_%X"))
  # tree_fn <- paste0(save_dir,datatag,'_',added_time,'_tree.newick')
  tree_fn <- paste0(save_dir,datatag,'_tree.newick')
  ape::write.tree(phy = tree, file = tree_fn, tree.names = F)
  
  # saveRDS(tree, paste0(save_dir,datatag,'_',added_time,'_tree.rds'))
  
}

# remove odd cells
# Does not count whole-chromosome events
compute_jump_cells <- function(mat) {
  # Cells that have too many CN up and downs may represent douplet? or s-phases
  # Exclude them after a threshold
  # dat <- load_new_cn_data(datatag)
  dat <- sort_mat_by_bins(mat)
  
  # Bins X cells
  mat_delta <- abs(dat[1:(nrow(dat)-1), ] - dat[2:nrow(dat), ])
  mat_delta[mat_delta > 1] <- 1
  njumps <- colSums(mat_delta)
  avgCNA <- colMeans(dat)
  # res <- sort(njumps, decreasing = T)
  stopifnot(colnames(mat_delta) == colnames(dat))
  data.frame(cell_id = colnames(dat), njumps = njumps, avgCNA = avgCNA, stringsAsFactors = F) %>% dplyr::arrange(desc(njumps))
}

### Matrix functions
################################################################################################################################################
sort_mat_by_bins <- function(the_mat) {
  # prevent scientific notation
  options(scipen=999)
  options(stringsAsFactors=FALSE) 
  cnv_txt <- parse_bin_names(rownames(the_mat), as_factors = F)
  the_mat <- cbind(cnv_txt, the_mat)
  
  # Sort the matrix by their chromosome (just inside the chromosome)
  the_mat$chr[the_mat$chr == 'X'] <- '40' 
  the_mat$chr[the_mat$chr == 'Y'] <- '55' 
  the_mat$chr <- as.numeric(the_mat$chr)
  the_mat <- the_mat[order(the_mat$chr, the_mat$start), ]
  the_mat$chr[the_mat$chr == '40'] <- 'X' 
  the_mat$chr[the_mat$chr == '55'] <- 'Y' 
  
  # Remove chr, start, end
  the_mat$chr <- NULL
  the_mat$start <- NULL
  the_mat$end <- NULL
  
  the_mat
}

parse_bin_names <- function(bin_names, as_factors = FALSE) {
  # Remove corrupt_tree locus tag if it's there
  bin_names <- gsub('locus_', '', bin_names)
  chr <- gsub('([0-9]+|X|Y)_[0-9]+_[0-9]+', '\\1', bin_names)
  start <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_[0-9]+', '\\2', bin_names))
  end <- as.numeric(gsub('([0-9]+|X|Y)_([0-9]+)_([0-9]+)', '\\3', bin_names))
  data.frame(chr = chr, start = start, end = end, stringsAsFactors = as_factors)
}


