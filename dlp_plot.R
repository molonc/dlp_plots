require('ape')
require('argparse')
require('dplyr')
require('ggplot2')
require('ggtree')
require('gtools')
require('RColorBrewer')
require('devtools')
require('ComplexHeatmap')
require("data.table")

options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)

input_directory <- "/projects/molonc/scratch/sbeatty/SCY-298/data/A144168B"
metadata <- "/projects/molonc/scratch/sbeatty/SCY-298/metadata/leap114groupings.csv"


hdbscran_viz <- function(filtered_CNV_fn, reads_fn,
                         save_dir, datatag, grouping_file, 
                         pctcells = 0.05, min_dist_grain_level = 0.05, 
                         remove_outliers=TRUE, probJump=0.995, probAVG=0.995){

### identify input files  
hmmcopy_reads_filtered_states_path <- paste0(input_directory,"/","hmmcopy/reads_filtered.csv.gz")


### read in input files 
cnv <- data.table::fread(hmmcopy_reads_filtered_states_path) %>% as.data.frame()
rownames(cnv) 
chr_segment <-  paste(cnv$chr,cnv$start, cnv$end,sep="_") 

 ##   Note: format of cell id is: sampleId-libraryId-rowId-colId, ex: "AT11391-A98166B-R45-C08"
filtered_cells <- colnames(cnv)
  
  ##   Note: format of chr region is: chr_start_end, ex: '1_2000001_2500000'
filtered_chr_regions <- chr_segment


#### might need to add chr_segement to copynumber/cnv
list_of_filtered_hmmcopy_cnv_data <- list(copynumber=cnv,
              filtered_cells=filtered_cells, 
              filtered_chr_regions=filtered_chr_regions)
copynumber <- list_of_filtered_hmmcopy_cnv_data$copynumber
cells_used <- colnames(copynumber)


compute_jump_cells <- function(mat) {
    # remove odd cells
    # Does not count whole-chromosome events
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


remove_outliers <- TRUE
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
