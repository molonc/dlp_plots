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



main <- function(){
  ##   Note: please change the directory to your script dir here. 
  script_dir <- '/home/htran/Projects/hakwoo_project/git_repos/dlp_hdbscan_clustering/hdbscan_clustering/'
  source(paste0(script_dir,'hdbscan_clustering_utils.R')) # clustering function, load data func
  source(paste0(script_dir,'make_cell_copynumber_tree_heatmap.R')) # all visualization func
  
  ##  Note: you can change the input directory to point to hdbscan_clustering/testing_data/ 
  # results_dir <- '/home/htran/Projects/hakwoo_project/git_repos/dlp_hdbscan_clustering/hdbscan_clustering/testing_data/'
  # save_dir <- paste0(results_dir,'output/')
  # grouping_file <- paste0(results_dir,'library_grouping.csv') # meta samples
  # reads_fn <- paste0(results_dir,'A98166B_reads.csv.gz') # input data from hmmcopy output
  # datatag <- 'A98166B'
  # filtered_CNV_fn <- paste0(results_dir,'A98166B_filtered_states.csv.gz')
  # 
  # hdbscran_viz(filtered_CNV_fn, reads_fn,
  #              save_dir, datatag, grouping_file, pctcells = 0.05)
  
  download_dir <- '/home/htran/storage/raw_DLP/LEAP_DLP/'
  results_dir <- '/home/htran/storage/LEAP_dlp_results/'
  save_dir <- paste0(results_dir,'hdbscan_clustering_results_200K/')
  
  libs <- list.files(paste0(results_dir,'filtered_data_200K/'))
  libs <- gsub('_filtered_states.csv.gz','',libs)
  
  length(libs)
  
  ## Note: 
  # Library id: A144171A, sample id: AT23026 (not id A144171B as original one) due to pooled samples
  # Library id: A144171A, sample id: AT23025
  # → Building 2 heatmaps for 2 different sample ids
  libs <- libs[!libs %in% c('A144171A-AT23025','A144171A-AT23026')] 
  
  
  # TO DO: add filtering outliers function here, filtered by medianCNV
  # datatag <- 'A138848B'
  
  ## For this library, can we get sub set of cells in a clone, and dividing cells into 2 clusters by repeat hdbscan
  datatag <- 'A138979A'
  
  ## To Do: check clustering grain level for this library: 'A138850A'
  datatag <- 'A138848B' # Done
  datatag <- 'A138850A'# Done partially, need to recut the None clone
  datatag <- 'A144169A' # check if min_dist can help to divide data into smaller clones
  datatag <- 'A144168A'
  libs <- libs[!libs %in% c('A138848B','A138850A')] 
  libs
  for(datatag in libs){
    print(datatag)
    grouping_file <- paste0(results_dir,'metadata/','library_groupings_',datatag,'.csv') # meta samples
    # reads_fn <- paste0(raw_data_dir,datatag,'/hmmcopy/',datatag,'_reads.csv.gz') # input data from hmmcopy output
    filtered_CNV_fn <- paste0(results_dir,'filtered_data_200K/',datatag,'_filtered_states.csv.gz')
    
    l <- datatag
    lib_fn1 <- paste0(download_dir,l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    lib_fn2 <- paste0(download_dir,l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
    lib_fn <- ''
    if(file.exists(lib_fn1)){
      lib_fn <- lib_fn1
    } else{
      if(file.exists(lib_fn2)){
        lib_fn <- lib_fn2
      } else{
        print(paste0('***ERROR: library: ', l))
        stop('Do not exist hmmcopy reads for library, double check input data')
      } 
    } 
    if(lib_fn!=''){
      hdbscran_viz(filtered_CNV_fn, lib_fn,
                   save_dir, datatag, grouping_file, 
                   pctcells = 0.05, min_dist_grain_level = 0.1,
                   remove_outliers=TRUE, probJump=0.995, probAVG=0.995)
    }  
  }
  
  
  
  ## Very special cases
  libs <- c('A144171A-AT23025','A144171A-AT23026')
  
  for(datatag in libs){
    print(datatag)
    l <- 'A144171A'
    grouping_file <- paste0(results_dir,'metadata/','library_groupings_',datatag,'.csv') # meta samples
    # reads_fn <- paste0(raw_data_dir,datatag,'/hmmcopy/',datatag,'_reads.csv.gz') # input data from hmmcopy output
    filtered_CNV_fn <- paste0(results_dir,'filtered_data_200K/',datatag,'_filtered_states.csv.gz')
    
    
    lib_fn1 <- paste0(download_dir,l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
    lib_fn2 <- paste0(download_dir,l,'/hmmcopy/',l,'_reads.csv.gz') # csv.gz or csv
    lib_fn <- ''
    if(file.exists(lib_fn1)){
      lib_fn <- lib_fn1
    } else{
      if(file.exists(lib_fn2)){
        lib_fn <- lib_fn2
      } else{
        print(paste0('***ERROR: library: ', l))
        stop('Do not exist hmmcopy reads for library, double check input data')
      } 
    } 
    if(lib_fn!=''){
      hdbscran_viz(filtered_CNV_fn, lib_fn,
                   save_dir, datatag, grouping_file, 
                   pctcells = 0.05, min_dist_grain_level = 0.1,
                   remove_outliers=TRUE, probJump=0.995, probAVG=0.995)
    }  
  }
  
}
# main()
