
### Hdbscan clustering the most updated version
Installation: to install R packages, see the first part of *hdbscan_clustering_utils.R* file

Run demo: see the file: *hdbscan_clustering.R*:

- First, you need to load utility functions from file:

- Load clustering function, load data func: script_dir <- 'yourDir/dlp_hdbscan_clustering/hdbscan_clustering/'
source(paste0(script_dir,'hdbscan_clustering_utils.R')) 
- Load all visualization func: source(paste0(script_dir,'make_cell_copynumber_tree_heatmap.R')) 
- To run a demo in this repo, you need to set directory to testding_data/ folder: 
  results_dir <- 'yourDir/dlp_hdbscan_clustering/hdbscan_clustering/testing_data/'
  
  save_dir <- paste0(results_dir,'output/')
  
- Meta samples file, see file for column names format:   
  grouping_file <- paste0(results_dir,'library_grouping.csv') # meta samples
  
- Loading raw reads file from hmmcopy input data:   
  reads_fn <- paste0(results_dir,'A98166B_reads.csv.gz') 
  
- Set other params:   
  
  datatag <- 'A98166B'
  
  filtered_CNV_fn <- paste0(results_dir,'A98166B_filtered_states.csv.gz')
  
- And execute function from hdbscan_clustering.R file:  

  hdbscran_viz(filtered_CNV_fn, reads_fn,
               save_dir, datatag, grouping_file, pctcells = 0.05)


*Detail*: samples were clustered by first applying UMAP to the normalized signature probabilities for the HMMcopy copy number states/ or copy number reads with n_neighbors = 20 (or your input) and min_dist = 0 (or your input, i.e 0.1) to produce two-dimensional sample embeddings. Next, HDBSCAN was run on the sample embeddings with min_samples = 5, min_cluster_size = 5 and cluster_selection_ epsilon = 0.75 (or your input) to produce the sample clusters.

See more description and results of this method at: *Single-cell genomic variation induced by mutational processes in cancer* :https://pubmed.ncbi.nlm.nih.gov/36289342/

### Method references: 
- https://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html
- https://pberba.github.io/stats/2020/01/17/hdbscan/
- R package with documentation and list of params: https://cran.r-project.org/web/packages/dbscan/dbscan.pdf

Have fun +_+
