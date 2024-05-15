## Installation
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

cn_colours <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
  ),
  names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
)

clone_palette_20 <- c(
  "#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
  "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
  "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
  "#cb48cb", "#d74391"
)
clone_none_black <- "#1B1B1B"


## START hdbscan_clustering_utils below
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

datatag <- 'A144168A'
# reads_fn <- lib_fn

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
  

  tree_fn <- paste0(save_dir,datatag,'_tree.newick')
  ape::write.tree(phy = tree, file = tree_fn, tree.names = F)
  
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

## THE END of hdbscan_clustering_utils ^
## BEGIN make_cell_copynumber_heatmap below

get_args <- function() {
  p <- ArgumentParser(description="Plot cell copynumber heatmap with tree")
  
  p$add_argument("--tree", "-t", help="cell newick tree file")
  p$add_argument("--copynumber", "-cn", help="cell copynumber tsv file")
  p$add_argument("--pdf", "-o", help="output plot pdf")
  
  p$add_argument("-c", "--clones", help="cell clone tsv file")
  p$add_argument(
    "-n", "--normalize-ploidy", action="store_true",
    help="normalize all cell ploidy to 2"
  )
  p$add_argument(
    "-b", "--branch-lengths", type="integer",
    help="set all tree branch lengths to this value"
  )
  
  p$add_argument("--grouping_file", default=NULL, help="file with sample grouping")
  
  return(p$parse_args())
}

read_tsv <- function(fn, ...) {
  df <- read.delim(fn, check.names=FALSE, stringsAsFactors=FALSE, sep=",", ...)
  return(df)
}

calc_state_mode <- function(states) {
  state_levels <- unique(states)
  state_mode <- state_levels[
    which.max(tabulate(match(states, state_levels)))
  ]
  return(state_mode)
}

normalize_cell_ploidy <- function(copynumber) {
  cell_ids <- colnames(copynumber)
  cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]
  
  for(cell_id in cell_ids) {
    state_mode <- calc_state_mode(copynumber[[cell_id]])
    copynumber[[cell_id]] <- as.integer(ceiling(
      copynumber[[cell_id]] / (state_mode / 2)
    ))
  }
  return(copynumber)
}

format_tree <- function(tree, brlen) {
  locus_tips <- grep('locus', tree$tip.label, value=TRUE)
  tree <- drop.tip(tree, locus_tips)
  
  if (!is.null(brlen)) {
    tree <- compute.brlen(tree, brlen)
  }
  
  tree$tip.label <- gsub('cell_', '', tree$tip.label)
  
  return(tree)
}

get_clone_members <- function(clones) {
  clone_members <- list()
  for(c in unique(clones$clone_id)) {
    if(c != "None") {
      clone_members[[c]] <- clones[clones$clone_id == c, "cell_id"]
    }
  }
  return(clone_members)
}

# Just modify function for predefined clone colors. 
# cell_clones <- read.csv('/home/htran/storage/datasets/drug_resistance/dlp_results/SA535/SA535_total_v2/cell_clones.csv', check.names = F, stringsAsFactors = F)
# levels <- unique(cell_clones$clone_id)
make_clone_palette <- function(levels) {
  clone_names <- sort(levels)
  
  ## Note: not easy to install inlmisc package, so I modified scripts here
  # install.packages("inlmisc", dependencies = TRUE)
  # pal <- as.character(inlmisc::GetColors(length(clone_names)))
  if(length(levels) <= 8) {
    pal <- brewer.pal(max(length(levels), 3), "Dark2")
  }
  else if (length(levels) > 8 & length(levels) <= 12) {
    pal <- brewer.pal(max(length(levels), 3), "Set3")
  } else if (length(levels) > 12 & length(levels) <= 20) {
    pal <- clone_palette_20
  } else {
    pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
    print("WARNING: more clones than palette can accomodate!, using extra colors")
  }
  
  names(pal) <- clone_names
  pal <- pal[levels]
  return(pal)
}
# make_clone_palette <- function(levels) {
#     if (length(levels) <= 12) {
#         pal <- brewer.pal(max(length(levels), 3), "Set3")
#     } else if (length(levels) <= 20) {
#         pal <- clone_palette_20
#     } else {
#         pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels))
#         print("WARNING: more clones than palette can accomodate!")
#     }
#     names(pal) <- levels
#     pal <- pal[levels]
#     return(pal)
# }

make_tree_ggplot <- function(tree, clones) {
  if(!is.null(clones)) {
    clone_members <- get_clone_members(clones)
    tree <- groupOTU(tree, clone_members) #ggtree::
    
    clone_levels <- mixedsort(unique(clones$clone_id))
    clone_pal <- make_clone_palette(clone_levels)
    clone_pal[["0"]] <- clone_none_black
    
    tree_aes <- aes(x, y, colour=group)
  } else {
    tree_aes <- aes(x, y)
  }
  
  p <- ggplot(tree, tree_aes) +
    geom_tree(size=0.25) +
    coord_cartesian(expand=FALSE) +
    ylim(0.5, length(tree$tip.label) + 0.5) +
    theme_void()
  
  
  if(!is.null(clones)) {
    p <- p + scale_colour_manual(values=clone_pal)
  }
  
  return(p)
}

make_discrete_palette <- function(pal_name, levels) {
  if (length(levels) <= 8) {
    pal <- brewer.pal(max(length(levels), 3), pal_name)
  } else if (length(levels) <= 12) {
    pal <- brewer.pal(max(length(levels), 3), "Set3")
  } else if (length(levels) <= 20) {
    pal <- clone_palette_20
  } else {
    pal <- clone_palette_20
    print("WARNING: more clones than palette can accomodate!")
  }
  pal <- pal[1:length(levels)]  # in case levels contains less than 3 categories
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}


format_copynumber_values <- function(copynumber) {
  copynumber[copynumber > 11] <- 11
  for(col in colnames(copynumber)) {
    values <- as.character(copynumber[, col])
    values[values == "11"] <- "11+"
    copynumber[, col] <- values
  }
  return(copynumber)
}

space_copynumber_columns <- function(copynumber, spacer_cols) {
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  print("205: DEBUG")
  spacer <- as.data.frame(matrix(
    data=NA, nrow=nrow(copynumber), ncol=spacer_cols
  ))
  chrom_copynumber_dfs <- list()
  for(chrom in mixedsort(unique(chroms))) {
    chrom_copynumber <- copynumber[, chroms == chrom, drop=FALSE]
    chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(chrom_copynumber))
    chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(spacer))
  }
  chrom_copynumber_dfs[length(chrom_copynumber_dfs)] <- NULL
  copynumber <- do.call(cbind, chrom_copynumber_dfs)
  
  return(copynumber)
}

get_ordered_cell_ids <- function(tree_plot_dat) {
  return(rev(arrange(tree_plot_dat[tree_plot_dat$isTip, ], y)$label))
}

format_copynumber <- function(copynumber, tree_plot_dat, spacer_cols=20) {
  if (!("chr" %in% colnames(copynumber))) {
    loci <- sapply(rownames(copynumber), strsplit, "_")
    copynumber$chr <- unname(sapply(loci, '[[', 1))
    copynumber$start <- as.numeric(unname(sapply(loci, '[[', 2)))
    copynumber$end <- as.numeric(unname(sapply(loci, '[[', 3)))
    copynumber$width <- (copynumber$end - copynumber$start + 1)
  }
  copynumber$chr <- gsub("chr", "", copynumber$chr)
  # sum(is.na(copynumber$chr))
  
  chr_val <- seq(1:24)
  chrs <- c(as.character(seq(1:22)),'X','Y')
  desc_chrs <- data.frame(chr=chrs, chr_val=chr_val)
  copynumber <- copynumber %>%
    left_join(desc_chrs, by='chr')
  copynumber <- arrange(copynumber, as.numeric(chr_val), chr_val, start)
  
  rownames(copynumber) <- paste0(
    copynumber$chr, ":", copynumber$start, ":", copynumber$end
  )
  copynumber <- subset(copynumber, select=-c(chr, chr_val, start, end, width))
  copynumber <- as.data.frame(t(copynumber))
  
  copynumber <- copynumber[get_ordered_cell_ids(tree_plot_dat), ]
  
  copynumber <- format_copynumber_values(copynumber)
  copynumber <- space_copynumber_columns(copynumber, spacer_cols)
  
  return(copynumber)
}

format_clones <- function(clones, tree_plot_dat) {
  tree_cells <- get_ordered_cell_ids(tree_plot_dat)
  clones <- merge(clones, data.frame(cell_id=tree_cells), all=TRUE)
  clones[is.na(clones$clone_id), "clone_id"] <- "None"
  
  clone_counts <- clones %>% group_by(clone_id) %>% summarise(count=n())
  
  for(i in 1:nrow(clone_counts)) {
    clone <- unlist(clone_counts[i, "clone_id"], use.names=FALSE)
    clone_count <- unlist(clone_counts[i, "count"], use.names=FALSE)
    clone_label <- paste0(clone, " (", clone_count, ")")
    clones[clones$clone_id == clone, "clone_label"] <- clone_label
  }
  
  rownames(clones) <- clones$cell_id
  clones <- clones[tree_cells, ]
  
  return(clones)
}

make_corrupt_tree_heatmap <- function(tree_ggplot) {
  tree_annot_func = AnnotationFunction(
    fun=function(index) {
      pushViewport(viewport(height=1))
      grid.draw(ggplotGrob(tree_ggplot)$grobs[[5]])
      popViewport()
    },
    var_import=list(tree_ggplot=tree_ggplot),
    width=unit(4, "cm"),
    which="row"
  )
  tree_annot <- HeatmapAnnotation(
    tree=tree_annot_func, which="row", show_annotation_name=FALSE
  )
  
  n_cells <- sum(tree_ggplot$data$isTip)
  tree_hm <- Heatmap(matrix(nc=0, nr=n_cells), left_annotation=tree_annot)
  
  return(tree_hm)
}

find_largest_contiguous_group <- function(x) {
  starts <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1)
  ends <- c(starts[-1] - 1, length(x))
  largest <- which.max(ends - starts + 1)
  return(x[starts[largest]:ends[largest]])
}

get_clone_label_pos <- function(clones) {
  clone_label_pos <- list()
  for(clone in unique(clones$clone_id)) {
    if(!grepl("None", clone)) {
      clone_idx <- which(clones$clone_id == clone)
      clone_idx <- find_largest_contiguous_group(clone_idx)
      clone_label_pos[[as.character(clone)]] <-
        as.integer(round(mean(clone_idx)))
    }
  }
  return(clone_label_pos)
}

get_library_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}

get_groupings <- function(cell_ids, grouping_file) {
  library_labels <- sapply(strsplit(cell_ids, "-"), function(x) {return(x[2])})
  groupings <- read_tsv(grouping_file)
  
  tmp <- as.data.frame(library_labels, stringsAsFactors=FALSE)
  tmp <- left_join(tmp, groupings, by = "library_labels")
  return(tmp$grouping)
}

get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    # return(paste0(x[1],'-',x[2])) # some special cases with hyphen -
    return(x[1])
  })
  return(as.character(labels))
}


# library_grouping.csv should contain a column grouping or library_id, and sample
get_library_grouping <- function(cell_ids, grouping_file) {
  
  library_ids <- get_library_id(cell_ids)
  groupings <- read_tsv(grouping_file)
  
  # colnames(groupings)[which(names(groupings) == "grouping")] <- "library_id"
  sample_ids <- get_sample_id(cell_ids)
  tmp <- data.frame(cell_id=cell_ids, library_id=library_ids, sample_id=sample_ids, stringsAsFactors=FALSE)
  tmp <- left_join(tmp, groupings, by = c("library_id","sample_id"))  #composite key
  # print("DEBUG")
  print(nrow(tmp)==length(cell_ids))
  return(tmp)
}
# cell_ids <- rownames(copynumber)
get_library_grouping_v2 <- function(cell_ids, grouping_file) {
  
  ## Note: get id from cell id, ex: cell id is: AT11391-A98166B-R45-C08, 
  ## detect '-' character, the library id is: A98166B, and sample id is: AT11391 
  library_ids <- sapply(strsplit(cell_ids, "-"), function(x) {return(x[2])})
  sample_ids <- sapply(strsplit(cell_ids, "-"), function(x) {return(x[1])})
  # library_labels <- sapply(strsplit(cell_ids, "-"), function(x) {return(x[3])}) # special case
  # sample_ids <- sapply(strsplit(cell_ids, "-"), function(x) {return(paste0(x[1],'-',x[2]))}) # special case
  # groupings <- read_tsv(grouping_file)
  groupings <- data.table::fread(grouping_file) %>% as.data.frame()
  if('grouping' %in% colnames(groupings) & !'library_id' %in% colnames(groupings)){
    groupings <- groupings %>%
      dplyr::rename(library_id=grouping)
  }
  if('sample' %in% colnames(groupings) & !'sample_id' %in% colnames(groupings)){
    groupings <- groupings %>%
      dplyr::rename(sample_id=sample)
  }
  print("DEBUG")
  tmp <- data.frame(library_id=library_ids, sample_id=sample_ids, stringsAsFactors=FALSE)
  print(dim(tmp))
  tmp <- left_join(tmp, groupings, by = c("library_id","sample_id"))
  print(dim(tmp))
  return(tmp)
}

make_left_annot <- function(copynumber, clones, grouping_file) {
  annot_colours <- list()
  annot_cols <- 3
  library_labels <- get_library_id(rownames(copynumber))
  library_levels <- mixedsort(unique(library_labels))
  annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
  library_legend_rows <- 10
  print("DEBUG2")
  grouping_labels <- NULL
  annot_colours$Groupings <- NULL
  grouping_legend <- NULL
  if (!is.null(grouping_file)) {
    grouping_labels <- get_library_grouping(rownames(copynumber), grouping_file)
    print(unique(grouping_labels))
    grouping_levels <- mixedsort(unique(grouping_labels))
    annot_colours$Groupings <- make_discrete_palette("Set2", grouping_levels)
    # grouping_legend <- list(nrow=8)
    annot_cols <- 4
  }
  print("DEBUG3")
  if(!is.null(clones)) {
    clone_levels <- unique(clones$clone_label)
    clone_level_none <- clone_levels[grepl("None", clone_levels)]
    clone_levels <- mixedsort(clone_levels[!grepl("None", clone_levels)])
    
    clone_pal <- make_clone_palette(clone_levels)
    print(clone_pal)
    if(length(clone_level_none > 0)) {
      clone_pal[[clone_level_none]] <- clone_none_black
    }
    annot_colours$Clone <- clone_pal
    print("DEBUG4")
    clone_label_generator <- function(index) {
      clone_label_pos <- get_clone_label_pos(clones)
      y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
      grid.text(
        names(clone_label_pos), 0.5, y_pos,
        just=c("centre", "centre")
      )
    }
    
    clone_legend_rows <- 10
    if(length(clone_levels) > 10) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
    }
    print("DEBUG5")
    left_annot <- HeatmapAnnotation(
      Clone=clones$clone_label, clone_label=clone_label_generator,
      Sample=library_labels, Groupings=grouping_labels,
      col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
      which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
      annotation_legend_param=list(
        Clone=list(nrow=clone_legend_rows),
        Sample=list(nrow=library_legend_rows),
        Groupings=grouping_legend
      )
    )
  } else {
    left_annot <- HeatmapAnnotation(
      Sample=library_labels, col=annot_colours,
      which="row", simple_anno_size=unit(0.4, "cm"),
      annotation_legend_param=list(
        Sample=list(nrow=library_legend_rows)
      )
    )
  }
  
  return(left_annot)
}


# copynumber_backup <- copynumber
make_left_annot_v2 <- function(copynumber, clones, grouping_file) {
  annot_colours <- list()
  annot_cols <- 3
  print("DEBUG1")
  library_labels <- get_library_id(rownames(copynumber))
  library_levels <- mixedsort(unique(library_labels))
  # annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
  
  grouping_labels <- NULL
  annot_colours$Groupings <- NULL
  grouping_legend <- NULL
  if (!is.null(grouping_file)) {
    print("DEBUG2")
    # lib_group_meta <- get_library_grouping(rownames(copynumber), grouping_file)
    lib_group_meta <- get_library_grouping_v2(rownames(copynumber), grouping_file)
    # grouping_labels <- get_library_grouping(rownames(copynumber), grouping_file)
    grouping_labels <- lib_group_meta$library_id
    grouping_levels <- mixedsort(unique(lib_group_meta$library_id))
    annot_colours$Groupings <- make_discrete_palette("Set2", grouping_levels)
    grouping_legend <- list(nrow=4)
    annot_cols <- 5
    print(grouping_levels)
    print(annot_colours$Groupings)
    print("DEBUG3")
    library_legend_rows <- 4
    if(length(grouping_levels) > 10) {
      library_legend_rows <- round(sqrt(length(library_legend_rows) * 4))
    }
    
    # treatment_labels <- lib_group_meta$treatment_st
    # treatment_levels <- mixedsort(unique(treatment_labels))
    # annot_colours$TreatmentSt <- make_discrete_palette("Set1", treatment_levels)
    # print("DEBUG 2")
    # lib_group_meta <- get_library_grouping(rownames(copynumber), grouping_file)
    # lib_group_meta <- get_library_grouping_v2(rownames(copynumber), grouping_file) # for special case where library_id is not unique, use library_id and sample_id as composite key
    # mainsite_labels <- lib_group_meta$mainsite
    # mainsite_levels <- mixedsort(unique(mainsite_labels))
    # annot_colours$MainSite <- make_discrete_palette("Set1", mainsite_levels)
    # mainsite_legend_rows <- length(mainsite_levels)
    # 
    # library_legend_rows <- 10
    # print("DEBUG4")
    # celltype_labels <- lib_group_meta$celltype
    # celltype_levels <- mixedsort(unique(celltype_labels))
    # celltype_levels <- ifelse(is.na(celltype_levels),'None',celltype_levels) # in case 
    # annot_colours$celltype <- make_discrete_palette("Accent", celltype_levels)
    # print(celltype_levels)
    # print(annot_colours$celltype)
    # 
    # drug_labels <- lib_group_meta$drug
    # drug_levels <- mixedsort(unique(drug_labels))
    # drug_levels <- ifelse(is.na(drug_levels),'None',drug_levels) # in case 
    # annot_colours$drug <- make_discrete_palette("Set1", drug_levels)
    # print("DEBUG5")
    # print(drug_levels)
    # print(annot_colours$drug)
    # pdx_labels <- lib_group_meta$PDX
    # pdx_levels <- mixedsort(unique(pdx_labels))
    # print(pdx_levels)
    # # annot_colours$PDX <- make_discrete_palette("Pastel1", pdx_levels)
    # print("DEBUG5")
    # # @Sohrab: need to modify this color scheme, depend on your data, here I have SA535_cisplatin, SA535_CX, SA535_untreated
    # col_pdx <- c('#008000','#800080','#C0C0C0') # SA535, green cisplatin, purple: CX 5461, grey:untreated
    # # col_pdx <- c('#008000','#C0C0C0','#FFFF00')   # SA1035
    # names(col_pdx) <- pdx_levels
    # annot_colours$PDX <- col_pdx
  }
  # print(annot_colours)
  
  if(!is.null(clones)) {
    print("DEBUG 3")
    clone_levels <- unique(clones$clone_label)
    clone_level_none <- clone_levels[grepl("None", clone_levels)]
    clone_levels <- mixedsort(clone_levels[!grepl("None", clone_levels)])
    
    clone_pal <- make_clone_palette(clone_levels)
    print(clone_pal)
    if(length(clone_level_none > 0)) {
      clone_pal[[clone_level_none]] <- clone_none_black
    }
    annot_colours$Clone <- clone_pal
    
    clone_label_generator <- function(index) {
      clone_label_pos <- get_clone_label_pos(clones)
      y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
      grid.text(
        names(clone_label_pos), 0.5, y_pos,
        just=c("centre", "centre")
      )
    }
    
    clone_legend_rows <- 10
    if(length(clone_levels) > 10) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
    }
    # Hoa, modification
    # annot_cols <- 4  # remove lib idx
    annot_cols <- 3
    # annot_cols <- length(annot_colours) + 2
    print(names(annot_colours))
    left_annot <- HeatmapAnnotation(
      Clone=clones$clone_label, clone_label=clone_label_generator, 
      #drug=drug_labels, celltype=celltype_labels, 
      # MainSite=mainsite_labels, 
      Groupings=grouping_labels,
      col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE, TRUE),  #, TRUE
      which="row", annotation_width=unit(rep(0.7, annot_cols), "cm"),
      annotation_legend_param=list(
        Clone=list(nrow=clone_legend_rows),
        # MainSite=list(nrow=mainsite_legend_rows),
        #drug=list(nrow=library_legend_rows),
        #celltype=list(nrow=library_legend_rows),
        Groupings=list(nrow=library_legend_rows)
      )
    )
    # left_annot <- HeatmapAnnotation(
    #     Clone=clones$clone_label, clone_label=clone_label_generator,
    #     MainSite=mainsite_labels, Pdx=pdx_labels, Groupings=grouping_labels,
    #     col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
    #     which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
    #     annotation_legend_param=list(
    #         Clone=list(nrow=clone_legend_rows),
    #         MainSite=list(nrow=library_legend_rows),
    #         Pdx=list(nrow=library_legend_rows),
    #         Groupings=grouping_legend
    #     )
    # )
  } 
  # else {
  #     left_annot <- HeatmapAnnotation(
  #         TreatmentSt=mainsite_labels, col=annot_colours,
  #         which="row", simple_anno_size=unit(0.4, "cm"),
  #         annotation_legend_param=list(
  #             TreatmentSt=list(nrow=library_legend_rows)
  #         )
  #     )
  # }
  
  return(left_annot)
}

get_chrom_label_pos <- function(copynumber) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    chrom_label_pos[[chrom]] <- as.integer(round(mean(chrom_idx)))
  }
  return(chrom_label_pos)
}

# From ComplexHeatmap, needed for modified anno_mark
recycle_gp = function(gp, n = 1) {
  for(i in seq_along(gp)) {
    x = gp[[i]]
    gp[[i]] = c(rep(x, floor(n/length(x))), x[seq_len(n %% length(x))])
  }
  return(gp)
}

# From ComplexHeatmap, modified
anno_mark = function(at, labels, which = c("column", "row"),
                     side = ifelse(which == "column", "top", "right"),
                     lines_gp = gpar(), labels_gp = gpar(), padding = 0.5,
                     link_width = unit(5, "mm"), link_height = link_width,
                     link_gp = lines_gp,
                     extend = unit(0, "mm")) {
  
  which = match.arg(which)[1]
  
  if(!is.numeric(at)) {
    # stop_wrap(paste0("`at` should be numeric ", which, " index corresponding to the matrix."))
    print(paste0("`at` should be numeric ", which, " index corresponding to the matrix."))
  }
  
  n = length(at)
  link_gp = recycle_gp(link_gp, n)
  labels_gp = recycle_gp(labels_gp, n)
  labels2index = structure(seq_along(at), names = labels)
  at2labels = structure(labels, names = at)
  
  if(length(extend) == 1) extend = rep(extend, 2)
  if(length(extend) > 2) extend = extend[1:2]
  if(!inherits(extend, "unit")) extend = unit(extend, "npc")
  
  height = link_width + max_text_width(labels, gp = labels_gp)
  width = unit(1, "npc")
  
  .pos = NULL
  .scale = NULL
  
  column_fun = function(index) {
    n = length(index)
    
    # adjust at and labels
    at = intersect(index, at)
    if(length(at) == 0) {
      return(NULL)
    }
    labels = at2labels[as.character(at)]
    
    if(is.null(.scale)) {
      .scale = c(0.5, n+0.5)
    }
    pushViewport(viewport(yscale = c(0, 1), xscale = .scale))
    if(inherits(extend, "unit")) extend = convertWidth(extend, "native", valueOnly = TRUE)
    # text_height = convertWidth(grobHeight(textGrob(labels, gp = labels_gp))*(1+padding), "native", valueOnly = TRUE)
    text_height = convertWidth(grobWidth(textGrob(labels, gp = labels_gp))*(2+padding), "native", valueOnly = TRUE)
    if(is.null(.pos)) {
      i2 = which(index %in% at)
      pos = i2 # position of rows
    } else {
      pos = .pos[which(index %in% at)]
    }
    h1 = pos - text_height*0.5
    h2 = pos + text_height*0.5
    pos_adjusted = smartAlign(h1, h2, c(.scale[1] - extend[1], .scale[2] + extend[2]))
    h = (pos_adjusted[, 1] + pos_adjusted[, 2])/2
    
    n2 = length(labels)
    # grid.text(labels, h, rep(max_text_width(labels, gp = labels_gp), n2), default.units = "native", gp = labels_gp, rot = 0, just = "center")
    grid.text(labels, h, rep(grobHeight(textGrob(labels, gp = labels_gp)), n2), default.units = "native", gp = labels_gp, rot = 0, just = "center")
    link_height = link_height - unit(1, "mm")
    grid.segments(pos, unit(rep(1, n2), "npc"), pos, unit(1, "npc")-rep(link_height*(1/3), n2), default.units = "native", gp = link_gp)
    grid.segments(pos, unit(1, "npc")-rep(link_height*(1/3), n2), h, unit(1, "npc")-rep(link_height*(2/3), n2), default.units = "native", gp = link_gp)
    grid.segments(h, unit(1, "npc")-rep(link_height*(2/3), n2), h, unit(1, "npc")-rep(link_height, n2), default.units = "native", gp = link_gp)
    upViewport()
  }
  
  fun = column_fun
  
  anno = AnnotationFunction(
    fun = fun,
    fun_name = "anno_mark",
    which = which,
    width = width,
    height = height,
    n = -1,
    var_import = list(at, labels2index, at2labels, link_gp, labels_gp, padding, .pos, .scale,
                      side, link_width, link_height, extend),
    show_name = FALSE
  )
  
  # anno@subset_rule$at = subset_by_intersect
  
  anno@subsettable = TRUE
  return(anno)
}

make_bottom_annot <- function(copynumber) {
  chrom_label_pos <- get_chrom_label_pos(copynumber)
  bottom_annot <- HeatmapAnnotation(chrom_labels=anno_mark(
    at=chrom_label_pos,
    labels=names(chrom_label_pos),
    side="bottom",
    padding=0.5, extend=0.01
  ), show_annotation_name=FALSE)
  return(bottom_annot)
}

# Add origin columns
# make_left_annot_v3 <- function(copynumber, clones, grouping_file) {
#     annot_colours <- list()
#     annot_cols <- 3
#     
#     library_labels <- get_library_id(rownames(copynumber))
#     library_levels <- mixedsort(unique(library_labels))
#     # annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
#     # library_legend_rows <- 10
#     
#     grouping_labels <- NULL
#     annot_colours$Groupings <- NULL
#     grouping_legend <- NULL
#     if (!is.null(grouping_file)) {
#         grouping_labels <- get_groupings(rownames(copynumber), grouping_file)
#         grouping_levels <- mixedsort(unique(grouping_labels))
#         # annot_colours$Groupings <- make_discrete_palette("Set2", grouping_levels)
#         grouping_legend <- list(nrow=10)
#         annot_cols <- 5
#         
#         lib_group_meta <- get_library_grouping(rownames(copynumber), grouping_file)
#         mainsite_labels <- lib_group_meta$mainsite
#         mainsite_levels <- mixedsort(unique(mainsite_labels))
#         annot_colours$MainSite <- make_discrete_palette("Set1", mainsite_levels)
#         
#         origin_labels <- lib_group_meta$origin
#         origin_levels <- mixedsort(unique(origin_labels))
#         annot_colours$Origin <- make_discrete_palette("Pastel1", origin_levels)
#         
#         library_legend_rows <- 10
#         pdx_labels <- lib_group_meta$pdxid
#         pdx_levels <- mixedsort(unique(pdx_labels))
#         annot_colours$Pdx <- make_discrete_palette("Accent", pdx_levels)
#         # library_legend_rows <- 10
#     }
#     
#     if(!is.null(clones)) {
#         clone_levels <- unique(clones$clone_label)
#         clone_level_none <- clone_levels[grepl("None", clone_levels)]
#         clone_levels <- mixedsort(clone_levels[!grepl("None", clone_levels)])
#         
#         clone_pal <- make_clone_palette(clone_levels)
#         if(length(clone_level_none > 0)) {
#             clone_pal[[clone_level_none]] <- clone_none_black
#         }
#         annot_colours$Clone <- clone_pal
#         
#         clone_label_generator <- function(index) {
#             clone_label_pos <- get_clone_label_pos(clones)
#             y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
#             grid.text(
#                 names(clone_label_pos), 0.5, y_pos,
#                 just=c("centre", "centre")
#             )
#         }
#         
#         clone_legend_rows <- 10
#         if(length(clone_levels) > 10) {
#             clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
#         }
#         # Hoa, modification
#         # annot_cols <- 4  # remove lib idx
#         annot_cols <- 5
#         left_annot <- HeatmapAnnotation(
#             Clone=clones$clone_label, clone_label=clone_label_generator,
#             MainSite=mainsite_labels, Origin=origin_labels, Pdx=pdx_labels,
#             col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE),
#             which="row", annotation_width=unit(rep(0.5, annot_cols), "cm"),
#             annotation_legend_param=list(
#                 Clone=list(nrow=clone_legend_rows),
#                 MainSite=list(nrow=library_legend_rows),
#                 Origin=list(nrow=library_legend_rows),
#                 Pdx=list(nrow=library_legend_rows)
#             )
#         )
#         # left_annot <- HeatmapAnnotation(
#         #     Clone=clones$clone_label, clone_label=clone_label_generator,
#         #     MainSite=mainsite_labels, Pdx=pdx_labels, Groupings=grouping_labels,
#         #     col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
#         #     which="row", annotation_width=unit(rep(0.4, annot_cols), "cm"),
#         #     annotation_legend_param=list(
#         #         Clone=list(nrow=clone_legend_rows),
#         #         MainSite=list(nrow=library_legend_rows),
#         #         Pdx=list(nrow=library_legend_rows),
#         #         Groupings=grouping_legend
#         #     )
#         # )
#     } else {
#         left_annot <- HeatmapAnnotation(
#             MainSite=mainsite_labels, col=annot_colours,
#             which="row", simple_anno_size=unit(0.4, "cm"),
#             annotation_legend_param=list(
#                 MainSite=list(nrow=library_legend_rows)
#             )
#         )
#     }
#     
#     return(left_annot)
# }

make_copynumber_heatmap <- function(copynumber, clones, grouping_file) {
  copynumber_hm <- Heatmap(
    name="Copy Number",
    as.matrix(copynumber),
    col=cn_colours,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    bottom_annotation=make_bottom_annot(copynumber),
    left_annotation=make_left_annot_v2(copynumber, clones, grouping_file),  #Hoa modified
    # heatmap_legend_param=list(nrow=10),
    heatmap_legend_param=list(nrow=4, by_row=FALSE),
    use_raster=TRUE,
    raster_quality=5
  )
  return(copynumber_hm)
}
# brlen <- NULL
make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
                                              brlen, grouping_file) {
  print('format tree')
  tree <- format_tree(tree, brlen)
  print('tree_ggplot')
  tree_ggplot <- make_tree_ggplot(tree, clones)
  print('format_copynumber')
  copynumber <- format_copynumber(copynumber, tree_ggplot$data)
  # copynumber <- format_copynumber(copynumber)
  if(!is.null(clones)) {
    # clones <- format_clones(clones)
    print('format_clones')
    clones <- format_clones(clones, tree_ggplot$data)
    # copynumber <- copynumber[,clones$cell_id]
  }
  print('tree_hm')
  tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
  print('copynumber_hm')
  copynumber_hm <- make_copynumber_heatmap(copynumber, clones, grouping_file)
  print('visualize tree')
  h <- tree_hm + copynumber_hm
  # h <- copynumber_hm
  print('Draw tree')
  draw(
    h,
    padding=unit(c(10, 2, 2, 2), "mm"),
    annotation_legend_side="right",
    heatmap_legend_side="bottom"
  )
}

main <- function() {
  argv <- get_args()
  
  tree <- read.tree(argv$tree)
  print(argv$tree)
  print(argv$copynumber)
  copynumber <- read.csv(argv$copynumber, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  # copynumber <- read_tsv(argv$copynumber)
  
  if(argv$normalize_ploidy) {
    print(paste(
      "WARNING: ceiling() will be applied to non-integer copynumber",
      "during normalization"
    ))
    copynumber <- normalize_cell_ploidy(copynumber)
  }
  
  if(!is.null(argv$clones)) {
    clones <- read_tsv(argv$clones)
  } else {
    clones <- NULL
  }
  
  # pdf(argv$pdf, width=10)
  png(argv$pdf, height = 2*800, width=2*1400, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, argv$branch_lengths, argv$grouping_file
  )
  dev.off()
}


main_debug <- function() {
  
  tree <- read.tree(tree_fn)
  
  copynumber <- read.csv(copynumber_fn, header=T, row.names=1, check.names = F,stringsAsFactors = FALSE)
  # copynumber <- read_tsv(argv$copynumber)
  
  if(argv$normalize_ploidy) {
    print(paste(
      "WARNING: ceiling() will be applied to non-integer copynumber",
      "during normalization"
    ))
    copynumber <- normalize_cell_ploidy(copynumber)
  }
  
  if(!is.null(argv$clones)) {
    clones <- read_tsv(clones_fn)
  } else {
    clones <- NULL
  }
  
  # pdf(argv$pdf, width=10)
  png(output_fn, height = 2*800, width=2*1400, res = 2*72)
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  dev.off()
  
  make_cell_copynumber_tree_heatmap(
    tree, copynumber, clones, NULL, grouping_file
  )
  
}

## START hdbscan_clustering
download_dir <- '~/dlp_hdbscan_clustering/hdbscan_clustering/testing_data/'
results_dir <- '~/dlp_hdbscan_clustering/hdbscan_clustering/testing_data/'
save_dir <- paste0(results_dir,'dlp_results/')

libs <- list.files(paste0(results_dir))
libs <- gsub('_filtered_states.csv.gz','',libs)


length(libs)

## Note: 
# Library id: A144171A, sample id: AT23026 (not id A144171B as original one) due to pooled samples
# Library id: A144171A, sample id: AT23025
# → Building 2 heatmaps for 2 different sample ids
libs <- libs[!libs %in% c('A144171A-AT23025','A144171A-AT23026')] 



# TO DO: add filtering outliers function here, filtered by medianCNV
# datatag <- 'A138979A'

## For this library, can we get sub set of cells in a clone, and dividing cells into 2 clusters by repeat hdbscan
datatag <- 'A98166B'

## To Do: check clustering grain level for this library: 'A138850A'
datatag <- 'A98166B' # Done
libs <- libs[!libs %in% c('A138848B','A138850A')] 
libs
for(datatag in libs){
  print(datatag)
  grouping_file <- paste0(results_dir, 'library_grouping.csv') # meta samples
  # reads_fn <- paste0(raw_data_dir,datatag,'/hmmcopy/',datatag,'_reads.csv.gz') # input data from hmmcopy output
  filtered_CNV_fn <- paste0(results_dir, 'A98166B_filtered_states.csv.gz')
  
  l <- datatag
  lib_fn1 <- paste0(download_dir,l, '_reads.csv.gz') # csv.gz or csv
  # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
  lib_fn2 <- paste0(download_dir,l,'_reads.csv.gz') # csv.gz or csv
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
  grouping_file <- paste0(results_dir, 'library_grouping.csv') # meta samples
  # reads_fn <- paste0(raw_data_dir,datatag,'/hmmcopy/',datatag,'_reads.csv.gz') # input data from hmmcopy output
  filtered_CNV_fn <- paste0(results_dir,'A138848B_filtered_states.csv.gz')
  
  
  lib_fn1 <- paste0(download_dir,l,'reads.csv.gz') # csv.gz or csv
  # lib_fn <- paste0(download_dir,l,'/',l,'/hmmcopy/','reads.csv.gz') # csv.gz or csv
  lib_fn2 <- paste0(download_dir,l,'_reads.csv.gz') # csv.gz or csv
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
                 pctcells = 0.05)
  }  
}
# main()
