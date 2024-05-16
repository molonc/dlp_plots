mixedsort <- gtools::mixedsort

CNV_COLOURS <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
  ),
  names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
)

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

get_ordered_cell_ids <- function(tree_plot_dat) {
  return(rev(arrange(tree_plot_dat[tree_plot_dat$isTip, ], y)$label))
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

find_largest_contiguous_group <- function(x) {
  starts <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1)
  ends <- c(starts[-1] - 1, length(x))
  largest <- which.max(ends - starts + 1)
  return(x[starts[largest]:ends[largest]])
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


format_copynumber <- function(copynumber, tree_plot_dat = NULL, spacer_cols=20) {
  if (!("chr" %in% colnames(copynumber))) {
    loci <- sapply(rownames(copynumber), strsplit, "_")
    copynumber$chr <- unname(sapply(loci, '[[', 1))
    copynumber$start <- as.numeric(unname(sapply(loci, '[[', 2)))
    copynumber$end <- as.numeric(unname(sapply(loci, '[[', 3)))
    copynumber$width <- (copynumber$end - copynumber$start + 1)
  }
  copynumber$chr <- gsub("chr", "", copynumber$chr)
  copynumber <- arrange(copynumber, as.numeric(chr), chr, start)
  
  rownames(copynumber) <- paste0(
    copynumber$chr, ":", copynumber$start, ":", copynumber$end
  )
  copynumber <- subset(copynumber, select=-c(chr, start, end, width))
  copynumber <- as.data.frame(t(copynumber))
  
  if (!is.null(tree_plot_dat)) {
  	  copynumber <- copynumber[get_ordered_cell_ids(tree_plot_dat), ]
  }
  
  copynumber <- format_copynumber_values(copynumber)
  copynumber <- space_copynumber_columns(copynumber, spacer_cols)
  
  return(copynumber)
}

dlp_long_to_wide <- function(raw) {
	raw$pos <- with(raw, paste0(chr, "_", start, "_", end))
	long <- raw[, c("pos", "state", "cell_id")]
	wide <- tidyr::pivot_wider(long, values_from = state, names_from = cell_id)
	wide <- tibble::column_to_rownames(wide, "pos")
	return(wide)
}

make_clone_palette <- function(levels) {
    clone_names <- sort(levels)
    n <- length(clone_names)
    if (n < 3) {
        pal <- brewer.pal(3, "Set1")
        pal <- head(pal, n)
    } else {
        pal <- brewer.pal(n, "Set1")
        if (n > 9) {
            pal <- brewer.pal(9, "Set1")
            pal <- c(pal, brewer.pal(n - 9, "Set2"))
            if (n > 17) {
                pal <- rainbow(n)
            }
        }
    }

    pal <- as.character(pal)
    names(pal) <- clone_names

    pal <- pal[levels]
    return(pal)
}


make_annot_chr_labels <- function(cnv_matrix) {
  output <- c()
  chroms <- sapply(strsplit(colnames(cnv_matrix), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    output[[chrom]] <- as.integer(round(min(chrom_idx)))
  }
  return(output)
}

CLONE_LABEL_GENERATOR <- function(index) {
  clone_label_pos <- get_clone_label_pos(clones)
  y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
  grid.text(
    names(clone_label_pos), 0.5, y_pos,
    just=c("centre", "centre")
  )
}