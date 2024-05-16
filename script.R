library(ComplexHeatmap) # Bioconductor
library(vroom)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(tibble)
library(ggtree) # Bioconductor
library(stringr)
library(signals) #devtools::install_github("shahcompbio/signals")

source("~/Desktop/hdbscan_minimal/hdbscan_minimal.R")

blacklist <- read.delim("~/Desktop/hdbscan_minimal/blacklist_2018.10.23.txt", sep = " ")

# metrics_files <- Sys.glob("../data/*/annotation/*metrics.csv.gz")
metrics_files <- "/Users/jeschan/Desktop/hdb_scan/download/SC-8888/results/annotation/metrics.csv.gz"
# metrics_files <- 
names(metrics_files) <- str_extract(metrics_files, "SC-[0-9]+")

metrics_raw <- lapply(metrics_files, read.csv)
metrics_raw <- mapply(add_column, .data = metrics_raw, file = metrics_files, SIMPLIFY = FALSE)
metrics <- bind_rows(metrics_raw)

# reads_files <- Sys.glob("../data/*/*/*reads.csv.gz")
reads_files <- "/Users/jeschan/Desktop/hdb_scan/download/SC-8888/results/hmmcopy/reads.csv.gz"
names(reads_files) <- str_extract(reads_files, "SC-[0-9]+")

reads_raw <- lapply(reads_files, vroom)
reads_raw <- mapply(add_column, .data = reads_raw, file = reads_files, SIMPLIFY = FALSE)
reads <- bind_rows(reads_raw)


# TODO 1: ######### FILTER HERE ######### You could even subset set, sample, etc.
keep <- subset(metrics, quality > 0.75 & experimental_condition %in% c("A", "B", "C", "D"))
# hTERT <- subset(metrics, grepl("hTERT", metrics$experimental_condition) & quality > 0.75 & state_mode == 2)
# keep <- bind_rows(keep, hTERT)

# TODO END

# FIGURE 1 - Cells present
ggplot(metrics, aes(experimental_condition)) + geom_bar()

# FIGURE 2 - Sequencing results
ggplot(metrics, aes(total_reads, total_mapped_reads, col = experimental_condition)) + geom_point() + coord_fixed()
ggplot(metrics, aes(y = total_reads, x = experimental_condition)) + geom_boxplot()

select <- subset(reads, cell_id %in% keep$cell_id)

# Reads matrix
reads_df <- select
reads_df$chr_desc <- paste(reads_df$chr, reads_df$start, reads_df$end, sep = "_")

# CNV matrix
copynumber <- dlp_long_to_wide(reads_df)

bl <- read.delim("~/Desktop/hdbscan_minimal/blacklist_2018.10.23.txt", sep = " ")
bl$bins <- bl$width / 5e5
blbig <- data.frame()
for (i in 1:nrow(bl)) {
	row <- bl[i, ]
	blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
	blbig <- bind_rows(blbig, blbit)
}
blbig$blacklist <- TRUE
tmp <- subset(reads_df, cell_id == head(reads_df$cell_id, 1))
tmp <- merge(tmp, blbig, all = TRUE)

exclude <- subset(tmp, blacklist)$chr_desc
copynumber[rownames(copynumber) %in% exclude, ] <- NA

reads_df <- reads_df %>%
	dplyr::filter(cell_id %in% colnames(copynumber) & 
	chr_desc %in% rownames(copynumber))

# minimum number of cells to consider as a clone
# pctcells = 0.01
# ncells <- length(unique(reads_df$cell_id))
# clone_min_cells <- max(round(pctcells * ncells), 2) 

# TODO 2: DEFINE YOUR CLONE SIZE
clone_min_cells <- 100

# TODO 3: SET YOUR SEED
set.seed(100)
print(i)
#?umap_clustering
clusters <- signals::umap_clustering(reads_df, 
	minPts = clone_min_cells, 
	field = "copy")
#### END TODO


orig_tree <- clusters$tree
tree <- orig_tree
clones <- clusters$clustering
print(unique(clones$clone_id))

brlen <- NULL

# make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
#                                               brlen, grouping_file) {

# FORMAT TREE, cleans it up, removes 
# tree <- format_tree(tree, brlen)
# locus_tips <- grep('locus', tree$tip.label, value=TRUE)
# tree <- drop.tip(tree, locus_tips)

# if (!is.null(brlen)) {
# 	tree <- compute.brlen(tree, brlen)
# }

# print('tree_ggplot')
# tree_ggplot <- make_tree_ggplot(tree, clones)
tmp <- group_by(clones, clone_id) %>% summarise(clone_members = list(cell_id))
names(tmp$clone_members) <- tmp$clone_id
clone_members <- tmp$clone_members

tree <- ggtree::groupOTU(tree, clone_members)

clone_levels <- mixedsort(unique(clones$clone_id))
clone_pal <- make_clone_palette(clone_levels)
clone_pal["0"] <- "#1B1B1B"

tree_aes <- aes(x, y, colour=group)

p <- ggplot(tree, tree_aes) +
	geom_tree(size=0.25) +
	coord_cartesian(expand=FALSE) +
	ylim(0.5, length(tree$tip.label) + 0.5) +
	theme_void() + scale_colour_manual(values=clone_pal)

tree_ggplot <- p

old_clones <- clones

# TREE ENDS HERE

orig_copynumber <- copynumber

copynumber <- format_copynumber(copynumber, tree_ggplot$data)
clones <- format_clones(clones, tree_ggplot$data) # ADDS CLONE_LABEL HERE
tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)

# TEST THIS
chrom_label_pos <- c()
chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
uniq_chroms <- c(as.character(1:22), "X", "Y")
for(chrom in uniq_chroms) {
	chrom_idx <- which(chroms == chrom)
	chrom_label_pos[[chrom]] <- as.integer(round(min(chrom_idx)))
}

bottom_annot <- HeatmapAnnotation(
	chrom_labels =
	anno_mark(
		at = unlist(chrom_label_pos),
		labels = names(chrom_label_pos),
		side = "bottom",
		padding = 0.5, extend = 0.01
	),
	show_annotation_name=FALSE)

# LEFT ANNOT
# TODO 4 customize annotation/legend if required (DEFAULT SHOULD BE FAIRLY SAFE MOST OF THE TIME)
annot_colours <- list()
# annot_cols <- 3

library_labels <- stringr::str_extract(rownames(copynumber), "A[0-9]+[A-D]?")
library_levels <- mixedsort(unique(library_labels))

tmp_exp <- subset(select(metrics, cell_id, experimental_condition, quality, state_mode), cell_id %in% clones$cell_id)
clones <- plyr::join(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
# clones3 <- merge(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
exp_levels <- mixedsort(unique(clones$experimental_condition))

n <- length(library_levels) + length(exp_levels)
six_col <- brewer.pal(n, "Set2")

lib_col <- six_col[1:length(library_levels)]
names(lib_col) <- library_levels
annot_colours$Library <- lib_col

clone_levels <- mixedsort(unique(clones$clone_label))
clo_col <- make_clone_palette(clone_levels)
annot_colours$Clone <- clo_col

exp_col <- six_col[(length(library_levels) + 1):(length(library_levels) + length(exp_levels))]
names(exp_col) <- exp_levels
annot_colours$Batch <- exp_col

left_annot <- HeatmapAnnotation(
	Clone=clones$clone_label, 
	clone_label=CLONE_LABEL_GENERATOR,
	Batch = clones$experimental_condition,
	# fixation=fix_labels,
	# sample=sample_labels, 
	# drug=drug_labels, 
	# celltype=celltype_labels, 
	Library=library_labels,
	col=annot_colours, 
	# show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE,TRUE),  #, TRUE
	show_annotation_name = c(TRUE, FALSE, TRUE),
	which="row", 
	annotation_width=unit(rep(0.5, 4), "cm"), # MATCH NUMBER OF ANNOTS
	annotation_legend_param=list(
		# Clone=list(nrow=10),
		# fixation=list(nrow=library_legend_rows),
		# sample=list(nrow=library_legend_rows),
		# drug=list(nrow=library_legend_rows),
		# celltype=list(nrow=library_legend_rows),
		# Groupings=list(nrow=10)
	)
)

# TODO 4 END HERE
# END LEFT ANNOT

copynumber_hm <- Heatmap(
	name="Copy Number",
	as.matrix(copynumber),
	col=CNV_COLOURS,
	na_col="white",
	show_row_names=FALSE,
	cluster_rows=FALSE,
	cluster_columns=FALSE,
	show_column_names=FALSE,
	bottom_annotation=bottom_annot,
	left_annotation= left_annot,  #Hoa modified
	heatmap_legend_param=list(nrow=4),
	use_raster=TRUE,
	raster_quality=5
)

print('visualize tree')
h <- tree_hm + copynumber_hm
# h <- copynumber_hm
print('Draw tree')


# TODO 5, give the proper title name + width (inches) + length (inches)
# You can also output PNG instead of PDF

# print(h)
# dev.off()
png("hdbscan_clustering_all_cells.png", 1600, 1600)
print(h)
dev.off()

# ggsave("hdbscan_clustering_all_cells.png", width=20, height=20)
pdf("hdbscan_clustering_all_cells.pdf", plot=h, 20, 20)
print(h)
dev.off()

# HEATMAP DONE!

clones$library <- str_extract(clones$cell_id, "A[^-T]+")

group_by(clones, clone_id) %>% summarise(total = length(cell_id), percentage = total / nrow(clones))

pick <- subset(clones, clone_id == "0")$cell_id
mpick <- subset(metrics, cell_id %in% pick)
slice <- subset(reads, cell_id %in% pick)

slice <- merge(slice, blbig, all.x = TRUE)
slice$state <- ifelse(is.na(slice$blacklist), slice$state, NA)

slice$chr <- factor(slice$chr, c(1:22, "X", "Y"))
slice$cell_id <- factor(slice$cell_id, mpick$cell_id[order(mpick$mean_copy)])

cols <- c(rev(brewer.pal(n = 3, "Blues"))[1:2], "#CCCCCC", tail(brewer.pal(n = 8, "OrRd"), 6))
cols <- c(cols, cols[c(9, 9, 9)])
names(cols) <- 1:12 - 1

g <- ggplot(slice, aes(start, cell_id, fill = as.factor(state))) + geom_tile() + scale_fill_manual("Copy", values = cols) + facet_grid(~chr, scales = "free", space = "free", switch = "x") + scale_x_continuous("Chromosome", expand = c(0, 0), breaks = NULL) + theme(panel.spacing = unit(0.1, "lines")) + scale_y_discrete("")

pdf("dlp_heatmap.pdf", width = 12, height = 3.1)
plot(g)
dev.off()

dir.create("png")
file.copy(paste0("/Users/dalai/research/BIOF-459/data/SC-8512/annotation/segs_pass/", pick, "_segments.png"), "png")

g2 <- ggplot(subset(slice, cell_id %in% c("AT22536-A138850A-R12-C38", "AT22536-A138850A-R09-C48", "AT22536-A138850A-R16-C64")), aes(start, copy, col = as.factor(state))) + geom_point(size = 0.5) + facet_grid(cell_id~chr, scales = "free", space = "free", switch = "x") + scale_x_continuous("Chromosome", expand = c(0, 0), breaks = NULL) + theme(panel.spacing = unit(0.1, "lines")) + scale_colour_manual("Copy", values = cols) + scale_y_continuous(breaks = seq(0, 10, by = 2), "Copy", lim = c(0, 9))

pdf("dlp_profiles.pdf", width = 10, height = 6)
plot(g2)
dev.off()
