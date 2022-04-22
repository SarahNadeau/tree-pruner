#!/usr/bin/env Rscript

# This script is to get a small, representative genome set from a large phylogenetic tree

suppressMessages(library(ape))
suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(phytools))
suppressMessages(library(tidytree))

# Take command line input
parser <- argparse::ArgumentParser()
parser$add_argument("-t", "--test", action="store_true", default=F,
                    help="do a test run")
parser$add_argument("-f", "--treefile", type="character",
                    help="tree to choose representative sequences from")
parser$add_argument("-n", "--nseqs", type = "integer", default=100,
                    help = "number of sequences to select [default %(default)s]")
parser$add_argument("-o", "--outdir", type = "character", default="results",
                    help = "directory for results [default \"%(default)s\"]")
args <- parser$parse_args()

if (args$test) {
  treefile <- "assets/test.tree"
  nseqs <- 4
} else {
  treefile <- args$treefile
  if (is.null(treefile)) {
    stop(paste("Invalid or missing argument --treefile", treefile))
  }
  nseqs <- args$nseqs
}
outdir <- args$outdir

##################################################
# Input/Output
##################################################

# Output directory
split_treefile <- strsplit(treefile, split = "\\.")[[1]]
treefile_prefix <- paste0(split_treefile[1:(length(split_treefile) - 1)], collapse = ".")  # use input tree prefix to name results
outdir2 <- paste(outdir, treefile_prefix, sep = "/")
clades_outdir <- paste0(outdir2, "/clades")
dir.create(path = clades_outdir, showWarnings = F, recursive = T)

# Load tree
tree <- ape::read.tree(treefile)
n_tips <- length(tree$tip.label)
message(paste("Loaded tree with", n_tips, "tips"))

# Check number of sequecnes to sample
if (nseqs >= n_tips) {
  stop(paste(
    "Number of sequences to sample", nseqs,
    "must be less than number of tips", n_tips,
    "(set with --nseqs)"
  ))
}

##################################################
# Cut tree to get desired number of clades
##################################################

tree_rooted <- phytools::midpoint.root(tree = tree)  # root at midpoint of longest branch
tree_ultrametric <- phytools::force.ultrametric(  # extend branches to match last sample
  tree = tree_rooted,
  method = "extend"
)
tree_binary <- ape::multi2di(phy = tree_ultrametric, random = T)  # randomly resolve polytomies
dendrogram <- ape::as.hclust.phylo(tree_binary)  # make into dendrogram

# Cut tree to get X clusters, where X = reduced sample size
clades <- cutree(tree = dendrogram, k = nseqs)
message(paste("Cut tree to form", nseqs, "clades"))

##################################################
# Get heights of each node in tree
##################################################

#' Calculates the height of a node (sum of edge lengths from root to the node)
#' @param tree_data tbl_tree object produced by tidytree's as_tibble function
#' @param node node number for node to get height of
#' @return the height of node
get_height <- function(tree_data, node, height = 0) {
  parent <- tree_data[[which(tree_data$node == node), "parent"]]
  # base case: reached root, return height
  if (parent == node) {
    return(height)  # is root
  }
  # recursive case: calculate height to parent node
  edge <- tree_data[[which(tree_data$node == node), "branch.length"]]
  height <- get_height(tree_data, parent, height + edge)
}

tbl_tree <- tidytree::as_tibble(tree)
tbl_tree$height <- unlist(lapply(
  X = tbl_tree$node,
  FUN = get_height,
  tree_data = tbl_tree
))
message("Calculated heights for each node in tree")

##################################################
# Take most basal seq from each clade (presumably closest to ancestral seq)
##################################################

# Merge clade information to tree data
clade_data <- data.frame(
  label = names(clades),
  clade = clades
)
tree_data <- tbl_tree %>% left_join(clade_data, by = "label")

# Select representative samples
tree_data_selected <- tree_data %>%
        group_by(clade) %>%
        top_n(n = 1, wt = -height) %>%
        group_by(clade) %>%
        top_n(n = 1, wt = label) %>%  # break ties by sample name
        filter(parent != node) %>%  # remove root node NA "clade"
        mutate(sampled = clade)

tree_data_full <- tree_data %>% left_join(
  tree_data_selected,
  by = c("parent", "node", "branch.length", "label", "height", "clade")
)
message("Selected most ancestral sequence from each clade")

##################################################
# Write out selected samples and the clades they represent
##################################################

f1 <- file(paste0(outdir2, "/representative_seqs.txt"), "w")

for (clade_idx in unique(tree_data_full$clade)) {
  if (is.na(clade_idx)) {
    next
  }
  clade_data <- tree_data_full %>% filter(clade == clade_idx, node <= n_tips)
  clade_labels <- clade_data$label
  clade_sample <- clade_data$label[which(!is.na(clade_data$sampled))]
  f2 <- file(paste0(clades_outdir, "/", clade_sample, ".txt"), 'w')
  writeLines(text = clade_labels, con = f2)
  close(f2)
  writeLines(text = clade_sample, con = f1)
}
close(f1)
message("Wrote out selected sequences and clades they represent")

##################################################
# Plot clade supports, get edge colors
# This is based on FastTree's SH-like local supports, which test node topology compared to 3 alternate NNIs at that node.
# Note: Missing SH-like local support values are from nodes without 3 valid NNI moves (a node with 2 descendent tips, both with branch length 0 from node).
##################################################

# Get MRCA of each clade and get edge colors for plotting
edge_colors <- rep(0, Nedge(tree))  # no color (color #0) for now, will adjust below
for (clade_idx in unique(tree_data_full$clade)) {
  if (is.na(clade_idx)) {
    next
  }
  clade_data <- tree_data_full %>% filter(clade == clade_idx)
  tips <- clade_data$node
  if (length(tips) == 1) {
    mrca <- tips[1]
  } else {
    mrca <- ape::getMRCA(phy = tree, tip = tips)
  }
  i <- which.edge(tree, tips)
  edge_colors[i] <- clade_idx
  tree_data_full$clade[which(tree_data_full$node == mrca)] <- clade_idx
}
# Get clade MRCA node support values
if (!args$test) {
  mrca_data <- tree_data_full %>% filter(node > n_tips, !is.na(clade))
  mrca_supports <- mrca_data %>%
          select(clade, label) %>%
          mutate(label = as.numeric(label)) %>%
          arrange(label) %>%
          rename("SH-like local support" = label)
  png(filename = paste0(outdir2, "/clade_support_values.png"))
  hist(mrca_supports$`SH-like local support`, breaks = 10)
  tmp <- dev.off()  # storing return value makes this quiet
}

##################################################
# Plot clade summary
##################################################

clade_sizes <- tree_data_full %>%
  filter(node <= n_tips) %>%
  group_by(clade) %>%
  summarise(n_samples = n())
png(filename = paste0(outdir2, "/clade_sizes.png"))
hist(clade_sizes$n_samples, breaks = 20)
tmp <- dev.off()

##################################################
# Plot clades and sampled representative sequences on tree
##################################################

longest_tip_label <- max(nchar(tree$tip.label))  # to scale figure width by
label_cex <- 0.25
png(
  filename = paste0(outdir2, "/tree_with_clades.png"),
  height = n_tips / 4, width = 2 + longest_tip_label * 0.25, units = "cm",
  res = 300
)
par(mar = rep(0, 4))  # remove figure margins

tip_colors <- tree_data_full$sampled %% 7 + 2  # avoid black (color #1)
edge_colors_new <- case_when(
  edge_colors == 0 ~ 1,  # make edges above clades black
  T ~ edge_colors %% 7 + 2  # avoid black in clades
)
plot(
  tree,
  tip.color = tip_colors,
  cex = label_cex,
  edge.color = edge_colors_new,
  show.node.label = T
)
tmp <- dev.off()
message("Generated results figures")
