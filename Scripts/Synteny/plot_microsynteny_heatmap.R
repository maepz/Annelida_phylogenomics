#!/usr/bin/env Rscript

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("The R package 'pheatmap' is required.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  "poisson_microsynteny_distance_matrix.tsv"
}
output_pdf <- if (length(args) >= 2) {
  args[[2]]
} else {
  "microsynteny_heatmap.pdf"
}

distance_table <- read.delim(
  input_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
distance_matrix <- as.matrix(distance_table)
storage.mode(distance_matrix) <- "double"

if (!identical(rownames(distance_matrix), colnames(distance_matrix))) {
  stop("Row and column taxon labels do not match.", call. = FALSE)
}
if (!isTRUE(all.equal(distance_matrix, t(distance_matrix)))) {
  stop("The distance matrix is not symmetric.", call. = FALSE)
}

upgma <- hclust(as.dist(distance_matrix), method = "average")
color_palette <- colorRampPalette(c("#0B0405", "#357BA2", "#DEF5E5"))(100)

pdf(output_pdf, width = 10, height = 8)
pheatmap::pheatmap(
  distance_matrix,
  cluster_rows = upgma,
  cluster_cols = upgma,
  color = color_palette,
  border_color = "white",
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Poisson-corrected synteny decay",
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)
invisible(dev.off())
