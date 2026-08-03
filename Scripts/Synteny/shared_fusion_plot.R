library(scales)
library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

# Defaults reproduce the filenames used in the original analysis workspace.
input_file <- if (length(args) >= 1) args[[1]] else "shared.fusion.Annelida.txt"
out_pdf <- if (length(args) >= 2) args[[2]] else "ordered_intersection_shared_fusions.Annelida.pdf"
out_csv <- if (length(args) >= 3) args[[3]] else "ordered_intersection_shared_fusions_matrix.Annelida.csv"

lines <- readLines(input_file, warn = FALSE)
lines <- lines[nzchar(trimws(lines))]

parse_fusion <- function(x) {
  vals <- suppressWarnings(as.integer(strsplit(x, "\\.")[[1]]))
  vals <- vals[!is.na(vals)]
  unique(vals)
}

species <- character()
fusion_list <- list()

for (ln in lines) {
  parts <- strsplit(ln, "\t")[[1]]
  sp <- parts[1]
  fusions <- parts[-1]
  fusions <- fusions[nzchar(fusions)]
  parsed <- list()
  if (length(fusions) > 0) {
    for (fx in fusions) {
      chroms <- parse_fusion(fx)
      if (length(chroms) >= 2) {
        parsed[[length(parsed) + 1]] <- chroms
      }
    }
  }
  species <- c(species, sp)
  fusion_list[[length(fusion_list) + 1]] <- parsed
}

n <- length(species)
mat <- matrix(0, nrow = n, ncol = n, dimnames = list(species, species))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    shared_keys <- character()
    fi_list <- fusion_list[[i]]
    fj_list <- fusion_list[[j]]
    if (length(fi_list) > 0 && length(fj_list) > 0) {
      for (fi in fi_list) {
        for (fj in fj_list) {
          inter <- sort(intersect(fi, fj))
          if (length(inter) >= 2) {
            shared_keys <- c(shared_keys, paste(inter, collapse = "."))
          }
        }
      }
    }
    mat[i, j] <- length(unique(shared_keys))
  }
}

write.csv(mat, out_csv, quote = FALSE)

df <- reshape2::melt(mat)
colnames(df) <- c("Species1", "Species2", "SharedCount")
df$Species1 <- factor(df$Species1, levels = rev(species))
df$Species2 <- factor(df$Species2, levels = species)

p <- ggplot(df, aes(x = Species2, y = Species1, fill = SharedCount)) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  scale_fill_distiller(
    palette = "YlOrBr",
    direction = 1,
    limits = c(0, 5),          # Sets the scale range from 0 to 5
    oob = squish               # Squishes values > 5 to the max color (requires 'scales' package)
  ) +
  labs(
    title = "Species-by-species shared fusions (R)",
    subtitle = "Original file order; a pair is counted when two fusions share >=2 chromosomes",
    x = "Species (original order)",
    y = "Species (original order)",
    fill = "Shared\ncount"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic", size = 12),
    axis.text.y = element_text(face = "italic", size = 12),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9)
  ) +
  coord_fixed()

ggsave(out_pdf, p, width = 24, height = 28, dpi = 300)
