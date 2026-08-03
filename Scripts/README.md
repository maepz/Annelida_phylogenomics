# Analysis scripts

This directory contains representative commands and cleaned custom scripts
used in the study. They are provided for methodological transparency rather
than as a turnkey reproducibility workflow. The analyses were run in a
server-side workspace whose directory structure and filenames differ from the
organization of this repository.

`Analysis_commands.Rmd` is the concise analysis-command archive, and
`Analysis_commands.pdf` is its rendered version. Complete custom programs are
kept as standalone files so that the command archive remains readable.

## Standalone scripts

| Directory | File | Purpose |
|---|---|---|
| `ILS_introgression/` | `QuIBL_wrapper.ipynb` | Representative QuIBL batch analysis for one deep branch |
| `ILS_introgression/` | `Dstat_from_alignment.ipynb` | D and f4 statistics from an amino-acid alignment, using L42 as the representative matrix |
| `Hypothesis_testing/` | `constrained_topology_tests.sh` | IQ-TREE SH and AU tests for the retained A, M, C, and S topology groups |
| `Hypothesis_testing/` | `DeltaL_analysis.ipynb` | Site-wise likelihood-deficit, opinionated-site, and gene-enrichment analyses for Topo A1–A4 |
| `Hypothesis_testing/` | `Parse_QuartetLikelihoodMapping_results.ipynb` | Parsing and plotting the Oweniidae and Sipuncula quartet-likelihood-mapping results |
| `Divergence_times/` | `MCMCTree_postprocessing.py` | Chain diagnostics and prior and posterior summaries for the MCMCTree analyses |
| `Synteny/` | `identify_microsynteny.py` | Identification of pairwise microsynteny blocks from BED and mutual-best-hit files |
| `Synteny/` | `poisson_distance_based_on_microsynteny.py` | Poisson-corrected distance matrix from retained pairwise microsynteny blocks |
| `Synteny/` | `plot_microsynteny_heatmap.R` | UPGMA clustering and heatmap of the Poisson-corrected microsynteny distances |
| `Synteny/` | `shared_fusion_plot.R` | Shared chromosomal-fusion matrix and heatmap |

Personal workstation and server paths have been removed. Scripts that access
the original analysis workspace use the `ANALYSIS_ROOT` environment variable
or an equivalent command-line option. Default placeholder paths must be
adapted before use.

The custom Python workflows use packages including pandas, NumPy, Matplotlib,
Seaborn, Biopython, ETE3, SciPy, statsmodels, and ArviZ. The synteny scripts
use the R packages ggplot2, reshape2, scales, and pheatmap. External analysis
programs and their versions are described in the manuscript and Supplementary
Methods.
