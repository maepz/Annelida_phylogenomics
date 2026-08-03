# Downstream analysis results

This directory contains the principal downstream result tables and trees used
in the manuscript and Supplementary Information. Related outputs are kept
together rather than divided into method-specific subdirectories. Filenames
retain the analysis or program name needed to identify their provenance.

## Incomplete lineage sorting and introgression

| File | Description |
|---|---|
| `phytop_combined_results.tsv` | PhyTop estimates of the fractions of gene-tree discordance attributed to incomplete lineage sorting and putative introgression |
| `quibl_combined_results.csv` | Combined QuIBL branch-length-distribution results across the reported datasets and triplets |
| `abba_baba_f4_crossclade_incremental_with_q.tsv` | Patterson's D and f4 results, including multiple-testing-adjusted q-values |

## Site-wise likelihood analysis

| File | Description |
|---|---|
| `deltaL_cladetopo_gene_enrichment.tsv` | Gene-level enrichment results for sites supporting the alternative early-branching topologies |

Eight cleaned site-wise likelihood tables used for Supplementary Figure 44 are
also included here. They cover `SOG10_tax_occ`, `L38_synt`, `L42_synt`, and
`NP_synt` under both the site-homogeneous (`EDGEPROP`) and site-heterogeneous
(`SITEHET2`) protein-evolution models. Filenames follow the pattern
`<model>.<matrix>.constrained_to_topo_d.sitelh.tsv.gz`; every table is limited
to the four reported early-branching hypotheses, Topo A1–A4.

## Quartet-likelihood mapping

| File | Description |
|---|---|
| `quartet_likelihood_mapping_results.tar.gz` | Compressed archive containing 83 IQ-TREE reports, 36 quartet-definition NEXUS files, an internal README, and a SHA-256 manifest |

The archive covers the Oweniidae, Sipuncula, Myzostomida, Clitellata,
Cirratuliformia, and Siboglinidae quartet-likelihood-mapping analyses reported
in the Supplementary Information.

## Genome architecture

| File | Description |
|---|---|
| `ordered_intersection_shared_fusions_matrix.Annelida.csv` | Pairwise counts of shared chromosomal fusions relative to reconstructed ancestral annelid linkage groups |
| `microsynteny_pairwise_blocks.tar.gz` | Pairwise microsynteny block files for the 28 genomes retained in the Poisson-distance analysis; both taxon orientations are included |
| `microsynteny_pairwise_statistics.tsv` | Block-count, block-length, anchor, score, and total-block-length summaries across the broader 94-taxon microsynteny screen |
| `microsynteny_taxa.txt` | Ordered list of the 28 taxa retained in the Poisson-distance matrix |
| `poisson_microsynteny_distance_matrix.tsv` | Symmetric 28 × 28 Poisson-corrected microsynteny-distance matrix used for UPGMA clustering and the supplementary heatmap |

The distance matrix was exported from the `microsynteny` worksheet of the
analysis workbook as a plain tab-delimited table. Its row and column labels
match exactly, its diagonal is zero, and the archived pairwise block files
reproduce its values to floating-point precision.

## Divergence-time analysis

| File | Description |
|---|---|
| `topo29.LG.SOG10_tax_occ.MCMCTree.v14.tree` | MCMCTree chronogram based on the site-homogeneous protein-evolution model |
| `topo29.SITEHET2.SOG10_tax_occ.MCMCTree.v14.tree` | MCMCTree chronogram based on the site-heterogeneous protein-evolution model |
| `topo29.MCMCTree.v14.aggregated_results.tsv` | Posterior mean node ages and 95% highest-posterior-density intervals pooled across chains |
| `topo29.MCMCTree.v14.convergence.tsv` | Per-node R-hat, bulk effective-sample-size, and tail effective-sample-size diagnostics |
| `topo29.MCMCTree.v14.dryrun.tsv` | Node-age distributions from the prior-only MCMCTree run |
| `topo29.MCMCTree.v14_results.distribution.tsv` | Posterior node-age samples used to summarize the two protein-evolution models |

These are analysis outputs, not a turnkey reconstruction of the original
server-side workflow. Representative commands and cleaned custom scripts are
provided under `Scripts/`.
