# Clarifying the annelid tree of life amid deep phylogenomic conflict

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23755892.svg)](https://doi.org/10.5281/zenodo.23755892)

This repository accompanies the manuscript **“Clarifying the annelid tree of
life amid deep phylogenomic conflict.”** It contains the principal sequence
datasets, metadata, matrix definitions, gene trees, species-tree results, and
analysis-command records used in the study.

## Scope

This is a data and methods archive, not a turnkey reproducibility workflow.
The analyses were run on a computing server whose directory structure,
environment names, and filenames differ from the organization used here.
Commands in `Scripts/Analysis_commands.Rmd` document the analyses but are not
expected to run unchanged after cloning the repository. Repository filenames
have been standardized to make the archived material easier to understand.

## Repository structure

- `Data/Assemblies_and_metadata/` contains the final predicted proteomes for
  the 144-species dataset, sequence and taxon maps, assembly summaries,
  synteny-dataset taxon metadata, chromosome-level assembly metadata, and
  fossil calibrations.
- `Data/Matrices/` contains matrix membership information, gene trees,
  gene-to-species maps where needed, taxon-occupancy summaries, and gene-level
  matrix statistics, grouped into `A144/`, `A31/`, `Synteny_informed/`, and
  `Summary_statistics/`.
- `Data/Species_trees/` contains collections of inferred species trees,
  harmonized branch-support summaries, and results from constrained-topology
  analyses.
- `Data/Analyses/` contains the principal downstream result tables for
  incomplete-lineage-sorting and introgression tests, site-wise likelihood
  analyses, genome-architecture analyses, and divergence-time estimation. The
  files are kept in one directory and indexed by `Data/Analyses/README.md`.
- `Scripts/Analysis_commands.Rmd` records representative commands used for
  data processing and phylogenomic analyses. A rendered copy is provided as
  `Scripts/Analysis_commands.pdf`.
- `Scripts/` also contains cleaned standalone programs for the
  introgression, topology-testing, divergence-time, and synteny analyses,
  organized by analysis family and described in `Scripts/README.md`.

## Final phylogenomic matrices

| Matrix | Taxon set | Orthology unit | Loci |
|---|---:|---|---:|
| `OG75` | 144 species | orthogroups | 4,404 |
| `SOGall` | 144 species | single-copy orthologs | 147,655 |
| `SOG10` | 144 species | single-copy orthologs | 9,902 |
| `SOG10_tax_occ` | 144 species | occupancy-filtered single-copy orthologs | 2,162 |
| `A31_OG` | 31 annelids | orthogroups | 8,000 |
| `A31_all` | 31 annelids | single-copy orthologs | 18,305 |
| `A31_fullcov` | 31 annelids | high-occupancy single-copy orthologs | 1,580 |
| `L38_synt` | 38 species | synteny-supported orthologs | 2,296 |
| `L42_synt` | 42 species | synteny-supported orthologs | 2,299 |
| `NP_synt` | 44 species | synteny-supported orthologs | 2,069 |

The full statistics underlying this summary are in
`Data/Matrices/Summary_statistics/matrices_stats.tsv`; per-gene values are
provided in the compressed long table
`Data/Matrices/Summary_statistics/matrices_stats.longtable.tsv.gz`.
Per-species occupancy for all three synteny-informed matrices is combined in
`Data/Matrices/Synteny_informed/synt.taxon_occupancy_per_species.tsv`.

Files are grouped by dataset family rather than placed in a separate
subdirectory for every matrix:

- `Data/Matrices/A144/` contains `OG75`, `SOGall`, `SOG10`, and
  `SOG10_tax_occ` files.
- `Data/Matrices/A31/` contains `A31_OG`, `A31_all`, and `A31_fullcov` files.
- `Data/Matrices/Synteny_informed/` contains `L38_synt`, `L42_synt`, and
  `NP_synt` files.
- `Data/Matrices/Summary_statistics/` contains matrix-level and gene-level
  statistics shared across dataset families.

Within each subdirectory, filenames follow the general form
`<matrix>.<content>.<format>`, for example
`A144/SOG10.site_homogeneous.gene_trees.nwk` and
`A31/A31_OG.locus_species_map.tsv`.

### Concatenated alignment availability

Compressed concatenated alignments and their partition definitions are
included for `SOG10_tax_occ`, `A31_fullcov`, `L38_synt`, `L42_synt`, and
`NP_synt`. The substantially larger concatenated alignments for `OG75`,
`SOGall`, `SOG10`, `A31_OG`, and `A31_all` are intentionally not included.
Their associated gene trees and matrix-membership information remain
available where applicable.

## Identifiers and corrected taxon names

Internal analysis identifiers are preserved so that sequence, gene-tree, and
species-tree labels remain traceable to the original computations. Display
names and current taxonomic assignments are given in
`Data/Assemblies_and_metadata/meta_species_map.txt`.

In particular, `Nereis_procera` is retained as the internal identifier for
assembly `JAUOBU01`, while its corrected display name is *Glycera* cf.
*americana* (Glyceridae). The Supplementary Information describes the evidence
for this corrected identification.

## Constrained-topology identifiers

Constraint trees and result tables use the final manuscript labels: Topo
A1–A4, M1–M5, C1–C4, and S1–S5. Exploratory alternatives not retained in the
manuscript are not included.

The complete `freqQuad`-derived quadripartitions and qPP values and the
coverage of the retained constrained analyses are documented in
`Data/Species_trees/Constrained/README.md`.

## Main-tree branch support

Branch-support values for all recoverable internal nodes in the 29 A144, A31,
L38_synt, L42_synt, and NP_synt main-tree analyses are compiled in
`Data/Species_trees/Branch_support/main_tree_branch_support.tsv`. The
corresponding collection of all 29 trees, with matching internal-node
identifiers, is provided as
`Data/Species_trees/Branch_support/main_tree_node_ids.nex`. Bootstrap,
ultrafast-bootstrap, SH-aLRT, posterior-probability, gene-concordance, and
likelihood-based site-concordance fields are labeled explicitly. Source
coverage, node matching, and metric availability are documented in
`Data/Species_trees/Branch_support/README.md`.

## File formats

- `.pep.gz`: gzip-compressed amino-acid FASTA files
- `.fa.gz`: gzip-compressed amino-acid alignments
- `.nwk`, `.tree`, `.tre`: Newick trees
- `.nex`: NEXUS tree collections or partitions
- `.tsv`, `.txt`: tab-delimited metadata and result tables
- `.gz`: gzip-compressed text files

## Citation

Please cite this archive when using these data, and also cite the associated article when available.

Perez, M., Li, Y., Hui, F., Li, Y.-X., Li, Y., Martín-Durán, J. M., Rouse, G. W., Qian, P.-Y., & Qiu, J.-W. (2026). *Clarifying the annelid tree of life amid deep phylogenomic conflict: data and analysis archive* (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.23755892

This DOI identifies version 1.0.0. The [DOI for all versions](https://doi.org/10.5281/zenodo.23755891) always resolves to the latest archived release.

Original code is licensed under MIT; original data under CC BY 4.0. Third-party material retains its original terms. See [LICENSE.md](LICENSE.md).
