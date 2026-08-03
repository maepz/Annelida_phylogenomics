# Constrained-topology analyses

This directory contains the constraint trees and result tables used to compare
alternative phylogenetic hypotheses.

## Directory layout

- `Constraint_trees/`: NEXUS collections of the alternative constrained
  topologies for the A144, A31, and L42 taxon sets
- `Quartet_support/`: compiled ASTRAL quadripartitions and qPP values
- `Topology_tests/`: consolidated likelihood scores and SH/AU topology-test
  results
- `Best_tree/`: the reference constrained tree and its quartet-support outputs

## Compiled qPP values

`Quartet_support/constrained_quadripartitions_qPP.tsv.gz` contains the complete
quartet output compiled directly from 161 raw ASTRAL `freqQuad` files. It
preserves all three resolutions (`t1`, `t2`, and `t3`) for every internal
branch, rather than only the posterior probability of the focal `t1`
resolution.

Constrained trees differ at only a few focal nodes, so most quadripartitions
are shared among several topologies. Equivalent quadripartitions were
canonicalized by alphabetizing taxa and treating the sides within each pair
and the two pairs as unordered. Identical quadripartition/qPP records within a
matrix and gene-tree model were then collapsed to one row. Their original
topology and resolution memberships are retained in `topology_resolutions`;
for example, `A1:t1;A2:t2` means that the same record occurred as `t1` under
Topo A1 and `t2` under Topo A2.

The compressed long-format table contains the following columns:

- `taxon_set`: A144, A31, L38, L42, or NP
- `matrix`: gene matrix used to infer the gene trees
- `gene_tree_model`: site-homogeneous or site-heterogeneous protein-evolution
  model class
- `quartet_method`: ASTRAL-IV or ASTRAL-Pro3
- `topology_resolutions`: semicolon-delimited manuscript topology and
  resolution memberships
- `quadripartition`: complete four-way taxon partition for the resolution
- `qPP`: quartet-based local posterior probability for the resolution

Numeric A144 and A31 tip identifiers and abbreviated L38/L42 tip identifiers
were expanded to the species-level analysis names used in the metadata files.
For every source branch, the three qPP values sum to one within the precision
of the raw output.

The original 36,120 topology-resolution records collapse to 5,040 unique
quadripartition/qPP rows, while all 36,120 memberships remain represented in
`topology_resolutions`.

### qPP coverage

- A144 SOG10: Topo A1–A4, C1–C4, M1–M5, and S1–S5 under both models
- A144 SOG10_tax_occ: the same topologies under the site-heterogeneous model
- A144 OG75: Topo M1–M2 under the site-homogeneous model
- A144 SOGall: Topo C1, M1–M2, and S2 under the site-homogeneous model
- A31_fullcov: Topo C1–C4 and S1–S5 under both models
- A31_OG: Topo C2 under the site-homogeneous model
- A31_all: Topo C2 under both models
- L38_synt and L42_synt: Topo A1–A4, C1–C4, M1–M5, and S1–S5 under both
  models
- NP_synt: Topo A1–A4 under both models

## Constraint trees

The three NEXUS files in `Constraint_trees/` provide the constraint collections
for their respective taxon sets. Topology labels follow the final manuscript
numbering: A1–A4, M1–M5, C1–C4, and S1–S5.

## Topology tests

`Topology_tests/constrained_topology_tests.tsv` combines the likelihood and
information-criterion scores with the RELL, SH, and approximately unbiased
topology-test results. Rows are identified by `SUPERMATRIX` and `topo`.
`SUPERMATRIX` combines the matrix name with the protein-evolution model:
`.MAT` denotes the site-homogeneous model and `.C60` the site-heterogeneous
model.

The `n_free_params` value of each constrained tree has been corrected to match
the corresponding unconstrained species tree. AICc and BIC were recalculated
from `LL` using the corrected parameter count. The table is limited to the
statistics presented in the Supplementary Information topology tables: `LL`,
`LL_unconstrained`, `n_free_params`, AICc, BIC, deltaL, bp-RELL, p-SH, and
p-AU. Unused tests, sign flags, and the duplicate topology-test log-likelihood
column were removed.

The table retains the species tree and all constrained topologies used in the
manuscript. Exploratory alternatives not retained in the final topology sets
are omitted.

The separate deltaL summary table is not included because it is not read by
the final supplementary-materials source. The site-wise likelihood-deficit
figures are based on their underlying per-site likelihood data, not on this
summary table.
