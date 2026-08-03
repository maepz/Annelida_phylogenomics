# Main-tree branch support

This directory contains branch-support values extracted directly from the
support-labeled trees for all recoverable internal branches in the 29 main
species-tree analyses (Supplementary Figures S7–S35).

## Files

- `main_tree_branch_support.tsv`: one row per internal branch, with all
  recovered scalar support measures, concordance/discordance vectors,
  quartet-frequency vectors, and the branch length from the qPP tree
- `main_tree_node_ids.nex`: the 29 source topologies with internal branches
  labeled as `node_<node_id>`

## Branch identifiers

`node_id` identifies the branch in the corresponding tree in
`main_tree_node_ids.nex`. It must be interpreted together with `taxon_set`,
`matrix`, `inference_method`, and `protein_evolution_model`.

`new_nodeid` is the display identifier used for the focal internal branches
shown in Supplementary Figures S7–S35. It ranges from 1 to 27 and is blank for
branches not labeled in those figures.

Branches and support values were matched by unrooted bipartition rather than
by row order or by reusing internal labels across topologically different
trees.

## Support fields

The principal scalar columns are:

- `BS`: standard bootstrap support for reconciliation-based trees
- `UFBoot`: ultrafast bootstrap support for concatenation-based trees
- `SH_aLRT`: SH-like approximate likelihood-ratio support
- `qPP`: the primary quartet-based local posterior probability (`pp1`) from
  the topology-matched qPP tree
- `sCFl`: likelihood-based site concordance factor
- `gCF`: gene concordance factor
- `qPP_tree_branch_length`: internal branch length from the qPP tree

The remaining columns retain the complete support vectors reported in the
source trees:

- gene concordance/discordance: `gCF_N`, `gDF1`, `gDF1_N`, `gDF2`,
  `gDF2_N`, `gDFP`, `gDFP_N`, `gN`
- site concordance/discordance: `sCF`, `sCF_N`, `sDF1`, `sDF1_N`, `sDF2`,
  `sDF2_N`, `sN`
- quartet support and counts: `pp1`, `pp2`, `pp3`, `q1`, `q2`, `q3`, `f1`,
  `f2`, `f3`

Blank fields indicate that the source tree did not report that measure.
Values have not been imputed. The `qPP` column is deliberately unified across
analysis types; earlier working tables used separate `pPP` and `qPP` labels
for the same primary quartet-posterior field.
