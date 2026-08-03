#!/usr/bin/env bash
set -euo pipefail

# Run IQ-TREE SH and AU tests for one retained topology group.
#
# Usage:
#   constrained_topology_tests.sh TEST SUPERMATRIX MODEL TREES [SCORES_FILE]
#
# Examples:
#   constrained_topology_tests.sh A SOG10_all MAT 1234
#   constrained_topology_tests.sh S L42_synt C60 12345
#
# TEST is one of A, M, C, or S. MODEL is MAT (edge-proportional
# partitioned model) or C60 (LG+C60+F+G). TREES is the continuous string
# of hypothesis numbers to test within that group.
#
# The script refers to the original analysis workspace, not the GitHub
# repository layout. Set ANALYSIS_ROOT before adapting it.

if (( $# < 4 || $# > 5 )); then
    echo \
        "Usage: $0 TEST SUPERMATRIX MODEL TREES [SCORES_FILE]" \
        >&2
    exit 2
fi

TEST=$1
SUPERMATRIX=$2
MODEL=$3
TREES=$4

ANALYSIS_ROOT=${ANALYSIS_ROOT:-/path/to/original_analysis_workspace}
DATA_DIR=${DATA_DIR:-"${ANALYSIS_ROOT}/data"}
SCRATCH_DIR=${SCRATCH_DIR:-"${ANALYSIS_ROOT}/scratch"}
TREE_BUILDING_DIR="${SCRATCH_DIR}/tree_building"
SCORES_FILE=${5:-"${TREE_BUILDING_DIR}/constrained_to_topo_d.tree_scores.tsv"}
CPUS=${CPUS:-12}

case "$TEST" in
    A|M|C|S) ;;
    *)
        echo "TEST must be A, M, C, or S (got: $TEST)" >&2
        exit 2
        ;;
esac

case "$MODEL" in
    MAT)
        base="EDGEPROP.${SUPERMATRIX}"
        ;;
    C60)
        base="SITEHET2.${SUPERMATRIX}"
        ;;
    *)
        echo "MODEL must be MAT or C60 (got: $MODEL)" >&2
        exit 2
        ;;
esac

if [[ ! -s "$SCORES_FILE" ]]; then
    echo "Missing constrained-tree score table: $SCORES_FILE" >&2
    exit 1
fi

work_dir="${TREE_BUILDING_DIR}/${SUPERMATRIX}"
if [[ ! -d "$work_dir" ]]; then
    echo "Missing analysis directory: $work_dir" >&2
    exit 1
fi
cd "$work_dir"

pattern="^topo${TEST}[${TREES}]d$"
score_dataset="${SUPERMATRIX}.${MODEL}"

reference_topology=$(
    awk -F'\t' -v dataset="$score_dataset" -v pattern="$pattern" '
        { sub(/\r$/, "", $0) }
        NR == 1 { next }
        $1 == dataset && $2 ~ pattern {
            value = $3 + 0
            if (!found || value > best) {
                best = value
                topology = $2
                found = 1
            }
        }
        END {
            if (found) {
                print topology
            }
        }
    ' "$SCORES_FILE"
)

if [[ -z "$reference_topology" ]]; then
    echo \
        "No matching topology in $SCORES_FILE for $score_dataset and $pattern" \
        >&2
    exit 1
fi
echo "Selected reference topology: $reference_topology"

prefix="${base}.constrained_to_topo_${TEST}d"
tree_list="${prefix}.treels"
reference_tree="${base}.constrained_to_${reference_topology}.treefile"
candidate_glob="${base}.constrained_to_topo${TEST}[${TREES}]d.treefile"

shopt -s nullglob
candidate_trees=( $candidate_glob )
shopt -u nullglob

if (( ${#candidate_trees[@]} == 0 )); then
    echo "No tree files found for pattern: $candidate_glob" >&2
    exit 1
fi
cat "${candidate_trees[@]}" > "$tree_list"

if [[ "$MODEL" == "MAT" ]]; then
    iqtree2 \
        -S "${SUPERMATRIX}_data" \
        -p "${SUPERMATRIX}.partition.nex" \
        --seqtype AA \
        --prefix "$prefix" \
        -te "$reference_tree" \
        -ntmax "$CPUS" \
        -nt "$CPUS" \
        -z "$tree_list" \
        -n 0 \
        -zb 10000 \
        -zw \
        -au \
        -redo
else
    reference_alignment=${REF_ALN:-"${work_dir}/${SUPERMATRIX}.concatenated.fa"}
    reference_site_frequencies=${REF_SITEFREQ:-"${work_dir}/SITEHET2.${SUPERMATRIX}.sitefreq"}
    iqtree2 \
        -s "$reference_alignment" \
        -m "LG+C60+F+G" \
        -fs "$reference_site_frequencies" \
        --seqtype AA \
        --prefix "$prefix" \
        -te "$reference_tree" \
        -ntmax "$CPUS" \
        -nt "$CPUS" \
        -z "$tree_list" \
        -n 0 \
        -zb 10000 \
        -zw \
        -au \
        --undo
fi

report="${prefix}.iqtree"
if [[ ! -s "$report" ]]; then
    report="${prefix}.log"
fi
if [[ ! -s "$report" ]]; then
    echo "Could not find IQ-TREE report for prefix: $prefix" >&2
    exit 1
fi

labels=""
for (( index=0; index<${#TREES}; index++ )); do
    topology_number=${TREES:index:1}
    if (( index > 0 )); then
        labels+=","
    fi
    labels+="topo${TEST}${topology_number}d"
done

output="${prefix}.AUtest.tsv"
awk -v labels="$labels" -v dataset="$score_dataset" '
    BEGIN {
        OFS = "\t"
        in_table = 0
        printed_header = 0
        number_labels = split(labels, topology_labels, /,/)
    }
    function print_header() {
        if (!printed_header) {
            print "SUPERMATRIX", "topo", "logL", "deltaL",
                  "bp_RELL", "bp_RELL_sign",
                  "p_KH", "p_KH_sign",
                  "p_SH", "p_SH_sign",
                  "p_WKH", "p_WKH_sign",
                  "p_WSH", "p_WSH_sign",
                  "c_ELW", "c_ELW_sign",
                  "p_AU", "p_AU_sign"
            printed_header = 1
        }
    }
    $1 == "Tree" && $2 == "logL" && $3 == "deltaL" {
        in_table = 1
        print_header()
        next
    }
    in_table && /^[[:space:]]*-{5,}[[:space:]]*$/ { next }
    in_table && /^[[:space:]]*$/ { in_table = 0; next }
    in_table && /^[[:space:]]*deltaL[[:space:]]*:/ {
        in_table = 0
        next
    }
    in_table && /^[[:space:]]*ALISIM COMMAND/ {
        in_table = 0
        next
    }
    in_table && /^[[:space:]]*[0-9]+[[:space:]]/ {
        tree_index = $1 + 0
        topology = (
            tree_index >= 1 && tree_index <= number_labels
            ? topology_labels[tree_index]
            : "Tree" tree_index
        )
        log_likelihood = $2
        delta_likelihood = $3
        column = 4

        for (metric=1; metric<=7; metric++) {
            value[metric] = (column <= NF ? $(column) : "")
            sign = (column + 1 <= NF ? $(column + 1) : "")
            if (sign == "+" || sign == "-") {
                significance[metric] = sign
                column += 2
            } else {
                significance[metric] = ""
                column += 1
            }
        }

        print dataset, topology, log_likelihood, delta_likelihood,
              value[1], significance[1],
              value[2], significance[2],
              value[3], significance[3],
              value[4], significance[4],
              value[5], significance[5],
              value[6], significance[6],
              value[7], significance[7]
    }
' "$report" > "$output"

echo "Wrote AU-test table: ${work_dir}/${output}"

# Set COMBINE_AU_RESULTS=1 to rebuild the cross-dataset result table after
# the individual test finishes.
if [[ ${COMBINE_AU_RESULTS:-0} == 1 ]]; then
    shopt -s nullglob
    result_tables=( "${TREE_BUILDING_DIR}"/*/[ES]*AUtest.tsv )
    shopt -u nullglob

    if (( ${#result_tables[@]} == 0 )); then
        echo "No AU-test tables found under $TREE_BUILDING_DIR" >&2
        exit 1
    fi

    combined="${TREE_BUILDING_DIR}/constrained_to_topo_d.AUtest.tsv"
    sort -u -r "${result_tables[@]}" > "$combined"
    results_dir="${DATA_DIR}/trees/key_results_topo_d"
    mkdir -p "$results_dir"
    cp "$combined" "$results_dir/"
    echo "Wrote combined table: $combined"
fi
