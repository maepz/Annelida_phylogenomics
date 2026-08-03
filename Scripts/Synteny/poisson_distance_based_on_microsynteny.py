#!/usr/bin/env python3
"""Calculate Poisson-corrected pairwise distances from microsynteny blocks."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


def load_taxa(path: Path) -> list[str]:
    """Read one taxon identifier per non-comment line, preserving file order."""
    taxa = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                taxa.append(value.split()[0])
    if len(taxa) != len(set(taxa)):
        raise ValueError(f"Duplicate taxon identifiers in {path}")
    return taxa


def parse_synteny_file(path: Path, minimum_block_size: int) -> tuple[int, set[str], set[str]]:
    """Return retained anchors and gene sets for the two taxa in one block file."""
    shared_gene_count = 0
    taxon1_genes: set[str] = set()
    taxon2_genes: set[str] = set()
    current_taxon1_genes: list[str] = []
    current_taxon2_genes: list[str] = []

    def retain_current_block() -> None:
        nonlocal shared_gene_count
        if len(current_taxon1_genes) >= minimum_block_size:
            shared_gene_count += len(current_taxon1_genes)
            taxon1_genes.update(current_taxon1_genes)
            taxon2_genes.update(current_taxon2_genes)

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#block_id"):
                continue
            if line.startswith("##"):
                fields = line.split("\t")
                if len(fields) >= 5:
                    current_taxon1_genes.append(fields[1])
                    current_taxon2_genes.append(fields[3])
            else:
                retain_current_block()
                current_taxon1_genes = []
                current_taxon2_genes = []

    retain_current_block()
    return shared_gene_count, taxon1_genes, taxon2_genes


def calculate_distances(
    input_dir: Path,
    taxa: list[str],
    minimum_block_size: int,
    zero_similarity_distance: float,
) -> pd.DataFrame:
    """Build the Poisson-corrected distance matrix for the requested taxa."""
    syntenic_genes = {taxon: set() for taxon in taxa}
    shared_counts = pd.DataFrame(0.0, index=taxa, columns=taxa)

    for index, taxon1 in enumerate(taxa):
        for taxon2 in taxa[index + 1 :]:
            forward = input_dir / f"{taxon1}.{taxon2}.microsynteny"
            reverse = input_dir / f"{taxon2}.{taxon1}.microsynteny"

            if forward.exists():
                source = forward
                reversed_order = False
            elif reverse.exists():
                source = reverse
                reversed_order = True
            else:
                raise FileNotFoundError(
                    f"No microsynteny file for {taxon1} and {taxon2}"
                )

            shared, first_genes, second_genes = parse_synteny_file(
                source, minimum_block_size
            )
            if reversed_order:
                first_genes, second_genes = second_genes, first_genes

            shared_counts.at[taxon1, taxon2] = shared
            shared_counts.at[taxon2, taxon1] = shared
            syntenic_genes[taxon1].update(first_genes)
            syntenic_genes[taxon2].update(second_genes)

    base_sizes = {taxon: len(genes) for taxon, genes in syntenic_genes.items()}
    distances = pd.DataFrame(0.0, index=taxa, columns=taxa)

    for taxon1 in taxa:
        for taxon2 in taxa:
            if taxon1 == taxon2:
                continue

            denominator = base_sizes[taxon1] + base_sizes[taxon2]
            if denominator == 0:
                distances.at[taxon1, taxon2] = zero_similarity_distance
                continue

            normalized_similarity = (
                2.0 * shared_counts.at[taxon1, taxon2] / denominator
            )
            distances.at[taxon1, taxon2] = (
                -math.log(normalized_similarity)
                if normalized_similarity > 0
                else zero_similarity_distance
            )

    return distances


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate D_poisson = -ln(2*N_shared/(N_A+N_B)) from pairwise "
            "microsynteny block files."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("microsynteny_files"),
        help="Directory containing <taxon1>.<taxon2>.microsynteny files",
    )
    parser.add_argument(
        "--taxa",
        type=Path,
        default=Path("microsynteny_taxa.txt"),
        help="Text file containing one included taxon identifier per line",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("poisson_microsynteny_distance_matrix.tsv"),
        help="Output tab-delimited distance matrix",
    )
    parser.add_argument(
        "--minimum-block-size",
        type=int,
        default=5,
        help="Minimum number of anchors in a retained microsynteny block [5]",
    )
    parser.add_argument(
        "--zero-similarity-distance",
        type=float,
        default=5.0,
        help="Distance assigned when normalized similarity is zero [5.0]",
    )
    args = parser.parse_args()

    taxa = load_taxa(args.taxa)
    distances = calculate_distances(
        args.input_dir,
        taxa,
        args.minimum_block_size,
        args.zero_similarity_distance,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    distances.to_csv(args.output, sep="\t", index_label="taxon")


if __name__ == "__main__":
    main()
