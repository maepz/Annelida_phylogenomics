#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from collections import defaultdict, Counter


def read_bed(bedfile):
    """
    Read BED-like file and build:
    gene -> (chrom, order, start, end)

    Expected columns (at least 4):
    chrom  start  end  gene
    """
    chrom_genes = defaultdict(list)
    with open(bedfile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                continue
            chrom, start, end, gene = fields[:4]
            chrom_genes[chrom].append((int(start), int(end), gene))

    gene_map = {}
    for chrom, genes in chrom_genes.items():
        genes.sort(key=lambda x: (x[0], x[1], x[2]))
        for order, (start, end, gene) in enumerate(genes):
            gene_map[gene] = (chrom, order, start, end)
    return gene_map


def read_mbh(mbhfile, bed1_map, bed2_map):
    """
    Read MBH file.

    Expected columns:
    gene1 gene2 [score]

    If score is absent, score=1.0.
    """
    anchors = []
    missing1 = 0
    missing2 = 0

    with open(mbhfile) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue

            g1, g2 = fields[:2]
            score = 1.0
            if len(fields) >= 3:
                try:
                    score = float(fields[2])
                except ValueError:
                    score = 1.0

            if g1 not in bed1_map:
                missing1 += 1
                continue
            if g2 not in bed2_map:
                missing2 += 1
                continue

            chr1, idx1, start1, end1 = bed1_map[g1]
            chr2, idx2, start2, end2 = bed2_map[g2]
            anchors.append({
                "g1": g1,
                "chr1": chr1,
                "idx1": idx1,
                "start1": start1,
                "end1": end1,
                "g2": g2,
                "chr2": chr2,
                "idx2": idx2,
                "start2": start2,
                "end2": end2,
                "score": score,
            })

    return anchors, missing1, missing2


def finalize_block(block_id, chrpair, anchors, orientation_hint=None):
    idx1s = [a["idx1"] for a in anchors]
    idx2s = [a["idx2"] for a in anchors]
    score = sum(a["score"] for a in anchors)
    n_anchors = len(anchors)
    block_len1 = max(idx1s) - min(idx1s) + 1
    block_len2 = max(idx2s) - min(idx2s) + 1

    diffs2 = []
    for a, b in zip(anchors[:-1], anchors[1:]):
        d2 = b["idx2"] - a["idx2"]
        if d2 > 0:
            diffs2.append("+")
        elif d2 < 0:
            diffs2.append("-")

    if diffs2:
        orientation = Counter(diffs2).most_common(1)[0][0]
    elif orientation_hint in ("+", "-"):
        orientation = orientation_hint
    else:
        orientation = "+"

    inversion_noise = sum(1 for a, b in zip(anchors[:-1], anchors[1:])
                          if ((b["idx2"] - a["idx2"]) > 0 and orientation == "-")
                          or ((b["idx2"] - a["idx2"]) < 0 and orientation == "+"))

    return {
        "block_id": block_id,
        "chr1": chrpair[0],
        "chr2": chrpair[1],
        "orientation": orientation,
        "n_anchors": n_anchors,
        "score": round(score, 3),
        "sp1_start_idx": min(idx1s),
        "sp1_end_idx": max(idx1s),
        "sp2_start_idx": min(idx2s),
        "sp2_end_idx": max(idx2s),
        "block_len1": block_len1,
        "block_len2": block_len2,
        "inversion_noise": inversion_noise,
        "anchors": anchors[:],
    }


def cluster_blocks(anchors, max_gap=20, min_anchors=4, max_inv=1, max_inv_span=3):
    """
    Single-linkage style chaining inspired by jcvi/MCscan.

    Allow limited inversion noise:
    - max_inv: maximum number of inversion-noise steps allowed per block
    - max_inv_span: maximum |d2| for an inversion-noise step

    Block orientation is inferred from the dominant trend on species2 gene order.
    """
    grouped = defaultdict(list)
    for a in anchors:
        grouped[(a["chr1"], a["chr2"])].append(a)

    all_blocks = []
    block_id = 1

    for chrpair, arr in grouped.items():
        arr.sort(key=lambda x: (x["idx1"], x["idx2"]))
        if not arr:
            continue

        current = [arr[0]]
        direction = None
        inv_count = 0

        for curr in arr[1:]:
            prev = current[-1]
            d1 = curr["idx1"] - prev["idx1"]
            d2 = curr["idx2"] - prev["idx2"]

            break_block = False
            append_current = False
            next_direction = direction
            next_inv_count = inv_count

            if d1 <= 0 or abs(d1) > max_gap or abs(d2) > max_gap:
                break_block = True
            else:
                step_dir = "+" if d2 >= 0 else "-"
                if direction is None:
                    next_direction = step_dir
                    append_current = True
                elif step_dir == direction:
                    append_current = True
                else:
                    if abs(d2) <= max_inv_span and (inv_count + 1) <= max_inv:
                        append_current = True
                        next_inv_count = inv_count + 1
                    else:
                        break_block = True

            if break_block:
                if len(current) >= min_anchors:
                    all_blocks.append(finalize_block(block_id, chrpair, current, orientation_hint=direction))
                    block_id += 1
                current = [curr]
                direction = None
                inv_count = 0
            elif append_current:
                current.append(curr)
                direction = next_direction
                inv_count = next_inv_count

        if len(current) >= min_anchors:
            all_blocks.append(finalize_block(block_id, chrpair, current, orientation_hint=direction))
            block_id += 1

    return all_blocks


def write_blocks(blocks, outfile):
    with open(outfile, "w") as out:
        out.write("#block_id\tchr1\tchr2\torientation\tn_anchors\tscore\tblock_len1\tblock_len2\tsp1_start_idx\tsp1_end_idx\tsp2_start_idx\tsp2_end_idx\tinversion_noise\n")
        for b in blocks:
            out.write(
                f"{b['block_id']}\t{b['chr1']}\t{b['chr2']}\t{b['orientation']}\t{b['n_anchors']}\t{b['score']}\t"
                f"{b['block_len1']}\t{b['block_len2']}\t{b['sp1_start_idx']}\t{b['sp1_end_idx']}\t"
                f"{b['sp2_start_idx']}\t{b['sp2_end_idx']}\t{b['inversion_noise']}\n"
            )
            for a in b["anchors"]:
                out.write(
                    f"## {b['block_id']}\t{a['g1']}\t{a['idx1']}\t{a['g2']}\t{a['idx2']}\t{a['score']}\n"
                )


def write_simple(blocks, outfile):
    with open(outfile, "w") as out:
        out.write("block_id\tstart_gene1\tend_gene1\tstart_gene2\tend_gene2\tscore\torientation\tn_anchors\tblock_len1\tblock_len2\tinversion_noise\n")
        for b in blocks:
            anchors = b["anchors"]
            start_gene1 = min(anchors, key=lambda x: x["idx1"])["g1"]
            end_gene1 = max(anchors, key=lambda x: x["idx1"])["g1"]
            start_gene2 = min(anchors, key=lambda x: x["idx2"])["g2"]
            end_gene2 = max(anchors, key=lambda x: x["idx2"])["g2"]
            out.write(
                f"{b['block_id']}\t{start_gene1}\t{end_gene1}\t{start_gene2}\t{end_gene2}\t"
                f"{b['score']}\t{b['orientation']}\t{b['n_anchors']}\t{b['block_len1']}\t{b['block_len2']}\t{b['inversion_noise']}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Identify pairwise microsynteny blocks from BED + MBH, allowing limited inversion noise."
    )
    parser.add_argument("--bed1", required=True, help="BED file for species 1")
    parser.add_argument("--bed2", required=True, help="BED file for species 2")
    parser.add_argument("--mbh", required=True, help="MBH file: gene1 gene2 [score]")
    parser.add_argument("-o", "--output", required=True, help="Main output block file")
    parser.add_argument("--simple", default=None, help="Optional simple summary output")
    parser.add_argument("--max-gap", type=int, default=20,
                        help="Maximum allowed gap in gene order between neighboring anchors [20]")
    parser.add_argument("--min-anchors", type=int, default=4,
                        help="Minimum anchors to define a block [4]")
    parser.add_argument("--max-inv", type=int, default=1,
                        help="Maximum allowed inversion-noise steps per block [1]")
    parser.add_argument("--max-inv-span", type=int, default=3,
                        help="Maximum |d2| for one inversion-noise step [3]")
    args = parser.parse_args()

    bed1_map = read_bed(args.bed1)
    bed2_map = read_bed(args.bed2)
    anchors, missing1, missing2 = read_mbh(args.mbh, bed1_map, bed2_map)

    if not anchors:
        sys.stderr.write("No usable anchors found. Check MBH gene IDs and BED gene IDs.\n")
        sys.exit(1)

    blocks = cluster_blocks(
        anchors,
        max_gap=args.max_gap,
        min_anchors=args.min_anchors,
        max_inv=args.max_inv,
        max_inv_span=args.max_inv_span,
    )

    write_blocks(blocks, args.output)
    if args.simple:
        write_simple(blocks, args.simple)

    sys.stderr.write(f"BED1 genes loaded: {len(bed1_map)}\n")
    sys.stderr.write(f"BED2 genes loaded: {len(bed2_map)}\n")
    sys.stderr.write(f"Usable anchors: {len(anchors)}\n")
    sys.stderr.write(f"Missing anchors (species1 IDs absent in BED1): {missing1}\n")
    sys.stderr.write(f"Missing anchors (species2 IDs absent in BED2): {missing2}\n")
    sys.stderr.write(f"Total blocks: {len(blocks)}\n")
    sys.stderr.write(f"Main output: {args.output}\n")
    if args.simple:
        sys.stderr.write(f"Simple output: {args.simple}\n")


if __name__ == "__main__":
    main()
