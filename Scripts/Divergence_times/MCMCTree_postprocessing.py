#!/usr/bin/env python3
"""Summarize MCMCTree priors, posterior ages, and chain diagnostics.

This script records the post-processing used for the two protein-evolution
models in the manuscript. It refers to the original analysis workspace, not
the organization of the GitHub repository.
"""

import argparse
import os
import re
from pathlib import Path

DEFAULT_MODELS = [
    "LG.SOG10_tax_occ",
    "SITEHET2.SOG10_tax_occ",
]


def read_chain(path):
    """Read one MCMCTree sample file and retain node-age columns."""
    table = pd.read_csv(path, sep=r"\s+", header=0, comment="#")
    node_columns = [column for column in table.columns if column.startswith("t_")]
    if not node_columns:
        raise ValueError(f"No node-age columns beginning with 't_' in {path}")
    return table[node_columns]


def common_chain_arrays(paths, extra_burn_fraction=0.0):
    """Return equally sized arrays with the node columns shared by all chains."""
    tables = [read_chain(path) for path in paths]
    common_columns = [
        column
        for column in tables[0].columns
        if all(column in table.columns for table in tables)
    ]
    if not common_columns:
        raise ValueError("The MCMCTree chains have no common node-age columns")

    arrays = []
    for table in tables:
        start = int(len(table) * extra_burn_fraction)
        arrays.append(table.loc[:, common_columns].iloc[start:].to_numpy(float))

    minimum_draws = min(array.shape[0] for array in arrays)
    arrays = [array[:minimum_draws] for array in arrays]
    return common_columns, arrays


def chain_diagnostics(paths, extra_burn_fraction=0.0):
    """Calculate rank-normalized R-hat and bulk and tail ESS."""
    try:
        import arviz as az
    except ImportError as error:
        raise ImportError(
            "MCMCTree diagnostics require the Python package 'arviz'"
        ) from error

    node_columns, arrays = common_chain_arrays(
        paths,
        extra_burn_fraction=extra_burn_fraction,
    )
    posterior = np.stack(arrays, axis=0)
    inference_data = az.from_dict(
        posterior={"theta": posterior},
        coords={"node": node_columns},
        dims={"theta": ["node"]},
    )

    metrics = {
        "rhat": az.rhat(inference_data)["theta"].to_pandas(),
        "ess_bulk": az.ess(
            inference_data,
            method="bulk",
        )["theta"].to_pandas(),
        "ess_tail": az.ess(
            inference_data,
            method="tail",
        )["theta"].to_pandas(),
    }
    return (
        pd.DataFrame(metrics)
        .rename_axis("node")
        .reset_index()
        .melt(id_vars="node", var_name="variable", value_name="value")
    )


def highest_posterior_density(samples, probability=0.95):
    """Return the shortest interval containing the requested sample mass."""
    values = np.asarray(samples, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("No finite posterior samples")
    if not 0 < probability < 1:
        raise ValueError("probability must be between zero and one")

    ordered = np.sort(values)
    interval_size = int(np.ceil(probability * ordered.size))
    widths = (
        ordered[interval_size - 1:]
        - ordered[: ordered.size - interval_size + 1]
    )
    start = int(np.argmin(widths))
    return (
        float(ordered[start]),
        float(ordered[start + interval_size - 1]),
        float(values.mean()),
    )


def posterior_summary(paths, calibration, model):
    """Pool four chains and calculate posterior means and 95% HPD intervals."""
    tables = [read_chain(path) for path in paths]
    common_columns = [
        column
        for column in tables[0].columns
        if all(column in table.columns for table in tables)
    ]
    pooled = pd.concat(
        [table[common_columns] for table in tables],
        ignore_index=True,
    )

    rows = []
    for node in common_columns:
        low, high, mean = highest_posterior_density(pooled[node])
        rows.append(
            {
                "node": node,
                "HPD_low": low,
                "HPD_high": high,
                "mean": mean,
                "calibration": calibration,
                "alignment": model,
            }
        )
    return pd.DataFrame(rows)


def extract_prior_section(mcmctree_output):
    """Extract the posterior-mean table printed during an MCMCTree dry run."""
    start_pattern = re.compile(r"^Posterior mean")
    end_pattern = re.compile(r"^mu")
    selected = []
    in_section = False

    with open(mcmctree_output, encoding="utf-8") as stream:
        for line in stream:
            if start_pattern.match(line):
                in_section = True
            if in_section:
                selected.append(line)
            if in_section and end_pattern.match(line):
                break

    if selected:
        selected = selected[:-1]
    if len(selected) > 2:
        selected = selected[2:]
    return selected


def parse_prior_summary(mcmctree_output, calibration):
    """Parse the dry-run node-age summary into a tabular form."""
    translation = str.maketrans({"(": " ", ")": " ", ",": " "})
    rows = []

    for line in extract_prior_section(mcmctree_output):
        parts = line.strip().translate(translation).split()
        if not parts:
            continue
        if len(parts) < 9:
            raise ValueError(f"Unexpected dry-run line: {line.rstrip()}")
        rows.append(
            {
                "node": parts[0],
                "mean": float(parts[1]),
                "equal_tail_low": float(parts[2]),
                "equal_tail_high": float(parts[3]),
                "HPD_low": float(parts[4]),
                "HPD_high": float(parts[5]),
                "HPD_CI": float(parts[6]),
                "jnode": f"{parts[7]} {parts[8]}",
                "calibration": calibration,
            }
        )

    table = pd.DataFrame(rows)
    if table.empty:
        raise ValueError(f"No dry-run node summaries found in {mcmctree_output}")
    table["equal_tail_CI"] = (
        table["equal_tail_high"] - table["equal_tail_low"]
    )
    return table


def chain_paths(time_tree_dir, topology, calibration, model, runs):
    """Construct the original server-side paths for a set of chains."""
    return [
        (
            time_tree_dir
            / model
            / calibration
            / f"run_{run}"
            / (
                f"run_{run}.{topology}_{calibration}."
                f"{model}.mcmctree.log"
            )
        )
        for run in runs
    ]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize the MCMCTree v14 dry run and four posterior chains "
            "for the two protein-evolution models."
        )
    )
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=Path(
            os.environ.get(
                "ANALYSIS_ROOT",
                "/path/to/original_analysis_workspace",
            )
        ),
    )
    parser.add_argument("--topology", default="topo29")
    parser.add_argument("--calibration", default="v14")
    parser.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    parser.add_argument("--runs", nargs="+", type=int, default=[1, 2, 3, 4])
    parser.add_argument(
        "--extra-burn-fraction",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--dryrun-output",
        type=Path,
        default=None,
        help="MCMCTree dry-run output; inferred from the analysis root if omitted",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; defaults to ANALYSIS_ROOT/data/trees",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    global np, pd
    try:
        import numpy as np
        import pandas as pd
    except ImportError as error:
        raise ImportError(
            "MCMCTree post-processing requires NumPy and pandas"
        ) from error

    time_tree_dir = (
        args.analysis_root / "scratch" / "tree_building" / "time_tree"
    )
    output_dir = args.output_dir or args.analysis_root / "data" / "trees"
    output_dir.mkdir(parents=True, exist_ok=True)

    diagnostic_tables = []
    posterior_tables = []
    for model in args.models:
        paths = chain_paths(
            time_tree_dir,
            args.topology,
            args.calibration,
            model,
            args.runs,
        )
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                "Missing MCMCTree chain files:\n" + "\n".join(missing)
            )

        diagnostics = chain_diagnostics(
            paths,
            extra_burn_fraction=args.extra_burn_fraction,
        )
        diagnostics["alignment"] = model
        diagnostics["calib"] = args.calibration
        diagnostic_tables.append(diagnostics)
        posterior_tables.append(
            posterior_summary(
                paths,
                calibration=args.calibration,
                model=model,
            )
        )

    prefix = f"{args.topology}.MCMCTree.{args.calibration}"
    diagnostics = pd.concat(diagnostic_tables, ignore_index=True)
    diagnostics.to_csv(
        output_dir / f"{prefix}.convergence.tsv",
        sep="\t",
        index=False,
    )

    posterior = pd.concat(posterior_tables, ignore_index=True)
    posterior.to_csv(
        output_dir / f"{prefix}.aggregated_results.tsv",
        sep="\t",
        index=True,
    )

    dryrun_output = args.dryrun_output
    if dryrun_output is None:
        dryrun_output = (
            time_tree_dir
            / f"{args.calibration}_dryrun"
            / (
                f"{args.topology}.with_time_calibrations_"
                f"{args.calibration}.SOG10_tax_occ_10occtrim.mcmctree.out"
            )
        )
    if not dryrun_output.is_file():
        raise FileNotFoundError(dryrun_output)

    prior = parse_prior_summary(
        dryrun_output,
        calibration=args.calibration,
    )
    prior.to_csv(
        output_dir / f"{prefix}.dryrun.tsv",
        sep="\t",
        index=True,
    )

    print(f"Wrote MCMCTree summaries to {output_dir}")


if __name__ == "__main__":
    main()
