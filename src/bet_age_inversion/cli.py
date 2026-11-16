"""Command-line input, fitting, and output for the Atlantic BET package."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

from .workflow import allocate_operations, configuration, fit_monthly_states

MONTHLY_FIELDS = (
    "CPUE_RW_lt15GG", "CPUE_RW_15to25GG", "CPUE_RW_25to40GG", "CPUE_RW_ge40GG",
)
OUTPUT_FIELDS = (
    "cpue_0_to_lt3", "cpue_3_to_lt5", "cpue_ge5",
    "proportion_0_to_lt3", "proportion_3_to_lt5", "proportion_ge5",
)


def load_monthly_input(path):
    with path.open(newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames != ["month", *MONTHLY_FIELDS]:
            raise ValueError("monthly CSV has incorrect columns or bin order")
        rows = list(reader)
    months = [row["month"] for row in rows]
    observations = np.array([[float(row[field]) for field in MONTHLY_FIELDS] for row in rows])
    return months, observations


def load_operations(path, bin_columns):
    required = ("Year", "Month", "Hooks", *bin_columns)
    with path.open(newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError("operation CSV must include " + ", ".join(required))
        rows = list(reader)
    return {field: [row[field] for row in rows] for field in required}


def check_reference(example_dir, months, result, indicators):
    reference = json.loads((example_dir / "reference_result.json").read_text(encoding="utf-8"))
    if not result["converged"] or not reference["converged"]:
        raise RuntimeError("expected a converged flat-reference fit")
    for key in ("coefficients", "objective", "projected_gradient_inf", "monthly_group_proportion_mean"):
        np.testing.assert_allclose(result[key], reference[key], rtol=1e-7, atol=1e-7, err_msg=key)
    with (example_dir / "reference_monthly_age_indicators.csv").open(newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames != ["month", *OUTPUT_FIELDS]:
            raise ValueError("reference monthly CSV has incorrect columns")
        expected_rows = list(reader)
    if months != [row["month"] for row in expected_rows]:
        raise ValueError("reference months do not match the synthetic input")
    expected = [[float(row[field]) for field in OUTPUT_FIELDS] for row in expected_rows]
    np.testing.assert_allclose(indicators, expected, rtol=1e-7, atol=1e-7)


def main(*, example_dir: Path | None = None):
    """Run supplied CSVs, or the frozen example when invoked by its script."""
    description = (
        "Run the frozen synthetic example or fit your own CSVs."
        if example_dir is not None
        else "Estimate Atlantic bigeye tuna age-group indicators from your CSVs."
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--selectivity", choices=("flat", "tropical", "southern"), default="flat")
    parser.add_argument("--output-dir", type=Path, default=Path("example_output"))
    require_inputs = example_dir is None
    parser.add_argument("--monthly-csv", type=Path, required=require_inputs, help="monthly weight CPUE CSV")
    parser.add_argument("--operations-csv", type=Path, required=require_inputs, help="operation-level catch and effort CSV")
    args = parser.parse_args()
    if (args.monthly_csv is None) != (args.operations_csv is None):
        parser.error("supply both --monthly-csv and --operations-csv, or neither for the synthetic example")

    cfg = configuration()
    synthetic = args.monthly_csv is None
    months, observations = load_monthly_input(
        example_dir / "synthetic_monthly_weights.csv" if synthetic else args.monthly_csv
    )
    fitted = fit_monthly_states(
        observations, months=months, bin_columns=cfg["weight_bin_columns_GG"],
        selectivity=args.selectivity,
        max_seconds=cfg["solver_time_limits_seconds"]["synthetic_example"] if synthetic else None,
    )
    fit, bundle, states = fitted["fit"], fitted["bundle"], fitted["states"]
    # Eq. (9) with the supplied monthly inputs already expressed as kg RW per 1,000 hooks.
    bin_moments = states @ bundle["kernel"]
    positive_bins = observations > 0
    if np.any(positive_bins & (bin_moments <= 0)):
        raise ValueError("a positive monthly weight bin has no support in the fitted state")
    conversion = np.divide(observations, bin_moments, out=np.zeros_like(observations), where=positive_bins)
    cpue_by_age = states * (conversion @ bundle["bin_probabilities"].T)
    ages = bundle["ages"]
    lower, upper = cfg["group_cutpoints_yr"]
    groups = np.column_stack((ages < lower, (ages >= lower) & (ages < upper), ages >= upper))
    cpue = cpue_by_age @ groups
    totals = cpue.sum(axis=1, keepdims=True)
    shares = np.divide(cpue, totals, out=np.full_like(cpue, np.nan), where=totals > 0)
    indicators = np.column_stack((cpue, shares))
    result = {
        "converged": fit.converged, "iterations": fit.iterations,
        "objective": fit.objective, "projected_gradient_inf": fit.projected_gradient_inf,
        "coefficients": fit.coefficients.tolist(),
        "monthly_group_proportion_mean": shares[totals[:, 0] > 0].mean(axis=0).tolist(),
        "software": {"python": sys.version.split()[0], "numpy": np.__version__},
        "synthetic_only": synthetic, "selectivity": args.selectivity,
    }
    if synthetic and args.selectivity == "flat":
        check_reference(example_dir, months, result, indicators)

    if synthetic:
        operations = json.loads((example_dir / "synthetic_operations.json").read_text(encoding="utf-8"))
    else:
        operations = load_operations(args.operations_csv, cfg["weight_bin_columns_GG"])
    output = allocate_operations(bundle, states, operations)
    np.testing.assert_allclose(output["group_cpue"].sum(axis=1), output["cpue_by_age"].sum(axis=1))
    np.testing.assert_allclose(output["group_proportions"][output["positive"]].sum(axis=1), 1)
    np.testing.assert_array_equal(output["group_counts"][~output["positive"]], 0)
    np.testing.assert_array_equal(output["group_cpue"][~output["positive"]], 0)
    if not np.isnan(output["group_proportions"][~output["positive"]]).all():
        raise RuntimeError("expected undefined proportions for the zero-catch operation")
    allocation = {
        "synthetic_only": synthetic, "selectivity": args.selectivity,
        "state_source": "fit to the bundled synthetic monthly CPUE" if synthetic else "fit to the supplied monthly CPUE",
        "group_order": ["0 to <3 yr", "3 to <5 yr", ">=5 yr"],
        "group_counts_fish_equivalents": output["group_counts"].tolist(),
        "group_cpue_fish_equivalents_per_1000_hooks": output["group_cpue"].tolist(),
        "group_proportions": [
            values.tolist() if positive else None
            for values, positive in zip(output["group_proportions"], output["positive"], strict=True)
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, record in (("fit_result.json", result), ("operation_results.json", allocation)):
        (args.output_dir / filename).write_text(
            json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    with (args.output_dir / "monthly_age_indicators.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["month", *OUTPUT_FIELDS])
        writer.writerows(
            [month, *(None if np.isnan(value) else value for value in row)]
            for month, row in zip(months, indicators, strict=True)
        )
    source = "synthetic" if synthetic else "supplied"
    print(
        f"{args.selectivity}: fitted {len(months)} {source} months and allocated "
        f"{len(operations['Year'])} {source} operations."
    )


if __name__ == "__main__":
    main()
