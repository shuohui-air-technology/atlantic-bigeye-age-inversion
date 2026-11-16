"""Monthly fitting and operation allocation for the Atlantic BET model."""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

import numpy as np

from .inversion import (
    build_pde_design, build_pde_state_basis, projected_gradient_nnls,
    standardize_design,
)
from .kernels import kernel_at_age
from .operations import group_cpue, project_operations
from .parameters import BET
from .selectivity import SELECTIVITY_CURVES, length_weighted_kernel
from .transfer import adjacent_misclassification_matrix

PACKAGE = Path(__file__).resolve().parent


def configuration() -> dict:
    return json.loads((PACKAGE / "model_configuration.json").read_text(encoding="utf-8"))


@lru_cache(maxsize=7)
def _observation_kernels(kernel_kind: str, epsilon: float, selectivity: str):
    """Cache the age-by-bin kernels, which are independent of the calendar."""
    if kernel_kind not in {"normal", "lognormal", "gamma"}:
        raise ValueError("kernel_kind must be normal, lognormal or gamma")
    cfg = configuration()
    transfer = cfg["transfer"]
    if epsilon not in (0.0, *transfer["epsilon_values"]):
        raise ValueError("epsilon must be 0, 0.05 or 0.10")
    if epsilon and kernel_kind != "normal":
        raise ValueError("the documented transfer configurations use the normal kernel")
    if selectivity not in ("flat", *SELECTIVITY_CURVES):
        raise ValueError("selectivity must be flat, tropical or southern")
    if selectivity != "flat" and (kernel_kind != "normal" or epsilon):
        raise ValueError("length selectivity uses the normal kernel without bin transfer")
    species = BET
    ages = np.asarray(species.age_grid_yr(), dtype=float)
    if selectivity == "flat":
        entries = [kernel_at_age(species, float(age), kernel_kind) for age in ages]
        probabilities = np.asarray([entry.probabilities for entry in entries])
        moments = np.asarray([entry.weight_moments_kg_rw for entry in entries])
    else:
        probabilities, moments = length_weighted_kernel(species, selectivity)
    if epsilon:
        boundaries = tuple(transfer["boundary_indices_zero_based"])
        matrix = np.asarray(adjacent_misclassification_matrix(4, epsilon, boundaries))
        probabilities = probabilities @ matrix.T
        moments = moments @ matrix.T
    for array in (ages, probabilities, moments):
        array.setflags(write=False)
    return ages, probabilities, moments


def _calendar(months):
    labels = tuple(months)
    if len(labels) < 2:
        raise ValueError("provide at least two observed months in YYYY-MM format")
    numbers = []
    for label in labels:
        if not isinstance(label, str) or not re.fullmatch(r"[0-9]{4}-(0[1-9]|1[0-2])", label):
            raise ValueError("months must be valid YYYY-MM labels")
        year, month = map(int, label.split("-"))
        if year == 0:
            raise ValueError("months must be valid YYYY-MM labels")
        numbers.append(12 * year + month - 1)
    if any(right <= left for left, right in zip(numbers, numbers[1:])):
        raise ValueError("months must be unique and in increasing order")
    offsets = tuple(number - numbers[0] for number in numbers)
    if offsets[-1] < 2:
        raise ValueError("the calendar must span at least two monthly transitions for the two entry endpoints")
    return labels, offsets


def observation_bundle(
    kernel_kind: str = "normal", epsilon: float = 0.0, selectivity: str = "flat",
    *, months=None,
) -> dict:
    """Build the BET kernels and design on the supplied observation calendar."""
    cfg = configuration()
    labels, offsets = _calendar(cfg["observation_months"] if months is None else months)
    knots = (0, offsets[-1] - 1)
    ages, probabilities, moments = _observation_kernels(kernel_kind, epsilon, selectivity)
    return {
        "months": labels,
        "ages": ages,
        "bin_probabilities": probabilities,
        "kernel": moments,
        "state_basis": build_pde_state_basis(BET, offsets, knots),
        "design": build_pde_design(BET, moments, offsets, knots),
    }


def states_from_coefficients(bundle: dict, coefficients) -> np.ndarray:
    values = np.asarray(coefficients, dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("coefficients must contain three finite nonnegative values")
    return np.einsum("mac,c->ma", bundle["state_basis"], values)


def fit_monthly_states(
    observations, *, months, bin_columns, kernel_kind="normal", epsilon=0.0,
    selectivity="flat", sample="full", max_seconds=None,
) -> dict:
    """Fit monthly RW CPUE on the supplied calendar, using only training rows."""
    cfg = configuration()
    months, bin_columns = tuple(months), tuple(bin_columns)
    if list(bin_columns) != cfg["weight_bin_columns_GG"]:
        raise ValueError("bin_columns must match the documented BET bin order")
    values = np.asarray(observations, dtype=float)
    n_months, n_bins = len(months), len(bin_columns)
    if values.shape != (n_months, n_bins):
        raise ValueError("observations must have one row per supplied month and four bins")
    if sample == "full":
        train = np.arange(n_months)
    elif sample in cfg["folds"]:
        if months != tuple(cfg["observation_months"]):
            raise ValueError("F1, F2 and F3 are defined only for the paper's observation months")
        train = np.asarray(cfg["folds"][sample]["train_indices"], dtype=int)
    else:
        raise ValueError("sample must be full, F1, F2 or F3")
    bundle = observation_bundle(kernel_kind, epsilon, selectivity, months=months)
    rows = np.concatenate([np.arange(index * n_bins, (index + 1) * n_bins) for index in train])
    design, observed, scales = standardize_design(bundle["design"][rows], values[train])
    limit = cfg["solver_time_limits_seconds"]["flat" if selectivity == "flat" else "length_selectivity"]
    if max_seconds is not None:
        limit = max_seconds
    fit = projected_gradient_nnls(design, observed, max_seconds=limit, **cfg["solver"])
    if not fit.converged or fit.stop_reason != "converged":
        raise RuntimeError("monthly fit did not meet the convergence criterion: " + fit.stop_reason)
    return {
        "fit": fit,
        "training_bin_scales": scales,
        "bundle": bundle,
        "states": states_from_coefficients(bundle, fit.coefficients),
    }


def allocate_operations(bundle: dict, states, operations) -> dict:
    """Apply equations (9)–(10) to operations in GG kg and positive hook counts.

    ``operations`` is a DataFrame or a column mapping with Year, Month, Hooks,
    wBET15D, wBET15U, wBET25 and wBET40. Blanks in the four catch fields mean zero.
    """
    cfg = configuration()
    months = bundle.get("months", cfg["observation_months"])
    projected = project_operations(bundle, states, operations, months)
    grouped = group_cpue(projected, bundle["ages"], cfg["group_cutpoints_yr"])
    hooks = np.asarray(operations["Hooks"], dtype=float)
    counts = grouped * hooks[:, None] / 1000.0
    totals = counts.sum(axis=1)
    proportions = np.full_like(counts, np.nan)
    np.divide(counts, totals[:, None], out=proportions, where=totals[:, None] > 0)
    return {
        **projected,
        "group_counts": counts,
        "group_cpue": grouped,
        "group_proportions": proportions,
    }
