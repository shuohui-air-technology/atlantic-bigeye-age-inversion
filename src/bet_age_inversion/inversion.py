"""Linear monthly PDE design and bounded nonnegative least squares."""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
import time

import numpy as np

from .parameters import DT_YR, SpeciesParameters


@dataclass(frozen=True)
class FitResult:
    """Deterministic termination record for a projected-gradient NNLS fit.

    ``projected_gradient_inf`` is NaN only when a time limit is reached before
    solver initialization; in that case it means "not computed", not convergence.
    """

    coefficients: np.ndarray
    converged: bool
    iterations: int
    objective: float
    projected_gradient_inf: float
    elapsed_seconds: float
    stop_reason: str


def _readonly(array: np.ndarray) -> np.ndarray:
    array.setflags(write=False)
    return array


def piecewise_linear_basis(
    n_transitions: int,
    knot_indices: tuple[int, ...],
) -> np.ndarray:
    """Return convex linear interpolation weights at every transition."""
    if isinstance(n_transitions, bool) or not isinstance(n_transitions, Integral) or n_transitions <= 0:
        raise ValueError("n_transitions must be a positive integer")
    try:
        knots = tuple(knot_indices)
    except TypeError as exc:
        raise ValueError("knot_indices must be a nonempty sequence") from exc
    if not knots:
        raise ValueError("knot_indices must not be empty")
    if any(isinstance(knot, bool) or not isinstance(knot, Integral) for knot in knots):
        raise ValueError("knot indices must be integers")
    knots = tuple(int(knot) for knot in knots)
    if any(right <= left for left, right in zip(knots, knots[1:])):
        raise ValueError("knot indices must be strictly increasing")
    if knots[0] != 0 or knots[-1] != int(n_transitions) - 1:
        raise ValueError("knots must span the first and last transition")
    if any(knot < 0 or knot >= n_transitions for knot in knots):
        raise ValueError("knot index is outside the transition range")

    basis = np.zeros((int(n_transitions), len(knots)), dtype=float)
    if len(knots) == 1:
        basis[0, 0] = 1.0
        return basis
    for column, (left, right) in enumerate(zip(knots, knots[1:])):
        width = right - left
        fractions = np.arange(width + 1, dtype=float) / float(width)
        basis[left : right + 1, column] = 1.0 - fractions
        basis[left : right + 1, column + 1] = fractions
    return basis


def _monthly_survival(species: SpeciesParameters) -> np.ndarray:
    ages = species.age_grid_yr()
    if not ages:
        raise ValueError("species age grid must not be empty")
    mortality = np.asarray([species.mortality_per_yr(age) for age in ages], dtype=float)
    if not np.isfinite(mortality).all() or (mortality < 0.0).any():
        raise ValueError("mortality must be finite and nonnegative in yr^-1")
    survival = np.exp(-mortality * DT_YR)
    if not np.isfinite(survival).all() or (survival < 0.0).any() or (survival > 1.0).any():
        raise ValueError("invalid monthly survival")
    return survival


def stable_initial_shape(species: SpeciesParameters) -> np.ndarray:
    """Return the finite stationary age shape induced by unit monthly entry."""
    survival = _monthly_survival(species)
    if survival.size < 2:
        raise ValueError("the age grid must contain a regular cell and a plus group")
    shape = np.empty(survival.size, dtype=float)
    shape[0] = 1.0
    for index in range(1, survival.size - 1):
        shape[index] = shape[index - 1] * survival[index - 1]
    terminal_denominator = 1.0 - survival[-1]
    if terminal_denominator <= 0.0:
        raise ValueError("terminal mortality must permit a finite stationary plus group")
    incoming = shape[-2] * survival[-2]
    shape[-1] = incoming / terminal_denominator
    if not np.isfinite(shape).all() or (shape < 0.0).any():
        raise ValueError("stationary shape is not finite and nonnegative")
    return _readonly(shape)


def build_pde_state_basis(
    species: SpeciesParameters,
    observation_month_offsets: tuple[int, ...],
    knot_indices: tuple[int, ...],
) -> np.ndarray:
    """Build coefficient-wise states on observed months of a complete process grid.

    The returned axes are ``(observed_month, biological_age_cell, coefficient)``.
    The final age cell is the terminal plus group whose coordinate is its lower
    bound, not a known within-group mean age.
    """
    try:
        observed_offsets = tuple(observation_month_offsets)
    except TypeError as exc:
        raise ValueError("observation_month_offsets must be a sequence") from exc
    if len(observed_offsets) < 2:
        raise ValueError("at least two observation month offsets are required")
    if any(isinstance(offset, bool) or not isinstance(offset, Integral) for offset in observed_offsets):
        raise ValueError("observation month offsets must be integers")
    observed_offsets = tuple(int(offset) for offset in observed_offsets)
    if observed_offsets[0] != 0:
        raise ValueError("the first observation month offset must be zero")
    if any(right <= left for left, right in zip(observed_offsets, observed_offsets[1:])):
        raise ValueError("observation month offsets must be strictly increasing")

    n_transitions = observed_offsets[-1]
    entry_basis = piecewise_linear_basis(n_transitions, knot_indices)
    survival = _monthly_survival(species)
    if survival.size < 2:
        raise ValueError("the age grid must contain a regular cell and a plus group")
    n_coefficients = 1 + entry_basis.shape[1]
    states = np.zeros((n_coefficients, survival.size), dtype=float)
    states[0] = stable_initial_shape(species)
    state_basis = np.empty(
        (len(observed_offsets), survival.size, n_coefficients),
        dtype=float,
    )

    observation_index = 0
    for process_offset in range(n_transitions + 1):
        if process_offset == observed_offsets[observation_index]:
            state_basis[observation_index] = states.T
            observation_index += 1
        if process_offset == n_transitions:
            continue
        next_states = np.zeros_like(states)
        next_states[:, 1:] = states[:, :-1] * survival[:-1]
        next_states[:, -1] += states[:, -1] * survival[-1]
        next_states[1:, 0] = entry_basis[process_offset]
        states = next_states

    if not np.isfinite(state_basis).all() or (state_basis < 0.0).any():
        raise ValueError("PDE state basis is not finite and nonnegative")
    return _readonly(state_basis)


def build_pde_design(
    species: SpeciesParameters,
    kernel: np.ndarray,
    observation_month_offsets: tuple[int, ...],
    knot_indices: tuple[int, ...],
) -> np.ndarray:
    """Build observations on a complete monthly process grid.

    ``observation_month_offsets`` are zero-based calendar-month offsets from
    the first observation.  The PDE advances across every intervening month;
    rows are emitted only for the listed observation offsets.
    """
    try:
        observed_offsets = tuple(observation_month_offsets)
    except TypeError as exc:
        raise ValueError("observation_month_offsets must be a sequence") from exc
    if len(observed_offsets) < 2:
        raise ValueError("at least two observation month offsets are required")
    if any(isinstance(offset, bool) or not isinstance(offset, Integral) for offset in observed_offsets):
        raise ValueError("observation month offsets must be integers")
    observed_offsets = tuple(int(offset) for offset in observed_offsets)
    if observed_offsets[0] != 0:
        raise ValueError("the first observation month offset must be zero")
    if any(right <= left for left, right in zip(observed_offsets, observed_offsets[1:])):
        raise ValueError("observation month offsets must be strictly increasing")
    n_transitions = observed_offsets[-1]
    ages = species.age_grid_yr()
    try:
        observation_kernel = np.asarray(kernel, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("kernel must be a numeric matrix") from exc
    if observation_kernel.ndim != 2 or observation_kernel.shape[0] != len(ages):
        raise ValueError("kernel shape must be (n_age_cells, n_bins)")
    if observation_kernel.shape[1] == 0:
        raise ValueError("kernel must contain at least one observation bin")
    if not np.isfinite(observation_kernel).all() or (observation_kernel < 0.0).any():
        raise ValueError("kernel entries must be finite and nonnegative")

    state_basis = build_pde_state_basis(species, observed_offsets, knot_indices)
    n_bins = observation_kernel.shape[1]
    n_coefficients = state_basis.shape[2]
    design = np.empty((len(observed_offsets) * n_bins, n_coefficients), dtype=float)
    for observation_index in range(len(observed_offsets)):
        row_start = observation_index * n_bins
        for coefficient_index in range(n_coefficients):
            # one biological-age state vector is observed at a time.  A single
            # matrix contraction over every coefficient is algebraically equal
            design[
                row_start : row_start + n_bins,
                coefficient_index,
            ] = DT_YR * (
                np.ascontiguousarray(
                    state_basis[observation_index, :, coefficient_index],
                    dtype=float,
                )
                @ observation_kernel
            )

    if not np.isfinite(design).all() or (design < 0.0).any():
        raise ValueError("PDE design is not finite and nonnegative")
    return _readonly(design)


def standardize_design(
    design: np.ndarray,
    observations: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Scale each bin by its median positive observation in month-major order."""
    try:
        matrix = np.asarray(design, dtype=float)
        observed = np.asarray(observations, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("design and observations must be numeric arrays") from exc
    if matrix.ndim != 2:
        raise ValueError("design must be two-dimensional")
    if observed.ndim != 2:
        raise ValueError("observations must have shape (n_months, n_bins)")
    if matrix.shape[0] != observed.size:
        raise ValueError("design rows must match flattened observations")
    if not np.isfinite(matrix).all():
        raise ValueError("design must be finite")
    if not np.isfinite(observed).all() or (observed < 0.0).any():
        raise ValueError("observations must be finite and nonnegative")

    scales = np.empty(observed.shape[1], dtype=float)
    for bin_index in range(observed.shape[1]):
        positive = observed[:, bin_index][observed[:, bin_index] > 0.0]
        if positive.size == 0:
            raise ValueError("each bin must contain at least one positive observation")
        scales[bin_index] = float(np.median(positive))
    row_scales = np.tile(scales, observed.shape[0])
    standardized_design = np.array(matrix / row_scales[:, None], copy=True)
    standardized_observations = np.array(observed.ravel(order="C") / row_scales, copy=True)
    return (
        _readonly(standardized_design),
        _readonly(standardized_observations),
        _readonly(scales),
    )


def _projected_gradient(coefficients: np.ndarray, gradient: np.ndarray) -> np.ndarray:
    return np.where((coefficients <= 0.0) & (gradient > 0.0), 0.0, gradient)


def projected_gradient_nnls(
    design: np.ndarray,
    observations: np.ndarray,
    max_seconds: float,
    max_iterations: int,
    tolerance: float,
    *,
    column_scaling: str = "none",
) -> FitResult:
    """Minimize ``0.5 * ||design @ x - observations||^2`` for ``x >= 0``."""
    started = time.perf_counter()
    try:
        matrix = np.asarray(design, dtype=float)
        observed = np.asarray(observations, dtype=float)
        seconds = float(max_seconds)
        convergence_tolerance = float(tolerance)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("solver inputs and controls must be numeric") from exc
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("design must be a nonempty two-dimensional matrix")
    if observed.ndim != 1 or observed.shape[0] != matrix.shape[0]:
        raise ValueError("observations must be a vector matching the design rows")
    if not np.isfinite(matrix).all() or not np.isfinite(observed).all():
        raise ValueError("design and observations must be finite")
    if not math.isfinite(seconds) or seconds < 0.0:
        raise ValueError("max_seconds must be finite and nonnegative")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")
    if not math.isfinite(convergence_tolerance) or convergence_tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    if column_scaling not in {"none", "l2"}:
        raise ValueError("column_scaling must be 'none' or 'l2'")
    if not np.any(matrix != 0.0):
        raise ValueError("design must have a finite nonzero spectral norm")

    def time_limit_result() -> FitResult:
        coefficients = _readonly(np.zeros(matrix.shape[1], dtype=float))
        return FitResult(
            coefficients=coefficients,
            converged=False,
            iterations=0,
            objective=0.5 * float(observed @ observed),
            projected_gradient_inf=float("nan"),
            elapsed_seconds=time.perf_counter() - started,
            stop_reason="time_limit",
        )

    if seconds == 0.0 or time.perf_counter() - started >= seconds:
        return time_limit_result()

    column_norms: np.ndarray | None = None
    working_matrix = matrix
    if column_scaling == "l2":
        column_norms = np.linalg.norm(matrix, axis=0)
        if not np.isfinite(column_norms).all() or (column_norms <= 0.0).any():
            raise ValueError("l2 column scaling requires finite positive column norms")
        working_matrix = matrix / column_norms

    spectral_norm = float(np.linalg.norm(working_matrix, ord=2))
    if not math.isfinite(spectral_norm) or spectral_norm <= 0.0:
        raise ValueError("design must have a finite nonzero spectral norm")
    if time.perf_counter() - started >= seconds:
        return time_limit_result()
    step = 1.0 / (spectral_norm * spectral_norm)
    reference_gradient = matrix.T @ observed
    convergence_scale = convergence_tolerance * max(
        1.0,
        float(np.linalg.norm(reference_gradient, ord=np.inf)),
    )

    coefficients = np.zeros(matrix.shape[1], dtype=float)
    iterations = 0
    converged = False
    stop_reason = "max_iterations"
    for _ in range(int(max_iterations)):
        if time.perf_counter() - started >= seconds:
            stop_reason = "time_limit"
            break
        residual = working_matrix @ coefficients - observed
        gradient = working_matrix.T @ residual
        if column_norms is None:
            original_coefficients = coefficients
            original_gradient = gradient
        else:
            original_coefficients = coefficients / column_norms
            original_gradient = column_norms * gradient
        projected = _projected_gradient(original_coefficients, original_gradient)
        projected_inf = float(np.linalg.norm(projected, ord=np.inf))
        if projected_inf <= convergence_scale:
            converged = True
            stop_reason = "converged"
            break
        coefficients = np.maximum(0.0, coefficients - step * gradient)
        iterations += 1

    elapsed = time.perf_counter() - started
    if column_norms is not None:
        coefficients = coefficients / column_norms
    residual = matrix @ coefficients - observed
    gradient = matrix.T @ residual
    projected_inf = float(np.linalg.norm(_projected_gradient(coefficients, gradient), ord=np.inf))
    objective = 0.5 * float(residual @ residual)
    if stop_reason == "max_iterations" and projected_inf <= convergence_scale:
        converged = True
        stop_reason = "converged"
    coefficients = _readonly(np.array(coefficients, copy=True))
    return FitResult(
        coefficients=coefficients,
        converged=converged,
        iterations=iterations,
        objective=objective,
        projected_gradient_inf=projected_inf,
        elapsed_seconds=elapsed,
        stop_reason=stop_reason,
    )
