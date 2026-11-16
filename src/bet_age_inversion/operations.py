"""Operation-level fish-equivalent allocation and CPUE from monthly age states."""
from __future__ import annotations
import math
from typing import Any, Sequence
import numpy as np
from .parameters import GG_TO_RW

BET_BIN_COLUMNS = ("wBET15D", "wBET15U", "wBET25", "wBET40")

def project_operations(
    bundle: dict[str, Any],
    states: np.ndarray,
    frame: Any,
    months: Sequence[tuple[int, int] | str],
) -> dict[str, np.ndarray]:
    """Project four GG catch bins to model-conditioned fish-age CPUE.

    ``states`` must align row-for-row with ``months``. Blank bin cells mean
    zero catch, weights are converted from GG to RW once, and CPUE is fish
    equivalents per 1,000 hooks.
    """
    state_array = np.asarray(states, dtype=float)
    ages = np.asarray(bundle["ages"], dtype=float)
    probabilities = np.asarray(bundle["bin_probabilities"], dtype=float)
    moments = np.asarray(bundle["kernel"], dtype=float)
    if state_array.ndim != 2 or state_array.shape[1] != ages.size:
        raise ValueError("states must have shape (months, ages)")
    if len(months) != state_array.shape[0]:
        raise ValueError("months must align with the state rows")
    if probabilities.shape != moments.shape or probabilities.shape != (ages.size, len(BET_BIN_COLUMNS)):
        raise ValueError("bundle kernels must have shape (ages, 4 bins)")
    if not all(np.isfinite(x).all() for x in (state_array, ages, probabilities, moments)):
        raise ValueError("states, ages, and kernels must be finite")
    if any((x < 0.0).any() for x in (state_array, probabilities, moments)):
        raise ValueError("states and kernels must be nonnegative")

    def month_key(value: tuple[int, int] | str) -> str:
        if isinstance(value, str):
            return value
        return f"{int(value[0]):04d}-{int(value[1]):02d}"

    month_to_index = {month_key(value): index for index, value in enumerate(months)}
    if len(month_to_index) != len(months):
        raise ValueError("months contain duplicate keys")
    required = ("Year", "Month", "Hooks", *BET_BIN_COLUMNS)
    try:
        columns = {name: np.asarray(frame[name], dtype=object) for name in required}
    except KeyError as exc:
        raise ValueError(f"missing operation column: {exc.args[0]}") from exc
    if any(col.ndim != 1 for col in columns.values()):
        raise ValueError("operation columns must be one-dimensional")
    if len({len(col) for col in columns.values()}) != 1:
        raise ValueError("operation columns must have the same length")

    def date_column(name: str, maximum: int) -> np.ndarray:
        try:
            values = columns[name].astype(float)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{name} must be an integer from 1 to {maximum}") from exc
        if (not np.isfinite(values).all() or (values < 1).any()
                or (values > maximum).any() or (values != np.floor(values)).any()):
            raise ValueError(f"{name} must be an integer from 1 to {maximum}")
        return values.astype(int)

    years = date_column("Year", 9999)
    month_numbers = date_column("Month", 12)
    try:
        hooks = columns["Hooks"].astype(float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Hooks must be finite and positive") from exc

    def catch_column(name: str) -> np.ndarray:
        values = np.zeros(len(columns[name]), dtype=float)
        for index, value in enumerate(columns[name]):
            if value is None or (isinstance(value, str) and not value.strip()):
                continue
            try:
                numeric = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"{name} must contain numeric weights or blanks") from exc
            values[index] = 0.0 if math.isnan(numeric) else numeric
        return values

    gg_weights = np.column_stack([catch_column(name) for name in BET_BIN_COLUMNS])
    if not np.isfinite(hooks).all() or (hooks <= 0.0).any():
        raise ValueError("Hooks must be finite and positive")
    if not np.isfinite(gg_weights).all() or (gg_weights < 0.0).any():
        raise ValueError("BET catch bins must be finite and nonnegative")
    keys = [f"{year:04d}-{month:02d}" for year, month in zip(years, month_numbers, strict=True)]
    try:
        state_indices = np.asarray([month_to_index[key] for key in keys], dtype=int)
    except KeyError as exc:
        raise ValueError(f"operation month has no state: {exc.args[0]}") from exc

    n = len(keys)
    age_counts = np.zeros((n, ages.size), dtype=float)
    weights_rw = gg_weights * GG_TO_RW
    positive = weights_rw.sum(axis=1) > 0.0
    for row in np.flatnonzero(positive):
        month_state = state_array[state_indices[row]]
        state_mass = month_state.sum()
        if state_mass <= 0.0:
            raise ValueError(f"positive-catch operation month has no positive state: {keys[row]}")
        normalized_state = month_state / state_mass
        for bin_index, biomass in enumerate(weights_rw[row]):
            if biomass == 0.0:
                continue
            denominator = float(normalized_state @ moments[:, bin_index])
            if not math.isfinite(denominator) or denominator <= 0.0:
                raise ValueError("positive catch bin has no first-moment support")
            age_counts[row] += (
                biomass * normalized_state * probabilities[:, bin_index] / denominator
            )

    fish_count = age_counts.sum(axis=1)
    expected_age = np.full(n, np.nan, dtype=float)
    expected_age[positive] = (age_counts[positive] @ ages) / fish_count[positive]
    cpue_by_age = age_counts * (1000.0 / hooks[:, None])
    arrays = (age_counts, fish_count, cpue_by_age)
    if not all(np.isfinite(array).all() for array in arrays):
        raise ValueError("operation projection produced non-finite values")
    if any((array < 0.0).any() for array in arrays):
        raise ValueError("operation projection produced negative values")
    if np.any(positive & (np.abs(age_counts.sum(axis=1) - fish_count) > 1e-10)):
        raise ValueError("fish count does not equal summed age counts")
    if not np.allclose(cpue_by_age.sum(axis=1), fish_count * 1000.0 / hooks, rtol=1e-12, atol=1e-12):
        raise ValueError("age CPUE does not conserve total fish-equivalent CPUE")
    return {
        "age_counts": age_counts,
        "fish_count": fish_count,
        "expected_age": expected_age,
        "cpue_by_age": cpue_by_age,
        "positive": positive,
    }

def group_cpue(projected: dict[str, np.ndarray], ages: np.ndarray, cutpoints: Sequence[float]) -> np.ndarray:
    """Sum age-specific fish-equivalent CPUE into ordered half-open groups."""
    age_grid = np.asarray(ages, dtype=float)
    cpue = np.asarray(projected["cpue_by_age"], dtype=float)
    cuts = np.asarray(cutpoints, dtype=float)
    if age_grid.ndim != 1 or cpue.ndim != 2 or cpue.shape[1] != age_grid.size:
        raise ValueError("CPUE and age grid shapes are incompatible")
    if cuts.ndim != 1 or not np.isfinite(cuts).all() or (np.diff(cuts) <= 0.0).any():
        raise ValueError("cutpoints must be finite and strictly increasing")
    if not np.isfinite(age_grid).all() or not np.isfinite(cpue).all() or (cpue < 0.0).any():
        raise ValueError("ages and CPUE must be finite; CPUE must be nonnegative")
    boundaries = np.concatenate(([-np.inf], cuts, [np.inf]))
    memberships = np.vstack([
        (age_grid >= boundaries[i]) & (age_grid < boundaries[i + 1])
        for i in range(boundaries.size - 1)
    ])
    if not np.all(memberships.sum(axis=0) == 1):
        raise ValueError("age groups do not partition the registered age grid")
    grouped = cpue @ memberships.T
    if not np.allclose(grouped.sum(axis=1), cpue.sum(axis=1), rtol=1e-12, atol=1e-12):
        raise ValueError("grouped CPUE does not conserve total CPUE")
    return grouped
