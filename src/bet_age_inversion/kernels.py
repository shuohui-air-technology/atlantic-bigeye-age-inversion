"""Age-conditional weight-bin probabilities and first weight moments."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np

from .parameters import SpeciesParameters


@dataclass(frozen=True)
class KernelAtAge:
    """Observed-bin probabilities and unconditional first weight moments at one age."""

    age_yr: float
    mean_length_cm: float
    sd_length_cm: float
    probabilities: tuple[float, ...]
    probability_sum_observed_bins: float
    weight_moments_kg_rw: tuple[float, ...]


def normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def distribution_pdf(kind: str, mean: float, sd: float) -> tuple[Callable[[float], float], float]:
    """Return a positive-length density and a numerically adequate upper limit."""
    variance = sd * sd
    if kind == "normal":
        lower_mass = normal_cdf(-mean / sd)
        norm = 1.0 - lower_mass

        def pdf(value: float) -> float:
            if value < 0.0:
                return 0.0
            z = (value - mean) / sd
            return math.exp(-0.5 * z * z) / (sd * math.sqrt(2.0 * math.pi) * norm)

        upper = mean + 9.0 * sd
    elif kind == "lognormal":
        sigma2 = math.log1p(variance / (mean * mean))
        sigma = math.sqrt(sigma2)
        mu = math.log(mean) - 0.5 * sigma2

        def pdf(value: float) -> float:
            if value <= 0.0:
                return 0.0
            z = (math.log(value) - mu) / sigma
            return math.exp(-0.5 * z * z) / (value * sigma * math.sqrt(2.0 * math.pi))

        upper = math.exp(mu + 9.0 * sigma)
    elif kind == "gamma":
        shape = mean * mean / variance
        scale = variance / mean
        log_norm = math.lgamma(shape) + shape * math.log(scale)

        def pdf(value: float) -> float:
            if value <= 0.0:
                return 0.0
            log_pdf = (shape - 1.0) * math.log(value) - value / scale - log_norm
            return math.exp(log_pdf)

        upper = mean + 15.0 * sd
    else:
        raise ValueError(f"Unknown distribution: {kind}")
    return pdf, upper


def simpson(func: Callable[[float], float], lo: float, hi: float, intervals: int = 1200) -> float:
    if hi <= lo:
        return 0.0
    if intervals % 2:
        intervals += 1
    width = (hi - lo) / intervals
    total = func(lo) + func(hi)
    total += 4.0 * sum(func(lo + i * width) for i in range(1, intervals, 2))
    total += 2.0 * sum(func(lo + i * width) for i in range(2, intervals, 2))
    return total * width / 3.0


def bin_kernel(species: SpeciesParameters, age_yr: float, kind: str) -> KernelAtAge:
    """Calculate full-precision observed-bin probability and weight-moment terms."""
    mean = species.mean_length_cm(age_yr)
    sd = species.sd_length_cm(age_yr)
    if not (mean > 0.0 and sd > 0.0):
        raise AssertionError(f"Invalid length moments for {species.code} at age {age_yr}")
    pdf, upper = distribution_pdf(kind, mean, sd)
    total_probability = simpson(pdf, 0.0, upper)
    if abs(total_probability - 1.0) > 2e-4:
        raise AssertionError((species.code, age_yr, kind, total_probability))

    probabilities: list[float] = []
    weight_moments: list[float] = []
    for lo_weight, hi_weight in species.model_bins_kg_rw:
        lo_length = species.length_cm_from_weight_kg_rw(lo_weight)
        hi_length = upper if hi_weight is None else species.length_cm_from_weight_kg_rw(hi_weight)
        hi_length = min(upper, hi_length)
        probability = simpson(pdf, lo_length, hi_length)
        moment = simpson(lambda length: species.weight_kg_rw(length) * pdf(length), lo_length, hi_length)
        if probability < -1e-12 or moment < -1e-12:
            raise AssertionError("Negative kernel entry")
        probabilities.append(max(0.0, probability))
        weight_moments.append(max(0.0, moment))

    probability_sum = sum(probabilities)
    # The BET observation bins partition the positive weight axis.
    if abs(probability_sum - 1.0) > 3e-4:
        raise AssertionError((species.code, age_yr, kind, probability_sum))
    return KernelAtAge(
        age_yr=age_yr,
        mean_length_cm=mean,
        sd_length_cm=sd,
        probabilities=tuple(probabilities),
        probability_sum_observed_bins=probability_sum,
        weight_moments_kg_rw=tuple(weight_moments),
    )


def kernel_at_age(species: SpeciesParameters, age_yr: float, kind: str) -> KernelAtAge:
    """Return the first-moment observation kernel for one registered species age."""
    return bin_kernel(species, age_yr, kind)


def build_weight_moment_matrix(species: SpeciesParameters, kind: str = "normal") -> np.ndarray:
    """Build age-by-bin unconditional weight moments in registered grid/bin order."""
    return np.asarray(
        [kernel_at_age(species, age_yr, kind).weight_moments_kg_rw for age_yr in species.age_grid_yr()],
        dtype=float,
    )
