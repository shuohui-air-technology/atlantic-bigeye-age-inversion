"""Biological and observation parameters. Length is cm, round weight is kg per fish, age is yr, and mortality is per yr."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable


GG_TO_RW = 1.09
DT_YR = 1.0 / 12.0


@dataclass(frozen=True)
class SpeciesParameters:
    code: str
    min_age_yr: float
    max_age_yr: float
    coarse_age_groups_yr: tuple[tuple[float, float | None], ...]
    source_bins_kg_gg: tuple[tuple[float, float | None], ...]
    lw_a: float
    lw_b: float
    mean_length_cm: Callable[[float], float]
    sd_length_cm: Callable[[float], float]
    mortality_per_yr: Callable[[float], float]

    @property
    def model_bins_kg_rw(self) -> tuple[tuple[float, float | None], ...]:
        return tuple(
            (lo * GG_TO_RW, None if hi is None else hi * GG_TO_RW)
            for lo, hi in self.source_bins_kg_gg
        )

    def weight_kg_rw(self, length_cm: float) -> float:
        return self.lw_a * max(length_cm, 0.0) ** self.lw_b

    def length_cm_from_weight_kg_rw(self, weight_kg: float) -> float:
        if weight_kg <= 0.0:
            return 0.0
        return (weight_kg / self.lw_a) ** (1.0 / self.lw_b)

    def age_grid_yr(self) -> tuple[float, ...]:
        count = int(math.floor((self.max_age_yr - self.min_age_yr) / DT_YR + 1e-9))
        ages = [self.min_age_yr + i * DT_YR for i in range(count + 1)]
        if ages[-1] < self.max_age_yr - 1e-9:
            ages.append(self.max_age_yr)
        return tuple(ages)


def bet_length_cm(age_yr: float) -> float:
    return 161.21 * (1.0 - math.exp(-0.392 * (age_yr + 0.239)))


_BET_CALIBRATION_AGES = tuple(4.0 + i * DT_YR for i in range(73))
_BET_C = 0.32 / (
    sum(1.0 / bet_length_cm(age) for age in _BET_CALIBRATION_AGES)
    / len(_BET_CALIBRATION_AGES)
)


def bet_mortality(age_yr: float) -> float:
    return _BET_C / bet_length_cm(age_yr)


BET = SpeciesParameters(
    code="BET",
    min_age_yr=0.0,
    max_age_yr=17.0,
    coarse_age_groups_yr=((0.0, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, None)),
    source_bins_kg_gg=((0.0, 15.0), (15.0, 25.0), (25.0, 40.0), (40.0, None)),
    lw_a=2.396e-5,
    lw_b=2.9774,
    mean_length_cm=bet_length_cm,
    sd_length_cm=lambda _age: 8.75,
    mortality_per_yr=bet_mortality,
)
