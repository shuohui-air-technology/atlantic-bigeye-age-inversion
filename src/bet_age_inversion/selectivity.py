"""Published logistic size-selectivity curves for Atlantic bigeye tuna."""

import numpy as np
from scipy.special import expit, ndtr


# ICCAT SCRS/2025/151, Table 5. Columns are L50 and L95 - L50 in cm.
SELECTIVITY_CURVES = {"tropical": (120.92, 23.36), "southern": (108.92, 24.03)}


def length_weighted_kernel(sp, curve: str) -> tuple[np.ndarray, np.ndarray]:
    """Integrate s(length) times the positive-length density in each bin."""
    l50, d95 = SELECTIVITY_CURVES[curve]
    ages = np.array(sp.age_grid_yr())
    mu = np.array([sp.mean_length_cm(a) for a in ages])
    sd = np.array([sp.sd_length_cm(a) for a in ages])
    upper = mu + 9 * sd
    lo = np.array([sp.length_cm_from_weight_kg_rw(x[0]) for x in sp.model_bins_kg_rw])
    hi = np.array([np.inf if x[1] is None else sp.length_cm_from_weight_kg_rw(x[1]) for x in sp.model_bins_kg_rw])
    weights = np.ones(1201)
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    P = np.zeros((len(ages), 4))
    K = P.copy()
    t = np.linspace(0, 1, 1201)
    for b in range(4):
        high = np.minimum(upper, hi[b])
        width = np.maximum(high - lo[b], 0)
        L = lo[b] + width[:, None] * t
        z = (L - mu[:, None]) / sd[:, None]
        pdf = np.exp(-0.5 * z * z) / (sd[:, None] * np.sqrt(2 * np.pi) * (1 - ndtr(-mu / sd))[:, None])
        pdf *= expit(np.log(19) * (L - l50) / d95)
        P[:, b] = np.sum(pdf * weights, axis=1) * width / 3600
        K[:, b] = np.sum(pdf * weights * sp.lw_a * L ** sp.lw_b, axis=1) * width / 3600
    return P, K

