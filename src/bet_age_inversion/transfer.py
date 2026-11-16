"""Symmetric transfer across specified adjacent weight-bin boundaries."""
from __future__ import annotations

def adjacent_misclassification_matrix(
    n_bins: int,
    epsilon: float,
    boundary_indices: tuple[int, ...],
) -> tuple[tuple[float, ...], ...]:
    """Return P(recorded bin | true bin) for symmetric adjacent leakage.

    Boundary index ``j`` denotes the boundary between bins ``j`` and ``j + 1``.
    Each selected boundary moves ``epsilon`` of either adjacent true bin to the
    other recorded bin. Columns therefore represent true bins and sum to one.
    """
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2")
    if not 0.0 <= epsilon <= 1.0:
        raise ValueError("epsilon must be between 0 and 1")
    boundaries = tuple(dict.fromkeys(boundary_indices))
    if any(index < 0 or index >= n_bins - 1 for index in boundaries):
        raise ValueError("boundary index is outside the bin range")

    matrix = [[0.0 for _ in range(n_bins)] for _ in range(n_bins)]
    for index in range(n_bins):
        matrix[index][index] = 1.0
    for left in boundaries:
        right = left + 1
        matrix[left][left] -= epsilon
        matrix[right][left] += epsilon
        matrix[right][right] -= epsilon
        matrix[left][right] += epsilon

    if any(matrix[index][index] < -1e-12 for index in range(n_bins)):
        raise ValueError("epsilon is too large for the selected adjacent boundaries")
    return tuple(tuple(round(value, 15) for value in row) for row in matrix)
