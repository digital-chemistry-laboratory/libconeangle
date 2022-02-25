"""Cone angle code."""
from __future__ import annotations

from ctypes import byref, c_double, c_int

import numpy as np

from libconeangle.lib import lib
from libconeangle.typing import Array1D, Array2D, ArrayLike1D, ArrayLike2D


def cone_angle(
    coordinates: ArrayLike2D, radii: ArrayLike1D, idx_metal: int
) -> tuple[float, Array1D, list[int]]:
    """Calculates cone angle, cone axis and tangent atoms.

    Args:
        coordinates: Coordinates as n x 3 matrix (Å)
        radii: vdW radii (Å)
        idx_metal: Index of metal atom (1-indexed)

    Returns:
        cone_angle: Cone angle (degrees)
        axis: Cone normal axis (Å)
        tangent_atoms: Indices of atoms tangent to cone

    Raises:
        ValueError: If libconeangle exits with non-zero error code
    """
    # Set up input
    coordinates: Array2D = np.ascontiguousarray(coordinates)
    radii: Array1D = np.ascontiguousarray(radii, dtype=np.float64)
    cone_angle_ = c_double()
    axis = np.zeros(3)
    tangent_atoms: Array1D = np.zeros(3, dtype=np.int32)
    stat_ = c_int()
    n_atoms = c_int(coordinates.shape[0])

    # Run libconeangle
    lib.cone_angle(
        n_atoms,
        coordinates,
        radii,
        c_int(idx_metal),
        byref(cone_angle_),
        axis,
        tangent_atoms,
        byref(stat_),
    )

    # Check for non-zero exit code
    stat = int(stat_.value)
    if stat != 0:
        raise ValueError("libconeangle exited with non-zero exit code: {stat}")

    # Take out results
    cone_angle = float(cone_angle_.value)
    tangent_atoms = list(tangent_atoms[tangent_atoms.nonzero()])

    return cone_angle, axis, tangent_atoms
