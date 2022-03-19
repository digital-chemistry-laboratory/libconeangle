"""Cone angle code."""
from __future__ import annotations

from ctypes import byref, c_double, c_int, create_string_buffer

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
        ValueError: If libconeangle exits with non-zero error code or if input validation fails.
    """
    # Set up input
    coordinates: Array2D = np.ascontiguousarray(coordinates)
    radii: Array1D = np.ascontiguousarray(radii, dtype=np.float64)

    if len(coordinates) != len(radii):
        raise ValueError(
            f"Length of coordinates, {len(coordinates)}, "
            f"and radii, {len(radii)}, must be the same."
        )

    cone_angle_ = c_double()
    axis = np.zeros(3)
    tangent_atoms: Array1D = np.zeros(3, dtype=np.int32)
    stat_ = c_int()
    errmsg_ = create_string_buffer(100)
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
        errmsg_,
    )

    # Check for non-zero exit code
    stat = int(stat_.value)
    if stat != 0:
        raise ValueError(
            f"libconeangle exited with non-zero exit code: {stat}. "
            f"Error message: {errmsg_.value.decode()}"
        )

    # Take out results
    cone_angle = float(cone_angle_.value)
    tangent_atoms = list(tangent_atoms[tangent_atoms.nonzero()])

    return cone_angle, axis, tangent_atoms
