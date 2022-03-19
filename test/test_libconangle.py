"""Test libconeangle."""

import numpy as np
from numpy.testing import assert_allclose
from numpy.typing import NDArray
import pytest

from libconeangle import cone_angle


def test_PdCO():
    """Test PdCO with cone tangent to one atom."""
    coordinates: NDArray[np.float64] = np.array(
        [[0.0, 0.0, -0.52], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]]
    )
    radii: NDArray[np.float64] = np.array([2.1, 1.7, 1.52])
    c_angle, axis, tangent_atoms = cone_angle(coordinates, radii, 0)

    assert_allclose(c_angle, 96.423, atol=0.001)
    assert_allclose(axis, np.array([0.0, 0.0, 1.0]), atol=0.0001)
    assert set(tangent_atoms) == set([1])


def test_PdCO_close():
    """Test PdCO with one atom within vdW distance of Pd atom."""
    coordinates: NDArray[np.float64] = np.array(
        [[0.0, 0.0, 0.1], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]]
    )  # Pd atom too close
    radii: NDArray[np.float64] = np.array([2.1, 1.7, 1.52])
    with pytest.raises(ValueError):
        cone_angle(coordinates, radii, 0)


def test_PdCO_bounds():
    """Test PdCO with metal index out of bounds."""
    coordinates: NDArray[np.float64] = np.array(
        [[0.0, 0.0, 0.1], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]]
    )  # Pd atom too close
    radii: NDArray[np.float64] = np.array([2.1, 1.7, 1.52])
    with pytest.raises(ValueError):
        cone_angle(coordinates, radii, -1)


def test_PdCO_mismatch():
    """Test PdCO with metal index out of bounds."""
    coordinates: NDArray[np.float64] = np.array(
        [[0.0, 0.0, 0.1], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]]
    )  # Pd atom too close
    radii: NDArray[np.float64] = np.array([2.1, 1.7])
    with pytest.raises(ValueError):
        cone_angle(coordinates, radii, 0)


def test_Pdbpy():
    """Test Pd(bpy) with cone tangent to two atoms."""
    coordinates: NDArray[np.float64] = np.array(
        [
            [-1.899000e-03, -1.846004e00, -2.007000e-03],
            [-1.343609e00, -8.359200e-02, 4.197960e-01],
            [1.345137e00, -8.591400e-02, -4.193200e-01],
            [-7.399130e-01, 1.062540e00, 8.798000e-03],
            [-1.473982e00, 2.183787e00, -3.961990e-01],
            [-2.865881e00, 2.130461e00, -3.962120e-01],
            [-3.487476e00, 9.490540e-01, 1.374000e-02],
            [-2.690562e00, -1.192190e-01, 4.151720e-01],
            [-9.542920e-01, 3.073408e00, -7.377720e-01],
            [-3.451265e00, 2.985788e00, -7.207940e-01],
            [-4.569174e00, 8.579540e-01, 3.831000e-02],
            [-3.135734e00, -1.043925e00, 7.700380e-01],
            [2.691916e00, -1.229550e-01, -4.128390e-01],
            [3.489163e00, 9.451910e-01, -1.169700e-02],
            [2.868283e00, 2.127761e00, 3.960000e-01],
            [1.476591e00, 2.182376e00, 3.941950e-01],
            [7.424670e-01, 1.061320e00, -1.040300e-02],
            [3.136172e00, -1.048666e00, -7.660710e-01],
            [4.570699e00, 8.531320e-01, -3.487800e-02],
            [3.453586e00, 2.983188e00, 7.202120e-01],
            [9.563590e-01, 3.072406e00, 7.338950e-01],
        ]
    )
    radii: NDArray[np.float64] = np.array(
        [
            2.1,
            1.55,
            1.55,
            1.7,
            1.7,
            1.7,
            1.7,
            1.7,
            1.1,
            1.1,
            1.1,
            1.1,
            1.7,
            1.7,
            1.7,
            1.7,
            1.7,
            1.1,
            1.1,
            1.1,
            1.1,
        ]
    )
    c_angle, axis, tangent_atoms = cone_angle(coordinates, radii, 0)

    assert_allclose(c_angle, 190.799, atol=0.001)
    assert_allclose(axis, np.array([0.00200079, 0.99998474, 0.00514848]), atol=0.0001)
    assert set(tangent_atoms) == set([11, 17])


def test_PdPMe3():
    """Test PdPMe3 with cone tangent to three atoms."""
    coordinates: NDArray[np.float64] = np.array(
        [
            [-7.433970e-01, -9.920000e-04, -1.518770e00],
            [1.005896e00, 2.102000e-03, -1.636180e-01],
            [2.689042e00, -1.142000e-03, -9.464240e-01],
            [2.797104e00, -8.888730e-01, -1.576886e00],
            [3.486761e00, 3.481000e-03, -1.921520e-01],
            [2.795646e00, 8.800330e-01, -1.586258e00],
            [1.164898e00, 1.431506e00, 1.010324e00],
            [1.222204e00, 2.365419e00, 4.431760e-01],
            [2.057780e00, 1.338752e00, 1.642300e00],
            [2.773700e-01, 1.479634e00, 1.648431e00],
            [1.165320e00, -1.429939e00, 1.006886e00],
            [1.222850e00, -2.362601e00, 4.376840e-01],
            [2.776990e-01, -1.479572e00, 1.644733e00],
            [2.058112e00, -1.338384e00, 1.639153e00],
        ]
    )
    radii: NDArray[np.float64] = np.array(
        [2.1, 1.8, 1.7, 1.1, 1.1, 1.1, 1.7, 1.1, 1.1, 1.1, 1.7, 1.1, 1.1, 1.1]
    )
    c_angle, axis, tangent_atoms = cone_angle(coordinates, radii, 0)

    assert_allclose(c_angle, 117.110, atol=0.001)
    assert_allclose(
        axis, np.array([7.89524317e-01, 5.22544038e-05, 6.13719276e-01]), atol=0.0001
    )
    assert set(tangent_atoms) == set([5, 9, 12])
