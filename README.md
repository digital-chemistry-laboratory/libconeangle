# libconeangle

Library for calculating exact ligand cone angles according to the recipe of Allen and co-workers.[^1] This library is not meant as a standalone application but is rather meant for integration into other programs. An example is [ᴍᴏʀғᴇᴜs](https://github.com/kjelljorner/morfeus) by the same author. libconeangle is written in Fortran with a C interface.

## Installation

### pip

To install the Python API with the embedded libconeangle shared library you can use pip. You will need a Fortran compiler such as gfortran to compile the shared library but the whole process is automated via [scikit-build](https://github.com/scikit-build/scikit-build).

```shell
pip install git+https://github.com/kjelljorner/libconeangle
```

### conda

Another option to install the Python API that doesn't require a compiler is with conda.

```shell
conda install -c conda-forge libconeangle
```

### cmake

The shared library can be built and installed with cmake. An example worklfow is given below where you need to replace $PREFIX with the desired directory.

```shell
FC=gfortran cmake -B _build -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
cmake --build _build
cmake --install _build
```

### fpm

To use as a dependency in your Fortran project with [fpm](https://github.com/fortran-lang/fpm), add the following to your `fpm.toml`.

```toml
[dependencies]
libconeangle.git = "https://github.com/kjelljorner/libconeangle"
```

## Usage

### Python API

There is only one function: `cone_angle`. An example is given below for PdCO.

> ⚠️ All atoms are zero-index in the Python API.

```python
>>> from libconeangle import cone_angle
>>> coordinates =  np.array([[0.0, 0.0, -0.52], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]])
>>> radii = np.array([2.1, 1.7, 1.52])
>>> index_metal = 0 # Zero-indexed
>>> c_angle, axis, tangent_atoms = cone_angle(coordinates, radii, index_metal)
>>> c_angle
96.4237340645161
>>> axis
array([0., 0., 1.])
>>> tangent_atoms # Also zero-indexed
[1]
```

### Fortran API

The Fortran API exposes the function `cone_angle` with the follow signature.

> ⚠️ All atoms are one-index in the Fortran API.

```fortran
subroutine cone_angle(coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat)
  !! Calculate cone angle, cone axis and tangent atoms
  !> Coordinates (Å)
  real(wp), intent(in) :: coordinates(:, :)
  !> vdW radii (Å)
  real(wp), intent(in) :: radii(:)
  !> Index of metal atom
  integer, intent(in) :: index_metal
  !> Cone angle (degrees)
  real(wp), intent(out) :: alpha
  !> Cone axis (Å)
  real(wp), intent(out) :: axis(3)
  !> Indices of atoms tangent to cone
  integer, intent(out) :: tangent_atoms(3)
  !> Return code
  integer, intent(out) :: stat
  
  ...
  
end subroutine cone_angle
```

Here is one example of how it could be used as given in the demo [program](app/demo.f90).

```fortran
program demo
  use coneangle_main, only: cone_angle

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp) :: coordinates(3, 3), radii(3), alpha, axis(3)
  integer :: tangent_atoms(3), stat
  character(:), allocatable :: errmsg

  coordinates = reshape([0._dp, 0._dp, -0.52_dp, 0._dp, 0._dp, 1.76_dp, 0._dp, 0._dp, 2.86_dp], [3, 3])
  radii = [2.1_dp, 1.7_dp, 1.52_dp]
  call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)
  write (*, *) "Cone angle:", alpha
  write (*, *) "Cone axis:", axis
  write (*, *) "Tangent atoms:", tangent_atoms
end program demo
```

The `tangent_atoms` array has three elements. In the case of cones tangent to only one or two atoms, the rest of the elements are padded with zeros. In the case of an unsuccessful calculation, the return code `stat` will be non-zero and an error message is stored in `errmsg`.

A minimal [FORD](https://github.com/Fortran-FOSS-Programmers/ford) documentation can be built with `ford pages.md`

### C API

The C API exposes the subroutine `cone_angle_c` with the C name `cone_angle`. It's signature is the same as for the Fortran subroutine, but requires the specification of the number of atoms, `n_atoms`. 

> ⚠️ All atoms are zero-index in the C API.

```fortran
subroutine cone_angle_c(n_atoms, coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat) bind(c, name="cone_angle")
  !! Calculate cone angle, cone axis and tangent atoms
  !> Number of atoms
  integer(c_int), value, intent(in) :: n_atoms
  !> Coordinates (Å)
  real(c_double), intent(in) :: coordinates(3, n_atoms)
  !> vdW radii (Å)
  real(c_double), intent(in) :: radii(n_atoms)
  !> Index of metal atom
  integer(c_int), value, intent(in) :: index_metal
  !> Cone angle (degrees)
  real(c_double), intent(out) :: alpha
  !> Cone axis (Å)
  real(c_double), intent(out) :: axis(3)
  !> Indices of atoms tangent to cone
  integer(c_int), intent(out) :: tangent_atoms(3)
  !> Return code
  integer(c_int), intent(out) :: stat
  !> Error message
  character(c_char), intent(out) :: errmsg(*)
  
  call cone_angle(coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat)
end subroutine cone_angle_c
```

The C header file can be found [here](include/cone_angle.h). An [example](libconeangle/lib.py) is given for loading the shared library with ctypes in the Python API.

## Acknowledgements

Any published work derived from the use of libconeangle should cite the original publication for exact ligand cone angles.[^1]

- Cyrille Lavigne (@clavigne) for many discussions and for introducing me to modern Fortran
- Sebastian Ehlert (@avwk) for testing and many discussions, especially on packaging and distribution

## References

[^1]: Bilbrey, J. A.; Kazez, A. H.; Locklin, J.; Allen, W. D. Exact Ligand Cone Angles. *Journal of Computational Chemistry* **2013**, *34* (14), 1189–1197. https://doi.org/10.1002/jcc.23217.

