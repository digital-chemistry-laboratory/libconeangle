---
title: 'libconeangle: A library for calculating exact ligand cone angles'
tags:
  - Fortran
  - Python
  - ligands
  - molecular features
  - quantum chemistry
authors:
  - name: Kjell Jorner 
    orcid: 0000-0002-4191-6790
    affiliation: "1, 2, 3"
  - name: Alán Aspuru-Guzik^[Corresponding author]
    orcid: 0000-0002-8277-4434
    affiliation: "1, 2, 4, 5, 6, 7"
affiliations:
  - name: Department of Computer Science, University of Toronto, 40 St. George St, Toronto, Ontario M5S 2E4, Canada
    index: 1
  - name: Department of Chemistry, Chemical Physics Theory Group, 80 St. George St., University of Toronto, Ontario M5S 3H6, Canada
    index: 2
  - name: Department of Chemistry and Chemical Engineering, Chalmers University of Technology, Kemigården 4, SE-41258, Gothenburg, Sweden
    index: 3
  - name: Department of Chemical Engineering & Applied Chemistry, 200 College St., University of Toronto, Ontario M5S 3E5, Canada
    index: 4
  - name: Department of Materials Science & Engineering, 184 College St., University of Toronto, Ontario M5S 3E4, Canada
    index: 5
  - name: Vector Institute for Artificial Intelligence, 661 University Ave. Suite 710, Toronto, Ontario M5G 1M1, Canada
    index: 6
  - name: Lebovic Fellow, Canadian Institute for Advanced Research (CIFAR), 661 University Ave., Toronto, Ontario M5G 1M1, Canada
    index: 7
date: 30 March 2022
bibliography: paper.bib
---
# Summary

Applications of machine learning in chemistry require the representation of
molecules in a machine-readable format. Traditionally, this is done via
so-called expert features that capture target-relevant information, while modern
deep learning approaches allow implicit featurization via, *e.g.*, graph neural
networks [@duvenaud_2015; yang_2019]. Expert featurization of ligands for
homogenous transition metal catalysis has a long history [@durand_2019], with
the introduction of the ligand cone angle and electronic parameter by Tolman and
co-workers in the 1970s as an outstanding example [@tolman_1977]. Today there is
a range of expert features developed for different ligand types [@durand_2019],
but the cone angle remains one of the most important [@gensch_2022].

`libconeangle` is library for calculating exact ligand cone angles based on the
recipe of Allen and co-workers [@bilbrey_2013]. It is written in Modern Fortran
[@kedward_2022] with C and Python APIs, enabling easy incorporation into
featurization packages in different languages (Python, Julia, C++). It aims to
provide a fast reference implementation for cone angle calculations.

# Statement of need

Implementations of exact ligand cone angles include the original Mathematica
version by Allen and co-workers and two Python implementations by Jorner *et al.* 
in `ᴍᴏʀғᴇᴜs`(MIT license) [@morfeus], and Wheeler and co-workers in
`AaronTools.py` (GPL license) [@aarontools; @ingman_2021]. The Mathematica
version requires proprietary software to run, while the Python versions are slow
for large ligands and restricted to the Python ecosystem. For example, the
calculation of the AdBrettPhos ligand with `XX` heavy atoms takes ca 7 seconds
on a modern laptop computer with `ᴍᴏʀғᴇᴜs`. This calculation needs to be repeated
for each of the `XXX` low-lying conformers, leading to a total computational
time of XX seconds (single-core) to fully parametrize the ligand. For generation
of libraries with thousands of ligands, the computational cost becomes
significant.

We have written `libconeangle` as a fast reference implementation of the exact
ligand cone angle algorithm. A calculation of AdBrettPhos takes ca 40
milliseconds instead of 7 seconds. `libconeangle` is written in Modern Fortran,
with a C API that allows access from any language with a C Foreign Function
Interface (CFFI). A Python API is packaged for installation via `pip` and
`conda`, and it would be straightforward to package also for other languages
such as with [BinaryBuilder](https://binarybuilder.org) for Julia. The MIT
license further allows wide adoptation.

# Features and implementation

The original recipe for calculating cone angles by Tolman is incorporated into
various software packages, for example in Aarontools.py [@aarontools; @ingman_2021].
Tolman's definition was developed with phosphine ligands in mind and suffers
from a number of problems. To adress these shortcomings, Allen and co-workers
introduced the exact ligand cone angles in 2013 [@bilbrey_2013].

`libconeangle` provides the `cone_angle` Fortran subroutine to calculate the
exact ligand cone angle. The C API allows access via the similarly named
`cone_angle` C subroutine made possible via the Fortran-C interoperability
introduced in Fortran 2003. `libconeangle` takes advantage of recent advances in
Fortran tooling [@curcic_2021], relying on functions from the Fortran Standard Library [@stdlib] and it can be built with the Fortran Package Manager [@fpm].
Python wheels are automatically built on GitHub Actions runners for each release,
taking advantage of `cibuildwheel` [@cibuildwheel] and `scikit-build` [@scikit-build].
`libconeangle` is also built on [conda-forge](https://conda-forge.org) for installation with the `conda` package manager.

> UPDATE when the conda-forge recipe is cleared.

The Python API exposes the single function `cone_angle` that can be used
directly from `libconeangle`, or preferably, through some higher-level library
that also includes functions for reading ligand structures files and the
associated van der Waals radii. In fact, `ᴍᴏʀғᴇᴜs` now calculates exact ligand
cone angles using `libconeangle` by default. Under the hood, the Python API uses
`ctypes` and `numpy.ctypeslib` [@harris_2020] to load the shared library and
handle passing of arrays back and forward between Python and Fortran. These
technical details are completely hidden to the end user. 

# Acknowledgements

Cyrille Lavigne and Sebastian Ehlert are acknowledged for discussions during the
development of the library. K.J. acknowledges funding through an International
Postdoc grant from the Swedish Research Council (No. 2020-00314). A. A.-G.
thanks Dr. Anders G. Frøseth for his generous support. A. A.-G. also
acknowledges the generous support of Natural Resources Canada and the Canada 150
Research Chairs program.

> TODO: Add Alán additional acknowledgments

# References