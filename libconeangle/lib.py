"""Interaction with shared library using ctypes."""
import ctypes
from ctypes import c_char_p, c_double, c_int, POINTER
import os
from pathlib import Path
import sys
from typing import Optional

import numpy as np
from numpy.ctypeslib import ndpointer

# Load shared library
path = Path(__file__).parent
name = "libconeangle"
winmode: Optional[int] = None
if sys.platform == "win32":
    name = f"{name}.dll"
    os.add_dll_directory(path)
    winmode = 0
elif sys.platform.startswith("linux"):
    name = f"{name}.so"
elif sys.platform == "darwin":
    name = f"{name}.dylib"
else:
    raise ImportError(f"Your OS is not supported: {sys.platform}")

lib = ctypes.CDLL(str(path / name), winmode=winmode)

# Set argument types of subroutine
lib.cone_angle.argtypes = [
    c_int,
    ndpointer(dtype=np.float64, ndim=2, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    c_int,
    POINTER(c_double),
    ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),
    POINTER(c_int),
    c_char_p,
]
