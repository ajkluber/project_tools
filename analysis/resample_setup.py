''' Build script for cython resample extension

Run following on command line to build the cython extension used
in the bootstrapping script.

python resample_setup.py build_ext --inplace

The result of this command is a dynamic library (binary) resample_histo.so 
that can be imported into python:

>>> import resample_histo

'''
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("resample_histo.pyx")
)
