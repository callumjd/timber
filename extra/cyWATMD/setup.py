from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

setup(
	ext_modules = cythonize("example.pyx"),
	include_dirs=[numpy.get_include()]
)

