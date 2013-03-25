from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import os
import numpy

ext_modules = []
ext_modules.append(Extension(
    "chebyshev",
    ["chebyshev.pyx"],
    libraries=["m"],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
))

ext_modules.append(Extension(
    "_spline",
    ["spline.pyx"],
    libraries=["m"],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
))

ext_modules.append(Extension(
    "_pchip",
    ["_pchip.pyx"],
    libraries=["m"],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
))

ext_modules.append(Extension(
    "_rbf",
    ["_rbf.pyx"],
    libraries=["m"],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
))

setup(
    include_dirs=[numpy.get_include()],
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules, 
)
