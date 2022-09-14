from setuptools import setup
from Cython.Build import cythonize
import numpy
import sys

sys.setrecursionlimit(50000)

setup(
    name='RMF RBM Library',
    ext_modules=cythonize("rmf_rbm_hybrid.pyx"),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)