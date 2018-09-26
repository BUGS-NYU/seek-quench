from distutils.core import setup
from Cython.Build import cythonize

# Setup file to make cython package

setup(name='Cython Version',
      ext_modules=cythonize("src.pyx"))
