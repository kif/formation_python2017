from setuptools import setup, find_packages
from distutils.extension import Extension
from Cython.Build import cythonize



setup(
    name           = "nbody",
    packages       = find_packages(exclude=['tests*']),
#    ext_modules=cythonize([Extension("*", ["*/*/*.pyx"]), Extension("*", ["*/*.pyx"])]),

)
