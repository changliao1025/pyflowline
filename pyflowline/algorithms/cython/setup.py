from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize
import numpy as np
extensions = [Extension("kernel", ["kernel.pyx"], language="c++")]

setup(
    name='kernel',
    ext_modules=cythonize(extensions),
    zip_safe=False,
)