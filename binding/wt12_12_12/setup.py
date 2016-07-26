from os.path import abspath, curdir, join, pardir
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

setup(
    name="THETA12_12_12_CYTHON",
    ext_modules=cythonize(
        Extension("theta12_12_12_cython",
                  sources=["theta12_12_12_cython.pyx"],
                  include_dirs=[abspath(curdir), abspath(pardir)],
                  library_dirs=[join(abspath(pardir), "lib")],
                  libraries=["theta12_12_12_c", "e8vectors"]),
    ),
)
