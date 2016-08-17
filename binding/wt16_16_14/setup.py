from os.path import abspath, curdir, join, pardir
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

setup(
    name="THETA16_16_14_CYTHON",
    ext_modules=cythonize(
        Extension("theta16_16_14_cython",
                  sources=["theta16_16_14_cython.pyx"],
                  include_dirs=[abspath(curdir), abspath(pardir)],
                  library_dirs=[join(abspath(pardir), "lib")],
                  libraries=["theta16_16_14_c", "e8vectors", "rank16_vectors"]),
    ),
)
