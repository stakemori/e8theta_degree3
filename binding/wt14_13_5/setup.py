from os.path import abspath, curdir, join, pardir
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

setup(
    name="THETA14_13_5_CYTHON",
    ext_modules=cythonize(
        Extension("theta14_13_5_cython",
                  sources=["theta14_13_5_cython.pyx"],
                  include_dirs=[abspath(curdir), abspath(pardir)],
                  library_dirs=[join(abspath(pardir), "lib")],
                  libraries=["theta14_13_5_c", "e8vectors"]),
    ),
)
