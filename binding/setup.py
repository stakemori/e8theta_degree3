# -*- compile-command: "sage -c 'sh.eval(\"python setup.py build_ext -i\")'"-*-
import os
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

setup(
    name="E8Theta",
    ext_modules=cythonize(
        Extension("e8theta",
                  sources=["e8theta.pyx"],
                  include_dirs=[os.path.abspath(os.path.curdir)],
                  library_dirs=[os.path.abspath(os.path.curdir)],
                  libraries=["miyawaki_theta", "e8vectors"]),
    ),
)
