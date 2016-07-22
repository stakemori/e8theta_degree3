# -*- compile-command: "sage -c 'sh.eval(\"python setup.py build_ext -i\")'"-*-
from os.path import abspath, curdir, join
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

setup(
    name="E8Theta",
    ext_modules=cythonize(
        Extension("e8theta",
                  sources=["e8theta.pyx"],
                  include_dirs=[abspath(curdir)],
                  library_dirs=[join(abspath(curdir), "lib")],
                  libraries=["miyawaki_theta", "e8vectors"]),
    ),
)
