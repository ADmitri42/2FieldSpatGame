from distutils.core import setup
from distutils.extension import Extension

import numpy as np
from Cython.Distutils import build_ext
# from Cython.Build import cythonize

# OPT="-DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes" python3 setup.py build_ext
# export CFLAGS='-I/home/admitri/.local/lib/python3.6/site-packages/numpy/core/include'

# setup(ext_modules = cythonize(
#            ["meangame.pyx", "evolve.cc"],                 # our Cython source
#            include_path = [np.get_include(),],
#            language="c++",             # generate C++ code
#       ))

setup(
    ext_modules=[
        Extension(
            "spatgames",
            [
                "./games/meangame.pyx",
                "./games/cpp/games.cpp",
                "./games/cpp/utilities.cpp",
                "./games/cpp/spatgame.cpp",
                "./games/cpp/doublefield.cpp",
            ],
            language="c++",
            include_path=[
                np.get_include(),
            ],
            extra_compile_args=["-std=c++17"],
        )
    ],
    cmdclass={"build_ext": build_ext},
)
