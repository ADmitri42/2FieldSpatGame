import sys
from distutils.core import setup
from distutils.extension import Extension

import numpy as np
from Cython.Distutils import build_ext


if sys.platform == "darwin":
    comp_args = {
        "extra_compile_args": ["-std=c++17"]
    }
else:
    comp_args = {
        "extra_compile_args": ["-std=c++17", '-fopenmp'],
        "extra_link_args": ['-fopenmp']
    }

setup(
    ext_modules=[
        Extension(
            "spatgames",
            [
                "./games/meangame.pyx",
                "./games/cpp/spatgame.cpp",
                "./games/cpp/triangular_field_game.cpp",
            ],
            language="c++",
            include_path=[
                np.get_include(),
            ],
            **comp_args
        )
    ],
    cmdclass={"build_ext": build_ext},
)
