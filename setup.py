import numpy

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('cmesh',
              sources=['cmesh.pyx', 'mesh.c'],
              include_dirs = [numpy.get_include()],
              extra_compile_args = ['-Wno-unused-function',
                                    '-Wno-#warnings'])
    ]

setup(
    name = 'cmesh',
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext},
    )
