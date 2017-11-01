#!/usr/bin/env python3

from setuptools import setup, find_packages
from Cython.Build import cythonize
from distutils.extension import Extension

def make_ext(modname, pyxfilenames):
    import pysam
    import numpy
    return Extension(
        name=modname,
        sources=pyxfilenames,
        extra_link_args=pysam.get_libraries(),
        include_dirs=pysam.get_include() + [numpy.get_include()],
        define_macros=pysam.get_defines()
    )


extensions = cythonize([
    make_ext("metacov.pyfq", ["metacov/pyfq.pyx"]),
    make_ext("metacov.scan", ["metacov/scan.pyx"])
])


setup(
    name="MetaCov",
    use_scm_version=True,
    packages=find_packages(),
    ext_modules=extensions,
    setup_requires=[
        'setuptools_scm',
        'cython'
    ],
    tests_require=[
        'pytest'
    ],
    install_requires=[
        'Click',
        'pysam'
    ],
    entry_points='''
        [console_scripts]
        metacov=metacov.cli:main
    '''
)
