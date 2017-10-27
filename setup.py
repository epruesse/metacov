#!/usr/bin/env python3

from setuptools import setup, find_packages
from Cython.Build import cythonize
from distutils.extension import Extension


extensions = cythonize([
    Extension("metacov.pyfq", ["metacov/pyfq.pyx"])
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
