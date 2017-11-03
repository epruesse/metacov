#!/usr/bin/env python

from Cython.Build import cythonize

from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension


class BuildExt(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)

        import numpy
        self.include_dirs.append(numpy.get_include())

        import pysam
        self.include_dirs.extend(pysam.get_include())
        self.libraries.extend(pysam.get_libraries())


extensions = cythonize([
    Extension("metacov.pyfq", ["metacov/pyfq.pyx"]),
    Extension("metacov.scan", ["metacov/scan.pyx"])
])


setup(
    name="MetaCov",
    use_scm_version=True,
    packages=find_packages(),
    ext_modules=extensions,
    cmdclass={'build_ext': BuildExt},
    setup_requires=[
        'setuptools_scm',
        'cython',
        'pysam',
        'numpy'
    ],
    tests_require=[
        'pytest'
    ],
    install_requires=[
        'Click',
        'pysam',
        'numpy'
    ],
    entry_points='''
        [console_scripts]
        metacov=metacov.cli:main
    '''
)
