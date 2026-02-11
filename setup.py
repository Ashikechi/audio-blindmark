from setuptools import setup, Extension
from Cython.Build import cythonize
from numpy import get_include


cython_extensions: list[Extension] = cythonize([
        Extension(
            '*',
            ['src/**/*.pyx'],
            include_dirs = [get_include()],
            extra_compile_args = ['-O3'],
            language = 'c++',
        ),
    ],
    compiler_directives = {
        'embedsignature': True,
        'language_level': '3str',
        'warn.maybe_uninitialized': True,
        'warn.unused': True,
        'warn.unused_arg': True,
        'warn.unused_result': True
    }
)

setup(
    ext_modules = cython_extensions,
)
