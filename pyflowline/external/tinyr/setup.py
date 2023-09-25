
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


src_path = 'tinyr'

extensions = [Extension(
                'tinyr',
                [os.path.join(src_path, 'tinyr.pyx')],
                extra_compile_args = ['-O3', '-Wall'],
                extra_link_args = ['-g'],
                )]

setup(
    name='tinyr',
    version='0.1',
    description='Fast 2D R-Tree implementation in cython',
    author='Matthias Simon',
    url='',
    packages=['tinyr'],
    ext_modules=extensions,
    cmdclass=dict(build_ext=build_ext),
)

