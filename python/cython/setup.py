from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

setup(
	ext_modules = cythonize(Extension(name='cgal_partition',
				sources=['cgal_partition.pyx', 'cgal_partition_2_wrapper.cpp'],
				include_dirs=[np.get_include(), '/usr/include'],#really didn't like that '/usr/include/CGAL'],
				libraries=['gmp', 'mpfr'],
				library_dirs=['/usr/lib/x86_64/-linux-gnu/'],
				language='c++'))
)



