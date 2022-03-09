# from distutils.core import setup
# from distutils.extension import Extension
# from Cython.Distutils import build_ext
#
# module =[ Extension('panel_data', sources=['panel_data.pyx']),
#           Extension('dpd', sources=['dynamic_panel_data_model.pyx'])
#           ]
#
#


# setup(
#     name='cythonTest',
#     version='1.0',
#     author='jetbrains',
#     ext_modules=module
# )


from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(

    ext_modules = cythonize("*.pyx"),

    include_dirs=[numpy.get_include()]
)