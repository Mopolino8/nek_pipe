from distutils.core import setup, Extension

module1 = Extension('revolve', 
                     sources =
                     ['revolve_python.c','revolve.cpp','wrap_revolve.cpp'],
                     include_dirs=['./'],
                     library_dirs=['./'])

setup (name = 'Revolve',
               version = '0.1',
               description = 'This module implements the Revolve algorithm of Algorithmic Differentiation',
               ext_modules = [module1])
