"""
setup.py for package "entropy"
(using Cython to interface the C library) 
2014-09-30 first version wwith setuptools
2019-10-14 now using python3
2022-11-25 now using setuptools instead of distutils
2023-12-04 now the install is even more generic on Linux systems
2024-10-11 simplfied version
"""
from setuptools import Extension, setup
from Cython.Build import cythonize

# to have an annotated html:
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = False
DO_ANNOTATE = False

# Third-party modules - we depend on numpy for everything
import numpy
import sys
import os

# Obtain the numpy include directory:
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()
print("using NumPy from : ", numpy_include)
#cwd = os.getcwd(); print("current directory : ", cwd)
#local_dir_raw = os.path.realpath(__file__); print("source directory : ", local_dir_raw)
#local_dir = os.path.dirname(local_dir_raw); print("source dir extrc : ", local_dir)
#print(os.environ)
#local_dir = os.environ['BUILD_DIR'] # not passed from Makefile!
local_dir = os.environ['PWD'] # 2023-12-04: a trick...
print("source origin    : ", local_dir)

import subprocess 
#LIBS_PYTHON = subprocess.run(["python3-config","--libs"], capture_output=True).stdout # python3.7
#LIBS_PYTHON = subprocess.run(["python3-config","--libs"], stdout=subprocess.PIPE).stdout # python3.5
#LIBS_PYTHON = LIBS_PYTHON.rstrip().decode().split()
#LDFLAGS_PYTHON = subprocess.run(["python3-config", "--ldflags"], capture_output=True).stdout # python3.7
LDFLAGS_PYTHON = subprocess.run(["python3-config", "--ldflags"], stdout=subprocess.PIPE).stdout # python3.5
LDFLAGS_PYTHON = LDFLAGS_PYTHON.rstrip().decode().split()
#print("LIBS_PYTHON", LIBS_PYTHON)
print("LDFLAGS_PYTHON", LDFLAGS_PYTHON + ["-Wl,-undefined,error"]) 


print(sys.version)
import platform
if (platform.system()=='Darwin'):
		LIBS = ['c++']      # 2018-11-27: stdc++ replaced by c++ (for macos clang at least)
		LDFLAGS_PYTHON = ['-mmacosx-version-min=12'] # 2021-02-15: replaced 10.12 by 10.15 
		                                             # 2022-11-26: replaced 10.15 by 12 (Monterey)
		CFLAGS = ['-Wno-parentheses-equality']       # 2023-10-06: to suppress useless warnings from cython code
else:
		LIBS = ['stdc++', 'm']
		CFLAGS = []
		
print(LIBS)

INC_DIR  = ['.', '../include', '../../include', '/opt/local/include', numpy_include]
LIBS_DIR = ['../../lib', '/usr/lib', '/usr/local/lib', '/opt/local/lib']
LIBS     = ['nn_stats', 'ANN', 'pthread'] + LIBS #'m', # 2020-07-22: moved LIBS at the end (for stdc++ on linux)

if (platform.system()=='Linux'):
        INC_DIR.append(local_dir+'/../../include')
        LIBS_DIR.append(local_dir+'/../../lib')

nn_module = Extension(
                name = "nn_stats",
				sources = ['nn_stats/nn_statistics.pyx'], 
                include_dirs = INC_DIR,
                extra_compile_args = CFLAGS,
                libraries    = LIBS,
                library_dirs = LIBS_DIR,
                extra_link_args=LDFLAGS_PYTHON,
                define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],    
#               cython_directives = {"embedsignature": True}, 
                )
          
setup(name = 'nn_stats',
      version = '1.0', 
#      date ='2025.01.16',
      description = "nearest neighbors local estimates",
      author      = "Nicolas B. Garnier",
      author_email= "nicolas.garnier@ens-lyon.fr",
#      zip_safe=False,  # 2022-11-25, to work with cimport for pxd files when using them from a dependent package
                        # 2023-02-14 : line above commented
#	  cmdclass = {'build_ext': build_ext},
#	  ext_package = 'nn_stats',
	  ext_modules = cythonize([nn_module], annotate=DO_ANNOTATE)#, compiler_directives={'embedsignature': True}), 
	)

