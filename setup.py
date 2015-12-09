import os
from distutils.core import setup
# from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

def scandir(dir, files=[]):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            files.append(path.replace(os.path.sep, ".")[:-4])
        elif os.path.isdir(path):
            scandir(path, files)
    return files


def makeExtension(extName):
    extPath = extName.replace(".", os.path.sep)+".pyx"
    return Extension(
        extName,
        [extPath],
        include_dirs=['.',
                      '/Users/pstjohn/Packages/sundials_install/include',
                      numpy.get_include()],
        library_dirs=['/Users/pstjohn/Packages/sundials_install/lib'],
        libraries=['sundials_cvode', 'sundials_nvecserial', 'glpk'],
        )


extNames = scandir('sim')

extensions = [makeExtension(name) for name in extNames]

setup(
  name="dFBA",
  ext_modules=extensions,
  # cmdclass = {'build_ext': build_ext},
  # script_args = ['build_ext'],
  # options = {'build_ext':{'inplace':True, 'force':True}}
)
