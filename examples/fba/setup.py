from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy



extensions = [Extension('*', ['*.pyx'],
                        include_dirs=['/Users/pstjohn/Research/dFBA/sim',
                                      '/Users/pstjohn/Packages/sundials_install/include',
                                      numpy.get_include()],
                        # library_dirs=['/Users/pstjohn/Packages/sundials_install/lib'],
                        # libraries=['sundials_cvode', 'sundials_nvecserial', 'glpk'],
                        # extra_link_args=["-g"],
                        extra_compile_args=['-Wno-unused-function'],
              )]


setup(
    name='GlucoseKinetics',
    ext_modules=cythonize(extensions),
    cmdclass = {'build_ext': build_ext},
    script_args = ['build_ext'],
    options = {'build_ext':{'inplace':True, 'force':True}}
)



