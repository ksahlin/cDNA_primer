from setuptools import setup, Extension
import numpy

ext_modules = [
        #Extension("BioReaders", ["BioReaders.pyx"]),
        Extension("c_branch", ["c_branch.pyx"],
            include_dirs=[numpy.get_include()]),
        Extension("modified_bx_intervals.intersection_unique", ["modified_bx_intervals/intersection_unique.pyx"]),
]
setup(
		name = 'c_cDNA_primer',
        #cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules,
        author_email='etseng@pacificbiosciences.com',
        author='Elizabeth Tseng'
)
