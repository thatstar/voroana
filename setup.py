from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize


extensions = [
    Extension(
        name="voroana.voronoi",
        sources=[
            "voroana/voronoi.pyx",
            "voroana/voro++-0.4.6/src/voro++.cc"
        ],
        include_dirs=[
            "voroana/voro++-0.4.6/src"
        ],
        language="c++"
    )
]


setup(
    name='voroana',
    version='0.1',
    packages=['voroana'],
    ext_modules=cythonize(extensions),
    url='',
    license='GPL v3',
    author='surui',
    author_email='surui1986@gmail.com',
    description='A flexiable Voronoi analysis code for disordered materials.'
)
