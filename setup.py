from distutils.core import setup, Extension
import os
from Cython.Build import cythonize
import numpy

PACKAGE = "spyndle"
NAME = "spyndle"
DESCRIPTION = "Toolbox for analyzing sleep spindles."
AUTHOR = "Christian O'Reilly"
AUTHOR_EMAIL = "christian.oreilly@umontreal.ca"
VERSION = "0.4.0" #__import__(PACKAGE).__version__
URL = "https://bitbucket.org/christian_oreilly/spyndle/"
DOWNLOAD_URL = URL + "downloads/" + NAME + "-" + VERSION + ".zip"
extensions = ["spyndle/io/*.pyx", "spyndle/miscellaneous/*.pyx"]

def is_package(path):
    return (
        os.path.isdir(path) and
        os.path.isfile(os.path.join(path, '__init__.py'))
        )

def find_packages(path, base="" ):
    """ Find all packages in path """
    packages = {}
    for item in os.listdir(path):
        dir = os.path.join(path, item)
        if is_package( dir ):
            if base:
                module_name = "%(base)s.%(item)s" % vars()
            else:
                module_name = item
            packages[module_name] = dir
            packages.update(find_packages(dir, module_name))
    return packages

packages=find_packages(".")


setup(
    name=NAME,
    packages=packages.keys(),
    package_dir=packages,
    package_data={'spyndle.io.tests': ['*.SIG', '*.STS', '*.bdf'],
                  'spyndle.EEG': ['*.png', '*.svg'],
                  'spyndle.DevuystDB': ['*.txt', '*.edf']},
    version=VERSION,
    description=DESCRIPTION,
    long_description=open("README.md").read(),
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=AUTHOR,
    maintainer_email=AUTHOR_EMAIL,  
    license='LICENSE.txt',
    url=URL,
    download_url=DOWNLOAD_URL,
    requires=['pandas (>= 0.11.0)', 'lxml', 'comtypes', 'sqlalchemy'],
    classifiers=["Development Status :: 3 - Alpha",
			"Environment :: MacOS X",
			"Environment :: Win32 (MS Windows)",
			"Environment :: X11 Applications",
			"Intended Audience :: Developers",
			"Intended Audience :: Science/Research",
			"License :: Free for non-commercial use",
			"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
			"Natural Language :: English",
			"Programming Language :: Python :: 3.4",
			"Topic :: Scientific/Engineering"])
	
	
	
	
	
