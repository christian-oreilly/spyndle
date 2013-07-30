from setuptools import setup
import os

PACKAGE = "spyndle"
NAME = "spyndle"
DESCRIPTION = "Toolbox for analyzing sleep spindles."
AUTHOR = "Christian O'Reilly"
AUTHOR_EMAIL = "christian.oreilly@umontreal.ca"
VERSION = __import__(PACKAGE).__version__
URL = "https://bitbucket.org/christian_oreilly/spyndle/"
#DOWNLOAD_URL = URL + "downloads/" + NAME + "-" + VERSION + ".zip"


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
    package_dir=packages,
    version=VERSION,
    description=DESCRIPTION,
    long_description=open("README.txt").read(),
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license='LICENSE.txt',
    packages=packages.keys(),
    url=URL,
    setup_requires=['bitbucket-distutils >= 0.1.2', 'pandas >= 0.11.0', 'lxml', 'comtypes'])

