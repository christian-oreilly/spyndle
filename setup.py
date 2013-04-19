from setuptools import setup

PACKAGE = "spyndle"
NAME = "spyndle"
DESCRIPTION = "Toolbox for analyzing sleep spindles."
AUTHOR = "Christian O'Reilly"
AUTHOR_EMAIL = "christian.oreilly@umontreal.ca"
VERSION = __import__(PACKAGE).__version__
URL = "https://bitbucket.org/christian_oreilly/spyndle/"
#DOWNLOAD_URL = URL + "downloads/" + NAME + "-" + VERSION + ".zip"

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=open("README.txt").read(),
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license='LICENSE.txt',
    packages=['spyndle'],
    url=URL,
    setup_requires=['bitbucket-distutils >= 0.1.2'])

#, install_requires=[], download_url=DOWNLOAD_URL,
