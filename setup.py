from distutils.core import setup

PACKAGE = "spyndle"
NAME = "spyndle"
DESCRIPTION = "Toolbox for analyzing sleep spindles."
AUTHOR = "Christian O'Reilly"
AUTHOR_EMAIL = "christian.oreilly@umontreal.ca"
URL = "http://sourceforge.net/p/spyndle/"
VERSION = __import__(PACKAGE).__version__

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=open("README.txt").read(),
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    license='LICENSE.txt',
    download_url='https://github.com/USER/PROJECT/tarball/master',
    packages=['spyndle'],
    install_requires=[],
)