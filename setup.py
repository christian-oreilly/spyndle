from distutils.core import setup

PACKAGE = "spyndle"
NAME = "spyndle"
DESCRIPTION = "Toolbox for analyzing sleep spindles."
AUTHOR = "Christian O'Reilly"
AUTHOR_EMAIL = "christian.oreilly@umontreal.ca"
VERSION = __import__(PACKAGE).__version__
URL = "https://bitbucket.org/christian_oreilly/spyndle/"
DOWNLOAD_URL = URL + "downloads/" + NAME + "-" + VERSION + ".zip"

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=open("README.txt").read(),
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    license='LICENSE.txt',
    download_url=DOWNLOAD_URL,
    packages=['spyndle'],
    install_requires=[],
)