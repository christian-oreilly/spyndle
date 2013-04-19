pushd C:\Users\revestech\Dropbox\Spyndle

# register the package with PyPI, creates a source and
# an egg distribution. Do not upload it to PyPI because of 
# licence incompatibility. Using bitbucket-distutils, the toolbox 
# uploaded directly to the bitbucket source repository. 
#python setup.py sdist bdist upload -R christian_oreilly/spyndle  register
python setup.py release
