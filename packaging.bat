pushd C:\Users\revestech\Dropbox\code\Spyndle

# register the package with PyPI, creates a source and
# an egg distribution. Do not upload it to PyPI because of 
# licence incompatibility. Using bitbucket-distutils, the toolbox 
# uploaded directly to the bitbucket source repository. 
python setup.py release
