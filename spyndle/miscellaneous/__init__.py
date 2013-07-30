# __init__.py
__version__ = "0.1.1a.dev1"

#from propagation.XCST import computeXCST
from filters import Filter, channelType
from sleepCycles import cycleDefinitions, computeDreamCycles
from STransform import computeMST, computeST
from fastST import computeFastST, computeFastST_real