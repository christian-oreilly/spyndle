# __init__.py
__version__ = "0.1.1a.dev1"

#from propagation.XCST import computeXCST
from filters import Filter
from sleepCycles import cycleDefinitions, computeDreamCycles
from STransform import computeMST
from fastST import computeFastST, computeFastST_real