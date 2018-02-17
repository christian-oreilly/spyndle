# __init__.py
__version__ = "0.3.0"


from .miscellaneous.filters import Filter, channelType
from .miscellaneous.sleepCycles import CycleDefinitions, computeDreamCycles
from .miscellaneous.STransform import computeMST, computeST
from .miscellaneous.fastST import computeFastST, computeFastST_real

try:
    from .miscellaneous.line import Line, Point
except ImportError:
    pass
    #class Line
    #class Point

from .miscellaneous.utils import setUnbufferedPrint, Unbuffered, terminate_process, diff2

# For backward compatibility
from .miscellaneous.sleepCycles import cycleDefinitions

#import io, DevuystDB, detector, miscellaneous, propagation, EEG
