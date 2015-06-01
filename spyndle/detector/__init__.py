# __init__.py


"""
    Initialization of the spyndle.detector module.
    
    Copyright (C) 2013  Christian O'Reilly

    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author.
    
    If this code is used for research purpose, the reference [1] should be
    cited in the derived publication.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : November 1, 2013

"""

from spyndle.detector.transientDetectors import TransientDetector, ThresholdDetector, DetectedEvent

from spyndle.detector.kcomplexDetectors import TeagerKComplexDetector

from spyndle.detector.spindleDetectors import SpindleDetectorRMS, SpindleDetectorSigma, \
                                              DetectedSpindle, SpindleDetectorAmp, \
                                              SpindleDetectorTeager, SpindleDetectorRSP, SpindleDetectorMA

from spyndle.detector.detectorEvaluation import DetectorEvaluator

from spyndle.detector.slowWaveDetectors import MassiminiSlowWaveDetector

from spyndle.detector.muscularArtifactDetectors import MuscularArtifactDetector
from spyndle.detector.alphaArtifactDetectors import AlphaArtifactDetector
 