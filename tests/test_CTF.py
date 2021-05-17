# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:55:38 2019

@author: mlanger
"""

import pyphase.dataset as Dataset
import pyphase.phaseretrieval as Phaseretrieval
import matplotlib.pyplot as pyplot
from pyphase.config import *

# TEST CTF


# TEST real ID19 data
phantom = Dataset.ESRF('/data/id19/bones01/bones1/max/Holocalibration/', 'phantom_0.7um_20.5')
retriever=Phaseretrieval.CTF(phantom)

# TEST: Default constructor

# TEST: Get Alpha
# TEST: alpha single value
#Alpha = retriever.Alpha
#fig = pyplot.imshow(Alpha)
#pyplot.colorbar()

# TEST: Reconstruct
#retriever.ReconstructProjection(phantom, 0)

# TEST: alpha two values, a < b
retriever.alpha = [1e-8, 1e-3]

# TEST: get Alpha
#Alpha = retriever.Alpha
#fig = pyplot.imshow(Alpha)
#pyplot.colorbar()

# TEST: Reconstruct
retriever.ReconstructProjection(phantom, 0)
viewer.DisplayImage(phantom.GetImage(0))
# TEST: simulated ID19 data 
#phantom = Dataset.ESRFID19('/data/id19/bones01/bones1/max/HoloSim/max/HoloSim', 'holosim')
#phantom.preprocess=0
#phantom.correct_shifts=0
#phantom.correct_motion=0

#retriever.alpha=1e-99
#retriever.ReconstructProjection(phantom, 0)
