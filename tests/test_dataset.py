# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 18:13:18 2019

@author: mlanger
"""


import pyphase.dataset as Dataset
import numpy as np
import matplotlib.pyplot as plt
import pyphase.utilities as Utilities
from pyphase.config import *


# Open ID19 data set

data = Dataset.ESRF('phantom_0.7um_20.5', '/gpfs/offline1/staff/tomograms/users/maxlan/phantom07/')

# Align dataset
#data.AlignProjection(0)
data.Align()

# Display Shfits

#data.DisplayShifts(0)

# GetShift
#tmp=data.GetShift(0,3)
#print(tmp)
# Motion correction

#data.CalculateMotion()

# Initialise

#data.Initialise()


# Get projection

#im = data.GetImage(0)
#im = data.GetProjection(0,2)
#viewer.DisplayImage(im)

#TODO: Perhaps nicer would be phantom.GetProjection(distance=4, projection=0)

# Populate
#phantom.Populate()

# Initialise
#data.Initialise()
#

#data.Preprocess()

# Rotation axis
#data.CalculateAxisPosition(2) 

# GetSinogram
#sinogram=data.GetSinogram()
pass
