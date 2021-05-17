# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 17:11:16 2019

@author: mlanger
"""
import pyphase.dataset as Dataset
import pyphase.phaseretrieval as Phaseretrieval
import matplotlib.pyplot as pyplot
from pyphase.config import *

# data = Dataset.ESRF('/data/id19/bones01/bones1/max/Holocalibration/', 'phantom_0.7um_20.5')
data = Dataset.ESRF('/data/id19/bones01/bones3/max/holodata/knee/lbtoKneeWTOA8weekFeb15/', 'lbto_18_OA_8weeks')

retriever=Phaseretrieval.Mixed(data)
retriever.Lcurve(data, 0)
retriever.DisplayLcurve(data)
