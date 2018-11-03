#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 12:26:24 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
import matplotlib.pyplot as plt
#from cnPipeline import *
from helperFunctions import *
#from datetime import datetime
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from numpy.linalg import linalg

cw = 3934.3
desiredDispersion = 0.0174347 #nm/pixel
a,b = cnRefSpectrum(cw,"175um",desiredDispersion)

fig, ax=plt.subplots()
ax.plot(a,b)