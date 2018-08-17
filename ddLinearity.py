#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 12:11:00 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *

def linearity(data, threshold=66000.):
  """
  Returns a linearised ramp of the CryoNIRSP H2RG. Assuming only a quadratic 
  term is removed
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be linearity corrected.
    threshold : float
        threshold for the quadratic fit. Everything above will be excluded.

    Returns
    -------
    linearityCorrected : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that is linearity corrected.

    Raises
    ------

    Notes
    -----
    The linearized signal will be y' = y - c*t**2
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> t = np.arange(12,dtype=np.float32)
    >>> c = np.zeros((3,10,10),dtype=np.float32)
    >>> c[0,:,:] = .3
    >>> c[1,:,:] = 2
    >>> c[2,:,:] = 1
    >>> val = np.polyval([.3,2,1],t)
    >>> lin = np.polyval([2,1],t)
    >>> data = np.tile(val[:,None,None],(1,10,10))
    >>> res = linearity(data,c,t)
   """
  if len(data.shape) == 3:
    # only a single ramp is provided but matrix quadfit expects another axis
    data = np.expand_dims(data, axis=0)
    
  # in this case we use the NDR number as our 'time' axis
  dataTime = np.arange(data.shape[1])
  # fit quadratic to data
  a,b,c = matrixQuadfit(data,threshold=90,ignoreRef=False)
  
  # Use numpy's multiply capability to to the multiplication along the right axis
#  nonLinear = np.multiply(np.expand_dims(linearityCoefficients[0,:,:],axis=0),
#                          dataTime[:,None,None]**2.)
  nonLinear = np.multiply(c,dataTime[:,None,None]**2.)
  
  linearityCorrected = np.squeeze(data) - nonLinear
  
  return linearityCorrected

t = np.arange(12,dtype=np.float32)
c = np.zeros((3,10,10),dtype=np.float32)
c[0,:,:] = .3
c[1,:,:] = 2
c[2,:,:] = 1
val = np.polyval([.3,2,1],t)
val[2] = 100
lin = np.polyval([2,1],t)
data = np.tile(val[:,None,None],(1,10,10))
res = linearity(data)
#
fig, ax=plt.subplots()
ax.plot(t,val,t,lin, 'r', t, res[:,0,0],'g')
plt.show()