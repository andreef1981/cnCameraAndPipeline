#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 12:11:00 2018

@author: andreef
"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *

def linearity(data,
         linearityCoefficients,
         dataTime):
  """
  Returns a linearised ramp of the CryoNIRSP H2RG. Assuming only a quadratic 
  term is removed
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be linearity corrected.
    linearityCoefficients : (3, 2048, 2048) ndarray, float32
        3D data cube that contains the parameters to a quadratic curve 
        y = a*t**2 + b*t + c . Where the [0,:,:] contains the quadratic term.
        Units of [b]= ADU/second, [a]=ADU/second**2
    dataTime : (#NDRs,) ndarray, float32
        vector containing the timestamps for the NDRs up the ramp.
        Units [dataTime] = seconds    

    Returns
    -------
    linearityCorrected : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that is linearity corrected.

    Raises
    ------
    AssertationError
        If the shapes of the input arrays do not match.

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
  
  # Use numpy's multiply capability to to the multiplication along the right axis
  nonLinear = np.multiply(np.expand_dims(linearityCoefficients[0,:,:],axis=0),
                          dataTime[:,None,None]**2.)
  
  linearityCorrected = data - nonLinear
  
  return linearityCorrected

#t = np.arange(12,dtype=np.float32)
#c = np.zeros((3,10,10),dtype=np.float32)
#c[0,:,:] = .3
#c[1,:,:] = 2
#c[2,:,:] = 1
#val = np.polyval([.3,2,1],t)
#lin = np.polyval([2,1],t)
#data = np.tile(val[:,None,None],(1,10,10))
#res = linearity(data,c,t)
#
#fig, ax=plt.subplots()
#ax.plot(t,val,t,lin, 'r', t, res[:,0,0],'g')
#plt.show()