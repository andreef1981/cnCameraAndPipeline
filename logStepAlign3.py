# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 08:48:11 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
Revision history
----------------
"""
import numpy as np
#import matplotlib.pyplot as plt

def logStepAlign3(centerWavelength,
                  nrSteps,
                  fraction=0.1,
                  growthFactor=0.1759):
  #???: are angles the right units to return?
  
  """
  Function to generate a logarithmic step pattern starting at (0,0)
  
   Parameters
    ----------
    centerWavelength : float32
        selected center wavelength for the scan in nm
    nrSteps : uint16
        # of steps along the spiral.
    fraction : float16
        optional input specifying how much of the spectral range on the
        detector will be scanned over
    growthFactor : optional, float, default=0.1759
        optional growth factor for the steps
    
    
    Returns
    -------
    angle : (nrSteps,1) ndarray, float32
        vector containing the grating angles for the scan
    

    Raises
    ------

    Notes
    -----
    
    Logarithmic steps inspired by Nautilus
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> 

   """
  grating_blaze = 63.9
  n = 31.6  
  theta = -5.5
  Lb = 932. 
  pixel_size = 18. 
  spectral_pixels = 1024
  
  blorder = np.around(np.sin(np.radians(grating_blaze))*2. /
                                   (centerWavelength*n*1e-6))
  beta = np.degrees( np.arcsin( (1e-6*blorder*n*centerWavelength) / 
                             (2. * np.cos(-1.*np.radians(theta))) ) ) - theta
#  alpha = np.degrees( np.arcsin((1e-6*blorder*n*centerWavelength) -
#                               (np.sin(np.radians(beta))) ) )
  dlamdx = (1e3*np.cos(np.radians(beta))) / (blorder*n*Lb)
  rangeArray =  dlamdx*pixel_size*spectral_pixels
#  gratingAngle = alpha-(alpha-beta)/2.
#  gratingAngle = beta + theta
  
  
  maxW = centerWavelength+ fraction*rangeArray
  minW = centerWavelength- fraction*rangeArray
  maxBeta = np.degrees( np.arcsin( (1e-6*blorder*n*maxW) / 
                             (2. * np.cos(-1.*np.radians(theta))) ) ) - theta
  minBeta = np.degrees( np.arcsin( (1e-6*blorder*n*minW) / 
                             (2. * np.cos(-1.*np.radians(theta))) ) ) - theta
  maxGratingAngle = maxBeta + theta
  minGratingAngle = minBeta + theta
  
  t = np.arange(nrSteps)
  exSum = 0
  for i in range(nrSteps):
    exSum = exSum + np.exp(i*growthFactor)
  
  stepScaling = (maxGratingAngle-minGratingAngle)/exSum
#  print(minBeta,maxBeta,minGratingAngle,maxGratingAngle,stepScaling)
#  print(maxW,maxBeta,maxGratingAngle,stepScaling)
  angle = np.zeros(nrSteps+1)
  angle[0] = minGratingAngle
  
  for i in np.arange(1,nrSteps+1):
    angle[i] = angle[i-1] + stepScaling*np.exp(growthFactor*t[i-1])
  
    
  return np.float32(angle)

#a = logStepAlign3(1083.,5,0.1)
#print(a)
#fig, ax=plt.subplots()
#ax.plot(a)
#plt.show()