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
import matplotlib.pyplot as plt

def logStepAlign3(centerPosition,
                  centerWavelength,
                  nrSteps,
                  stepScaling=100,
                  growthFactor=0.1759):
  #???: are angles the right units to return?
  
  """
  Function to generate a logarithmic step pattern starting at (0,0)
  
   Parameters
    ----------
    centerPosition : float16
        position around which the scan pattern is centered
    centerWavelength : float16
        selected center wavelength for the scan in nm
    nrSteps : uint16
        # of steps along the spiral.
    stepScaling : optional, float, default=10
        optional inpput for scaling the step size
    growthFactor : optional, float, default=0.1759
        optional growth factor for the spiral
    
    
        
    Returns
    -------
    angle : (nrSteps,1) ndarray, float16
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
  # alpha = np.degrees( np.arcsin((1e-6*blorder*n*centerWavelength) -
  #                             (np.sin(np.radians(beta))) ) )
  dlamdx = (1e3*np.cos(np.radians(beta))) / (blorder*n*Lb)
  rangeArray =  dlamdx*pixel_size*spectral_pixels
  print(rangeArray)
  t = np.arange(nrSteps)
  angle = stepScaling*np.exp(growthFactor*t)
  
    
  return np.float16(angle)

a = logStepAlign3(63.3,1083.,15)

fig, ax=plt.subplots()
ax.plot(a)
plt.show()