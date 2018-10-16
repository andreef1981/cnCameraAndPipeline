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

def spiralGain1(nrSteps,
                stepScaling=100,
                growthFactor=0.1759,
                angleDivision=4.,
                scale=0.05):
  """
  Function to generate a logarithmic spiral pattern starting at (0,0)
  
   Parameters
    ----------
    nrSteps : uint16
        # of steps along the spiral.
    stepScaling : optional, float, default=10
        optional inpput for scaling the step size
    growthFactor : optional, float, default=0.1759
        optional growth factor for the spiral
    angleDivision : optional, float, default=4
        optional input defining how many points per pi radians (pitch)
    scale : optional, float, default=0.05 arcsec/pixel
        scale describing the mapping from arcsec to pixels
    
        
    Returns
    -------
    x : float16
        x position of beam steering mirror in arcsec
    y : float16
        y position of beam steering mirror in arcsec

    Raises
    ------

    Notes
    -----
    
    Logarithmic spiral inspired by Nautilus
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> x,y = spiralGain1(15)
    >>> fig, ax=plt.subplots()
    >>> ax.plot(x,y)
    >>> plt.show()

   """
  # scale = 0.05 #arcsec per pixel
  # stepScaling=100
  # growthFactor=0.1759
  # angleDivision=4.
  # b = 1.6180
  t = np.arange(nrSteps)*np.pi/angleDivision
  x = stepScaling*np.exp(growthFactor*t)*np.cos(t)
  y = stepScaling*np.exp(growthFactor*t)*np.sin(t)
    
  return np.float16(x*scale),np.float16(y*scale)

# x,y = spiralGain1(15)

# fig, ax=plt.subplots()
# ax.plot(x,y)
# plt.show()