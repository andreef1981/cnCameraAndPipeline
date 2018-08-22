#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 12:09:44 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""
import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *
from astropy.io import fits  # Reading/writing FITS data

def gain(data,
         gainTable):
  """
  Returns a gain corrected ramp of the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be gain corrected.
    gainTable : (2048, 2048) ndarray, float32
        output from GAIN3 calibration task
        normalized gain table to correct gain variation (multiplicative)

    Returns
    -------
    gainCorrected : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that is gain corrected.

    Raises
    ------

    Notes
    -----
    
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
  
  # broadcasting does the right thing by multiplying all frames with the gain
  # table
  gainCorrected = data*gainTable
  
  return gainCorrected

# sample code to run gain plugin
#a=cnH2rgRamps("data/gain/simGain*",readMode="SLOW",subArray=None,verbose=True)
#data = np.squeeze(a.read(dtype=np.uint16))
#
#gainTable = in_im = fits.open("data/gain/gain3Table.fits")[0].data.astype(np.float32)
#b = gain(data,gainTable)
#
#fig, ax=plt.subplots()
#plt.imshow(data[1,1,:,:])
#
#fig, ax=plt.subplots()
#plt.imshow(b[1,1,:,:])
#plt.show()
