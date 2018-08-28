#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------

28 August 2018:
    Changed file and function name to match other calibration functions.
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *

def calInstrumentDark(data, writeToFile=False, path=None, sequenceName=None,
                         fileFormat='fits'):
  
  """
  Returns the averaged instrument dark ramp of the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    writeToFile : bolean, optional, default=False
        writing to fits files for testing
    path: string
        directory to which optinal fits files will be written to
    sequenceName: string
        name of sequence to which optinal fits files will be written to
    fileFormat: string ('fits' | 'arr')
        writing fits files or little endian binary arrays

    Returns
    -------
    averagedDark : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube of master instrument dark ramp.

    Raises
    ------

    Notes
    -----
    
    CSS uses little endian (16 or 32 bits). BTD converts to big endian for 
    network transport. All simulated data should be little endian.
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> data = np.zeros((3,5,10,10),dtype='uint16')+6
    >>> InstrumentDark = calInstrumentDark(data, writeToFits=True,
                                              path='data/instrumentDark/',
                                              sequenceName='masterInstrumentDark')
   """
  
  
  #TODO: reference pixel correction should be done prior to averaging. Need to check when and if to do it.
  #TODO: Linearity correction prior or after averaging
  #print('mean variance and std of variance', np.mean(np.var(data,axis=0)),np.std(np.var(data,axis=0)))
  averagedDark = np.float32(np.average(data, axis=0))
  
  # for test purposes lets write the files to fits
  if writeToFile:
    for i in range(averagedDark.shape[0]):
      if fileFormat == 'fits':
        hdu = fits.PrimaryHDU(averagedDark[i])
        hdu.writeto(path+sequenceName+'-{0:05d}-{1:05d}'.format(0,i)+'.fits',overwrite=True)
      else:
        averagedDark[i].tofile(path+sequenceName+'-{0:05d}-{1:05d}'.format(0,i)+'.arr',sep="")
  
  return averagedDark

#
#a=cnH2rgRamps("data/instrumentDark/simInstrumentDark*",readMode="SLOW",subArray=None,verbose=True)
#data=a.read()
#b= calInstrumentDark(data, writeToFile=True, path='data/instrumentDark/',
#                        sequenceName='masterInstrumentDark', fileFormat='arr')