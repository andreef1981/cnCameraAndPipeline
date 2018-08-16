#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(Andre Fehlmann)s
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *

def masterBackgroundDark(data, masterInstrumentDark, writeToFile=False, path=None,
                      sequenceName=None, fileFormat='fits'):
  
  """
  Returns the averaged background dark ramp of the CryoNIRSP H2RG.
  
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
    AssertationError
        If the shapes of the input arrays do not match.

    Notes
    -----
    
    CSS uses little endian (16 or 32 bits). BTD converts to big endian for 
    network transport. All simulated data should be little endian.
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> data = np.zeros((3,5,10,10),dtype='uint16')+6
    >>> masterInstrumentDark = np.zeros((5,10,10),dtype='float32')+3
    >>> masterBackgroundDark = masterBackgroundDark(data, masterInstrumentDark,
                                              writeToFits=True,
                                              path='data/backgroundDark/',
                                              sequenceName='masterBackgroundDark')
   """
  
  
  #TODO: reference pixel correction should be done prior to averaging. Need to check when and if to do it.
  
  # make sure the data cubes have the same dimensions
  assert(data[0,:,:,:].shape  == masterInstrumentDark.shape),\
  'input array dimensions do not match'
  
  # first we need to subtract the master instrument Dark
  tempData = data - masterInstrumentDark
  
  #print('mean variance and std of variance', np.mean(np.var(data,axis=0)),np.std(np.var(data,axis=0)))
  averagedDark = np.average(tempData, axis=0)
  
  # for test purposes lets write the files to fits
  if writeToFile:
    for i in range(averagedDark.shape[0]):
      if fileFormat == 'fits':
        hdu = fits.PrimaryHDU(averagedDark[i])
        hdu.writeto(path+sequenceName+'-{0:05d}-{1:05d}'.format(0,i)+'.fits',overwrite=True)
      else:
        averagedDark[i].tofile(path+sequenceName+'-{0:05d}-{1:05d}'.format(0,i)+'.arr',sep="")
  
  return averagedDark


a=cnH2rgRamps("data/instrumentDark/masterInstrumentDark*",readMode="SLOW",subArray=None,verbose=True)
instrumentDark = np.squeeze(a.read(dtype=np.float32))
b=cnH2rgRamps("data/backgroundDark/simBackgroundDark*",readMode="SLOW",subArray=None,verbose=True)
data=b.read()
c= masterBackgroundDark(data, instrumentDark, writeToFile=True, path='data/backgroundDark/',
                     sequenceName='masterBackgroundDark',fileFormat='fits')