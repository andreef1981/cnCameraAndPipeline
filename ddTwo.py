#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:09:36 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
Revision history
----------------
  10 Feb 2019: - implement separate case for CI camera
"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
# from helperFunctions import *
from datetime import datetime
from astropy.io import fits

def ddTwo(data,
          modStates,
          modMode,
          demodMatrix,
          mode,
          badPixels,
          camera,
          waveCal=None,
          debug=False,
          logPath=None):
  
  
  #!!!: bad pixel mask should live in the bulk data base
  #!!!: demodulation matrix has wavelength dependency
  #!!!: Group model matrix is also needed
  """
  Detailed display plugin to perform beam subtraction and demodulation. This 
  assumes that ddOne has performed all background, linearity and gain correction.
  The left and right side of the beams have been mapped and can be subtracted.
  
   Parameters
    ----------
    data : (#modulationStates, 2048, 2048) ndarray, float32
        dD data cube that needs to be demodulated. 
    modStates : (#modulationStates) ndarray, uint8
        sequence of the modulation states from a given start postition, i.e 
        [4,5,6,7,0,1,2,3] for a modulation sequence with 8 states
    modMode : string
        defining the mode in which the modulator was running,
        stepped / continous / off / out
    demodMatrix : (2, 4, #modulationStates) ndarray, float32
        demodulation martices for the selected modulation states, camera mode
        and wavelength dependent. Sum and difference of beam have different
        demodulation
    mode : string
        defines the readout mode of the camera
        SLOW / FAST / RSTRDRD
    badPixels : (2048, 2048) ndarray, uint8
        array containing the good and bad pixels
    camera : string
        string describing which instrument is used "SP" or "CI"
    waveCal : (1024) ndarray, float32
        vector describing what the wavelength at each pixel is
        Only applicable to SP
    debug : boolean, default False
        printing debug messages 
    logPath : string, default None
        if None, printing to standard output, if path provided write to ddOne.log
        
        
    Returns
    -------
    result : (4,1024 , 2048) ndarray, float32
        4 Stokes images (I,Q,U,V) where the spectral direction is 1024 pixels
        wide
    
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

    >>> 
   """
  
  #TODO: Should we implement a check for exposure time, camera mode or is the
  #      check for the right calibration data instrumentDark and backgroundDark
  #      done before calling plugin.
  if not logPath is None:
#    file = open(logPath+"ddOne_"+datetime.now().strftime("%Y-%m-%d_%H:%M:%S")+".log", 'w')
    file = open(logPath+"ddOne.log", 'w')
#    file.write("your friendly logfile\n for the ddOne plugin\n")

  ################# 1. make some checks and conversions #######################
  # make sure the data cubes have the same dimensions
  #TODO: fast mode has shorter sequence 
#  assert(data.shape == backgroundDark.shape),\
#  "background dark array has not the same dimension as the data cube"
  # check data types
  
  
  ################# 2. reordering data to match demodulation matrix ###########
  nrModulationStates = np.size(modStates)
  back = np.mod(np.arange(nrModulationStates)+np.where(modStates==0)[0],nrModulationStates)
  sortedData = data[back,:,:]
    
  if debug:
    try:
      file.write("data has been sorted\n")
    except:
      print("data has been sorted\n")
    
  if camera == "SP":
    ################# 3. beam subtraction #######################################  
    left = sortedData[:,:,:1024]
    right = sortedData[:,:,1024:]
    summ = (left+right)/2.
    diff = (left-right)/2.
    # reshaping so that matrix multiplication can be used
    summ = np.reshape(summ,(nrModulationStates,2048*1024))
    diff = np.reshape(diff, (nrModulationStates,2048*1024))
    ################# 4. demodulation ###########################################  
    recovI = np.matmul(demodMatrix[0,:,:],summ)
    recovQUV = np.matmul(demodMatrix[1,:,:],diff)
    
    recovI = np.reshape(recovI,(4,2048,1024))
    recovQUV = np.reshape(recovQUV,(4,2048,1024))
    stokes = np.concatenate((recovI[0,:,:][None,],recovQUV[1:,:,:]))
    waveVector = waveCal
    
  elif camera == "CI":
    temp = np.reshape(data,(nrModulationStates,2048*2048))
    recov = np.matmul(demodMatrix,temp)
    stokes = np.reshape(recov,(4,2048,2048))
    waveVector = None
    
  return np.float32(stokes), waveVector
  

# # SP
# ## reading the data
# # cssStyle needs ramps and ndr information
# a=cnH2rgRamps("data/coronalObs-sensitivity/spObserve-ddOne",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=8)
# data = np.squeeze(a.read("fits",dtype=np.float32))

# waveCal = fits.open("data/calibrations/spMasterWavecal.000.fits")[0].data.astype(np.float32)

# badPixels = fits.open("data/calibrations/spBadPixels.000.fits")[0].data.astype(np.uint8)


# ######### run this after first part with first part commented
# modMode = "stepped"
# #modStates = np.arange(8)

# modStates = np.array([3,4,5,6,7,0,1,2])

# data = data[modStates,:,:]

# mode = "SLOW"
# camera = "SP"

# demodMatrix = np.load("data/spectra/demod-8-SiIX.npy")

# result, waveVector = ddTwo(data,
#                            modStates,
#                            modMode,
#                            demodMatrix,
#                            mode,
#                            badPixels,
#                            camera,
#                            waveCal=None,
#                            debug=False,
#                            logPath=None)

# #fig, ax=plt.subplots()
# #plt.imshow(data[0]-data[-1], vmin=0.,vmax=40000.)
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
# ax1.imshow(result[0,:,:], vmin=0.,vmax=3000.)
# ax2.imshow(result[1,:,:], vmin=0.,vmax=200.)
# ax3.imshow(result[2,:,:], vmin=0.,vmax=200.)
# ax4.imshow(result[3,:,:], vmin=0.,vmax=20.)
  

##########################################################################
#   # CI
# ## reading the data
# # cssStyle needs ramps and ndr information
# a=cnH2rgRamps("data/coronalObs-sensitivity/ciObserve-ddOne",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=8)
# data = np.squeeze(a.read("fits",dtype=np.float32))

# waveCal = None

# badPixels = fits.open("data/calibrations/ciBadPixels.000.fits")[0].data.astype(np.uint8)


# ######### run this after first part with first part commented
# modMode = "stepped"
# #modStates = np.arange(8)

# modStates = np.array([3,4,5,6,7,0,1,2])

# data = data[modStates,:,:]

# mode = "SLOW"
# camera = "CI"

# demodMatrix = np.load("data/spectra/demod-8-SiX-CI.npy")

# result, waveVector = ddTwo(data,
#                            modStates,
#                            modMode,
#                            demodMatrix,
#                            mode,
#                            badPixels,
#                            camera,
#                            waveCal=None,
#                            debug=False,
#                            logPath=None)
# #%%
# #fig, ax=plt.subplots()
# #plt.imshow(data[0]-data[-1], vmin=0.,vmax=40000.)
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
# ax1.imshow(result[0,:,:], vmin=0.,vmax=8000.)
# ax2.imshow(result[1,:,:], vmin=0.,vmax=200.)
# ax3.imshow(result[2,:,:], vmin=-200.,vmax=200.)
# ax4.imshow(result[3,:,:], vmin=-20.,vmax=20.)