#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:14:01 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
Revision history
----------------
20 September 2018
  created a single plugin for the first part of the detailed display reduction
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *
from datetime import datetime
from astropy.io import fits

def ddOne(data,
          linThreshold,
          mode,
          backgroundDark,
          gainTable,
          badPixels,
          camera,
          waveCal=None,
          beamMapping=None,
          debug=False,
          logPath=None):
  
  #!!!: linThreshold for SLOW, FAST and RSTRDRD mode should live in a property database
  #!!!: bad pixel mask should live in the bulk data base
  #!!!: background data should be linearity corrected
  """
  Unified detailed display plugin to perform linearity, background, gain and
  wavelength correction
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be corrected.
    linThreshold : float, default=0 for slow mode
        threshold for the quadratic fit. Everything above will be excluded.
    mode : string, default="SLOW"
        defines the readout mode of the camera
    backgroundDark : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that contains the background dark ramp with the exact same
        camera setup as the data cube. Has to be linearity corrected.
    gainTable : (2048, 2048) ndarray, float32
        output from GAIN calibration task
        normalized gain table to correct gain variation (multiplicative)
    badPixels : (2048, 2048) ndarray, uint8
        array containing the good and bad pixels
    camera : string
        string describing which instrument is used "SP" or "CI"
    waveCal : (2048, 2048) ndarray, float32, default None
        output from WAVECAL calibration task describing what the wavelength
        at each pixel is
        Only applicable to SP
    beamMapping : (2, 2048, 2048) ndarray, float32, default None
        output from ALIGN3 calibration task
        information describing how left and right side of beam match.
        Only applicable to SP
    debug : boolean, default False
        printing debug messages 
    logPath : string, default None
        if None, printing to standard output, if path provided write to ddOne.log
        
        
    Returns
    -------
    result : (2048, 2048) ndarray, float32
        Fully corrected image
    waveVector : (2048,) ndarray, float32
        vector containing the wavelength at each pixel

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

    >>> data = np.zeros((5,10,10),dtype='uint16')+6
    >>> instrumentDark = np.zeros((5,10,10),dtype='float32')+1
    >>> backgroundDark = np.zeros((5,10,10),dtype='float32')+2
    >>> darkSubtracted = dark(data,instrumentDark,backgroundDark)
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
  assert(backgroundDark.dtype == np.float32()),\
  "background dark data type is not float32"
  assert(gainTable.dtype == np.float32()),\
  "gain table data type is not float32"
  # turn into float for what is to come
  data = np.float32(data)
  
  ################# 2. reference pixel correction #############################
  # needed for slow mode only?
  if mode is "SLOW":
    #do something
    data = data
    
  if debug:
    try:
      file.write("camera mode is "+mode+"\n")
    except:
      print("camera mode is "+mode+"\n")
    
  
  ################# 3. make the linearity correction ##########################  
  
  linearityCorrected=cnNonLinearityCorrection(data,mode,linThreshold,multiRamp=False)
  
  ################# 4. subtract the background ################################  
  backgroundSubtracted = linearityCorrected-backgroundDark
  
  
  
  ################# 5. flat fielding ##########################################
  gainCorrected = backgroundSubtracted*gainTable
  
#  fig, ax=plt.subplots()
#  plt.imshow(backgroundSubtracted[-1,:,:])
#  fig, ax=plt.subplots()
#  plt.imshow(gainTable)
  
  # for slow mode signal is inverted
  if mode == "SLOW":
    result = gainCorrected[0,:,:] - gainCorrected[-1,:,:]
  else:
    result = gainCorrected[-1,:,:] - gainCorrected[0,:,:]
  

  ################# 6. wavecal and beam mapping ###############################
  if camera == "SP":
    # for now just use the center row
    waveVector = waveCal[1024,:]
  else:
    waveVector = 0
    
  if debug:
    try:
      file.write("result maximum is "+str(np.max(result))+"\n")
      file.write("result minimum is "+str(np.min(result))+"\n")
      file.write("result median "+str(np.median(result))+"\n")
    except:
      print("result maximum is "+str(np.max(result)))
      print("result minimum is "+str(np.min(result)))
      print("result median "+str(np.median(result)))
  
  return result, waveVector
  

######run this first with second part commented because spSpectrum sequence
# # is too long for a single call to ddOne

# # reading the data
# # cssStyle needs ramps and ndr information
# a=cnH2rgRamps("data/coronalObs-sensitivity/spObserve.",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=8, ndr=2)
# data = np.squeeze(a.read("fits",dtype=np.uint16))
# # reading the background data
# b=cnH2rgRamps("data/calibrations/spSlow2-masterBackgroundDark",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=2)
# backgroundDark = np.squeeze(b.read("fits",dtype=np.float32))

# gainTable = in_im = fits.open("data/calibrations/spMasterGain3-slow.000.fits")[0].data.astype(np.float32)

# waveCal = fits.open("data/calibrations/spMasterWavecal.000.fits")[0].data.astype(np.float32)

# badPixels = fits.open("data/calibrations/spBadPixels.000.fits")[0].data.astype(np.uint8)

# c=cnH2rgRamps("data/calibrations/spBeamMapping",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=2)
# beamMapping = np.squeeze(c.read("fits",dtype=np.float32))
# #data = np.squeeze(data[0,:,:,:])


# ######### run this after first part with first part commented
# linThreshold = 0
# mode = "SLOW"
# camera = "SP"
# fileFormat = "both"
# filePath = "data/coronalObs-sensitivity/"
# sequenceName = "spObserve-ddOneProcessed"

# for ii in range(8):
#   result,wavevector = ddOne(data[ii,:,:,:],
#                             linThreshold,
#                             mode,
#                             backgroundDark,
#                             gainTable,
#                             badPixels,
#                             camera,
#                             waveCal,
#                             beamMapping,
#                             debug=True,
#                             logPath=None)
#   if fileFormat == "fits":
#     hdu = fits.PrimaryHDU(result)
#     hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(ii)+'.fits',overwrite=True)
#   elif fileFormat == "raw":
#     result.tofile(filePath+sequenceName+'{0:03d}'.format(ii)+'.raw',sep="")
#   elif fileFormat == "both":
#     hdu = fits.PrimaryHDU(result)
#     hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(ii)+'.fits',overwrite=True)
#     result.tofile(filePath+sequenceName+'.{0:03d}'.format(ii)+'.raw',sep="")
#   else:
#     raise ValueError("Specify valid file type.")
  

################################## CI ########################################
#####run this first with second part commented because spSpectrum sequence
# is too long for a single call to ddOne

# # reading the data
# # cssStyle needs ramps and ndr information
# a=cnH2rgRamps("data/coronalObs-sensitivity/ciObserve.",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=72, ndr=4)
# data = np.squeeze(a.read("fits",dtype=np.uint16))
# # reading the background data
# b=cnH2rgRamps("data/calibrations/ciSlow4-masterBackgroundDark",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=4)
# backgroundDark = np.squeeze(b.read("fits",dtype=np.float32))

# gainTable = in_im = fits.open("data/calibrations/ciMasterGain1-slow.000.fits")[0].data.astype(np.float32)

# waveCal = None

# badPixels = fits.open("data/calibrations/ciBadPixels.000.fits")[0].data.astype(np.uint8)


# beamMapping = None
# #data = np.squeeze(data[0,:,:,:])


# ######### run this after first part with first part commented
# linThreshold = 0
# mode = "SLOW"
# camera = "CI"
# fileFormat = "both"
# filePath = "data/coronalObs-sensitivity/"
# sequenceName = "ciObserve-ddOneProcessed"

# for ii in np.arange(32, 40):
#   print(ii)
#   result,wavevector = ddOne(data[ii,:,:,:],
#                             linThreshold,
#                             mode,
#                             backgroundDark,
#                             gainTable,
#                             badPixels,
#                             camera,
#                             waveCal,
#                             beamMapping,
#                             debug=True,
#                             logPath=None)
#   if fileFormat == "fits":
#     hdu = fits.PrimaryHDU(result)
#     hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(ii)+'.fits',overwrite=True)
#   elif fileFormat == "raw":
#     result.tofile(filePath+sequenceName+'{0:03d}'.format(ii)+'.raw',sep="")
#   elif fileFormat == "both":
#     hdu = fits.PrimaryHDU(result)
#     hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(ii)+'.fits',overwrite=True)
#     result.tofile(filePath+sequenceName+'.{0:03d}'.format(ii)+'.raw',sep="")
#   else:
#     raise ValueError("Specify valid file type.")

