#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:14:01 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
Revision history
----------------
15 August 2018:
    Removed instrument dark as it is contained in the background dark already.
    Andre
"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *
from datetime import datetime

def ddOne(data,
          linThreshold,
          mode,
          backgroundDark,
          gainTable,
          waveCal,
          beamMapping,
          badPixels,
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
        output from GAIN3 calibration task
        normalized gain table to correct gain variation (multiplicative)
    waveCal : (2048, 2048) ndarray, float32
        output from WAVECAL calibration task describing what the wavelength
        at each pixel is
    beamMapping : (2, 2048, 2048) ndarray, float32
        output from ALIGN3 calibration task
        information describing how left and right side of beam match.
        
    Returns
    -------
    result : (2048, 2048) ndarray, float32
        Fully corrected image

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
  assert(data.shape == backgroundDark.shape),\
  "background dark array has not the same dimension as the data cube"
  # check data types
  assert(backgroundDark.dtype == np.float32()),\
  "background dark data type is not float32"
  assert(gainTable.dtype == np.float32()),\
  "gain table data type is not float32"
  # turn into float for what is to come
  data = np.float32(data)
  
  ################# 2. reference pixel correction #############################
  # needed for slow mode only
  if mode is "SLOW":
    #do something
    data = data
    
  if debug:
    try:
      file.write("camera mode is "+mode+"\n")
    except:
      print("camera mode is "+mode+"\n")
    
  
  ################# 3. make the linearity correction ##########################  
  # in this case we use the NDR number as our 'time' axis
  dataTime = np.arange(data.shape[0])
  # Do quadratic fit, use data up to threshold
  coef = cnPolyfit(data, 2, mode, linThreshold)
  nonLinear = np.multiply(coef[0,:,:],dataTime[:,None,None]**2.)
  linearityCorrected = data - nonLinear
  linearityCorrectedTwo = coef[1,:,:]*dataTime[-1]
  
  #NOTE: linearized last - first frame is not exactly the same as
  #      linear coefficient * time
#  print(linearityCorrected[0,0,0]-linearityCorrected[-1,0,0])
#  print(coef[1,0,0]*dataTime[-1])
  
  
  ################# 4. subtract the background ################################  
  backgroundSubtracted = linearityCorrected-backgroundDark
  
  
  
  ################# 5. flat fielding ##########################################
  gainCorrected = backgroundSubtracted*gainTable
  
#  fig, ax=plt.subplots()
#  plt.imshow(backgroundSubtracted[-1,:,:])
#  fig, ax=plt.subplots()
#  plt.imshow(gainTable)
  
  # for slow mode signal is inverted
  if mode is "SLOW":
    result = gainCorrected[0,:,:] - gainCorrected[-1,:,:]
  else:
    result = gainCorrected[-1,:,:] - gainCorrected[0,:,:]
  

  ################# 6. wavecal and beam mapping ###############################
  
  if debug:
    try:
      file.write("result maximum is "+str(np.max(result))+"\n")
      file.write("result minimum is "+str(np.min(result))+"\n")
      file.write("result median "+str(np.median(result))+"\n")
    except:
      print("result maximum is "+str(np.max(result)))
      print("result minimum is "+str(np.min(result)))
      print("result median "+str(np.median(result)))
  
  return result
  
linThreshold = 0
mode = "SLOW"


waveCal = np.ones((10,10),dtype="float32")
beamMapping = np.ones((2,10,10),dtype="float32")
badPixels = np.ones((10,10),dtype="uint8")

## reading the data
#a=cnH2rgRamps("data/spectra/simSpectrum*",
#              "fits",readMode="SLOW",subArray=None,verbose=True)
#data = np.squeeze(a.read("fits",dtype=np.uint16))
## reading the background data
#b=cnH2rgRamps("data/backgroundDark/masterBackgroundDark*",
#              "arr",readMode="SLOW",subArray=None,verbose=True)
#backgroundDark = np.squeeze(b.read("arr",dtype=np.float32))
#gainTable = in_im = fits.open("data/gain/gain3Table.fits")[0].data.astype(np.float32)

result = ddOne(data,
               linThreshold,
               mode,
               backgroundDark,
               gainTable,
               waveCal,
               beamMapping,
               badPixels,
               debug=True,
               logPath=None)

fig, ax=plt.subplots()
plt.imshow(data[0]-data[-1], vmin=0.,vmax=40000.)
fig, ax=plt.subplots()
plt.imshow(result, vmin=0.,vmax=40000.)