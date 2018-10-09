#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:09:36 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
Revision history
----------------

"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *
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
  """
  Unified detailed display plugin to perform linearity, background, gain and
  wavelength correction
  
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
    demodMatrix : (#modulationStates, 4) ndarray, float32
        demodulation martix for the selected modulation states, camera mode
        and wavelength
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
  
  # in this case we use the NDR number as our 'time' axis
  dataTime = np.arange(data.shape[0])
  
  if mode is not "SLOW":
    # fast mode does not use the first frame
    data = data[1:,:,:]
    # time vector is shorter but must maintain the values
    dataTime = dataTime[1:]
  
  # Check for order of correction
  if len(dataTime) == 2:
    order = 1
  elif len(dataTime) > 2:
    order = 2
  else:
    raise ValueError("sequence to short to apply polyfit")
    
  # Do quadratic fit, use data up to threshold
  coef = cnPolyfit(data, order, mode, linThreshold)
  if order == 2:
    nonLinear = np.multiply(coef[0,:,:],dataTime[:,None,None]**2.)
    linearityCorrected = data - nonLinear
  else:
    linearityCorrected = data 
#  linearityCorrectedTwo = coef[1,:,:]*dataTime[-1]
  
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
  if camera is "SP":
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
  

#######run this first with second part commented because spSpectrum sequence
# is too long for a single call to ddOne

## reading the data
# cssStyle needs ramps and ndr information
#a=cnH2rgRamps("data/coronalObs-sensitivity/ciObserve",
#              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#              ramps=72, ndr=4)
#data = np.squeeze(a.read("fits",dtype=np.uint16))
## reading the background data
#b=cnH2rgRamps("data/coronalObs-sensitivity/ciMasterBackgroundDark",
#              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#              ramps=1, ndr=4)
#backgroundDark = np.squeeze(b.read("fits",dtype=np.float32))
#
#gainTable = in_im = fits.open("data/coronalObs-sensitivity/ciMasterGain1.000.fits")[0].data.astype(np.float32)
#
#waveCal = fits.open("data/coronalObs-sensitivity/spWavecal.000.fits")[0].data.astype(np.float32)
#
#badPixels = fits.open("data/coronalObs-sensitivity/ciBadPixels.000.fits")[0].data.astype(np.uint8)
#
#c=cnH2rgRamps("data/coronalObs-sensitivity/spBeamMapping",
#              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#              ramps=1, ndr=2)
#beamMapping = np.squeeze(c.read("fits",dtype=np.float32))
#data = np.squeeze(data[0,:,:,:])


######### run this after first part with first part commented
linThreshold = 0
mode = "SLOW"
camera = "CI"
result,wavevector = ddOne(data,
                          linThreshold,
                          mode,
                          backgroundDark,
                          gainTable,
                          badPixels,
                          camera,
                          waveCal,
                          beamMapping,
                          debug=True,
                          logPath=None)

#fig, ax=plt.subplots()
#plt.imshow(data[0]-data[-1], vmin=0.,vmax=40000.)
fig, ax=plt.subplots()
plt.imshow(result, vmin=0.,vmax=7000.)