#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:46:37 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------

"""
import numpy as np
# from astropy.io import fits
# from cnPipeline import *
from helperFunctions import *

def calAlign1(data,
              stagePositionFM1A,
              stagePositionFM1B,
              camera,
              backgroundDark,
              gainTable,
              badPixels,
              oldAlignFM1A,
              oldAlignFM1B,
              oldStepSizeFM1A,
              oldStepSizeFM1B,
              oldRotationFM1A,
              oldRotationFM1B,
              linThreshold,
              mode,
              changeThreshold=None,
              simulateChange=np.zeros((2,3),dtype=bool),
              debug=False,
              logPath=None,
              writeToFile=False,
              filePath=None,
              sequenceName=None,
              fileFormat='fits'):

  #???: does the dark need to be a "GOS thermal dark"?
  
  """
  Returns the best center position of the beam steering mirror FM1A and FM1B 
  for a given slit mask, pickoff state. 
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that will be used to determine the best alignment.
        (raster scan)
    stagePositionFM1A : (#ramps,) ndarray, float32
        1D array that contains the positions of the FM1A stage for each ramp.
    stagePositionFM1B : (#ramps,) ndarray, float32
        1D array that contains the positions of the FM1B stage for each ramp.
    camera : string
        string describing which instrument is used "SP" or "CI"
    backgroundDark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gainTable: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldAlignFM1A: float32
        the previously determined best alignment for FM1A
    oldAlignFM1B: float32
        the previously determined best alignment for FM1B
    oldStepSizeFM1A: float32
        the previously determined step size for FM1A
    oldStepSizeFM1B: float32
        the previously determined step size for FM1B
    oldRotationsFM1A: float32
        the previously determined rotation for FM1A
    oldRotationFM1B: float32
        the previously determined rotation for FM1B
    linThreshold : float
        threshold for the quadratic fit. Everything above will be excluded.
    mode : string, "SLOW", "FAST", "LineByLine"
        defines the readout mode of the camera
    changeThreshold : (2,3) ndarray, float32/uint64
        change beyond which change flag is set
    simulateChange: (2,3) ndarray, boolean, default=False
        boolean array to simulate change in alginments
    debug : boolean, default False
        printing debug messages 
    logPath : string, default None
        if None, printing to standard output, if path provided write to 
        calAlign.log
    writeToFile : boolean, optional, default=False
        writing to fits files for testing
    path: string
        directory to which optinal fits files will be written to
    sequenceName: string
        name of sequence to which optinal fits files will be written to
    fileFormat: string ('fits' | 'arr')
        writing fits files or little endian binary arrays


    Returns
    -------
    newAlignFM1A : float32/int64
        the new best alignment for FM1A
    newAlignFM1B : float32/int64
        the new best alignment for FM1B
    newStepSizeFM1A: float32
        the previously determined step size for FM1A
    newStepSizeFM1B: float32
        the previously determined step size for FM1B
    newRotationsFM1A: float32
        the previously determined rotation for FM1A
    newRotationFM1B: float32
        the previously determined rotation for FM1B
    changeFlag: (2,3) ndarray, boolean
        boolean array to indicat changes in FM1A or FM1B

    Raises
    ------
    
    Notes
    -----
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> 
   """
   # initialize the change flag to all False
  changeFlag=np.zeros((2,3),dtype=bool)
   
  if not logPath is None:
    file = open(logPath+"calAlign1_"+datetime.now().strftime("%Y-%m-%d_%H:%M:%S")+".log", 'w')
    # file = open(logPath+"calAlign1.log", 'w')
 
  ################# 1. make some checks and conversions #######################
  
  
  # turn into float for what is to come
  data = np.float32(data)
  
  ################# 2. reference pixel correction #############################
   #TODO: reference pixel correction should be done prior to averaging. Need to check when and if to do it.
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
  if len(data.shape)==4:
    linearityCorrected=cnNonLinearityCorrection(data,mode,linThreshold,multiRamp=True)
  else:
    linearityCorrected=cnNonLinearityCorrection(data,mode,linThreshold,multiRamp=False)
  
  
  ################# 4. subtract the background ################################  
  backgroundSubtracted = linearityCorrected-backgroundDark
  
  
  ################# 5. flat fielding ##########################################
  gainCorrected = backgroundSubtracted*gainTable
  
  # if len(data.shape)==4:
  #   gainCorrected = np.average(gainCorrected, axis=0)
  # for slow mode signal is inverted
  if mode is "SLOW":
    result = gainCorrected[:,0,:,:] - gainCorrected[:,-1,:,:]
  else:
    result = gainCorrected[:,-1,:,:] - gainCorrected[:,0,:,:]
  
  # Fourth ignore/interpolate bad pixels
  
  # Fifth decide which instrument is used and loop through sequences to find
  # best alignment
  # TODO: helper function to find brightest spectra for SP (what is the GOS
  #       pinhole image size on the SP slit plane?) that is centered on slit
  #       image on array (data is raster scan)
  # TODO: for SP determine rotation of slit relative to steering mirror axes.
  # TODO: for CI; helper function to find FM1 positions that center
  #       GOS pinhole on array 
  # TODO: for CI determine rotation of array relative to steering mirror axes.
  # TODO: for both calibrate step size to pixels
  
  # Sixth determine if alignment has changed
  # TODO: what is the threshold for a change
  
  if simulateChange[0,0]:
    newAlignFM1A = oldAlignFM1A+0.05 # assuming user units simulate change by 0.05 degrees
    changeFlag[0,0] = True
  else:
    newAlignFM1A = oldAlignFM1A
  if simulateChange[1,0]:
    newAlignFM1B = oldAlignFM1B-0.05 # assuming user units simulate change by 0.05 degrees
    changeFlag[1,0] = True
  else:
    newAlignFM1B = oldAlignFM1B
  if simulateChange[0,1]:
    newStepSizeFM1A = oldStepSizeFM1A+3.0 # simulate change by 3 mu per pixel
    changeFlag[0,1] = True
  else:
    newStepSizeFM1A = oldStepSizeFM1A
  if simulateChange[1,1]:
    newStepSizeFM1B = oldStepSizeFM1B-3.0 # simulate change by 3 mu per pixel
    changeFlag[1,1] = True
  else:
    newStepSizeFM1B = oldStepSizeFM1B
  if simulateChange[0,2]:
    newRotationFM1A = oldRotationFM1A+0.8 # simulate change by 0.8 degrees
    changeFlag[0,2] = True
  else:
    newRotationFM1A = oldRotationFM1A
  if simulateChange[1,2]:
    newRotationFM1B = oldRotationFM1B-0.8 # simulate change by 0.8 degrees
    changeFlag[1,2] = True
  else:
    newRotationFM1B = oldRotationFM1B
  
  return newAlignFM1A, newAlignFM1B, newStepSizeFM1A, newStepSizeFM1B,\
         newRotationFM1A, newRotationFM1B, changeFlag

