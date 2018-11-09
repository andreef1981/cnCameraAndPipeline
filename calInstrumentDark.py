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

def calInstrumentDark(data,
                      linThreshold,
                      mode,
                      debug=False,
                      logPath=None,
                      writeToFile=False,
                      filePath=None,
                      sequenceName=None,
                      fileFormat='fits'):
  
  """
  Returns the averaged background dark ramp of the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    linThreshold : float, default=0 for slow mode
        threshold for the quadratic fit. Everything above will be excluded.
    mode : string, default="SLOW"
        defines the readout mode of the camera
    debug : boolean, default False
        printing debug messages 
    logPath : string, default None
        if None, printing to standard output, if path provided write to 
        calBackgroundDark.log
    writeToFile : bolean, optional, default=False
        writing to fits files for testing
    filePath: string
        directory to which optinal fits files will be written to
    sequenceName: string
        name of sequence to which optinal fits files will be written to
    fileFormat: string ('fits' | 'arr')
        writing fits files or little endian binary arrays

    Returns
    -------
    averagedDark : !!! difference for slow and fast mode
        slow (#NDRs, 2048, 2048) ndarray, float32
        fast (#NDRs-1, 2048, 2048) ndarray, float32
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

    >>>  data = np.zeros((3,5,10,10),dtype='uint16')+6
    >>> InstrumentDark = calInstrumentDark(data, writeToFits=True,
                                              path='data/instrumentDark/',
                                              sequenceName='masterInstrumentDark')
   """
  
  if not logPath is None:
#    file = open(logPath+"ddOne_"+datetime.now().strftime("%Y-%m-%d_%H:%M:%S")+".log", 'w')
    file = open(logPath+"calInstrumentDark.log", 'w')
 
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
      
      
  ################# 4. make the linearity correction ##########################
  averagedDark = np.mean(linearityCorrected, axis=0)
    
    
  if debug:
    try:
      file.write("mean variance is "+str(np.nanmean(np.var(linearityCorrected,axis=0)))+"\n")
      file.write("std of variance is "+str(np.nanstd(np.var(linearityCorrected,axis=0)))+"\n")
    except:
      print("mean variance is "+str(np.nanmean(np.var(linearityCorrected,axis=0))))
      print("std of variance is "+str(np.nanstd(np.var(linearityCorrected,axis=0))))
      
  # for test purposes lets write the files to fits
  if writeToFile:
    for i in range(averagedDark.shape[0]):
      if fileFormat == "fits":
        hdu = fits.PrimaryHDU(averagedDark[i])
        hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(i)+'.fits',overwrite=True)
      elif fileFormat == "raw":
        averagedDark[i].tofile(filePath+sequenceName+'{0:03d}'.format(i)+'.raw',sep="")
      elif fileFormat == "both":
        hdu = fits.PrimaryHDU(averagedDark[i])
        hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(i)+'.fits',overwrite=True)
        averagedDark[i].tofile(filePath+sequenceName+'.{0:03d}'.format(i)+'.raw',sep="")
      else:
        raise ValueError("Specify valid file type.")
  
  return averagedDark

#
#a=cnH2rgRamps("data/instrumentDark/simInstrumentDark*",readMode="SLOW",subArray=None,verbose=True)
#data=a.read()
#b= calInstrumentDark(data, writeToFile=True, path='data/instrumentDark/',
#                        sequenceName='masterInstrumentDark', fileFormat='arr')