#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
28 August 2018:
  - rename the function to match others
  - add inputs: stage position, dark, badPixels, threshold, simulateChange
  - add ouputs: changeFlag
24 September 2018:
  - 
    
"""
import numpy as np
#from astropy.io import fits
#import matplotlib.pyplot as plt
#from cnPipeline import *
from helperFunctions import *

def calGain3(data,
             stagePosition,
             backgroundDark,
             badPixels,
             oldGain,
             linThreshold,
             mode,
             changeThreshold=None,
             simulateChange=False,
             debug=False,
             logPath=None,
             writeToFile=False,
             filePath=None,
             sequenceName=None,
             fileFormat='fits'):
  #TODO: are the stage positions in user or raw units
  #      (float32/uint64 for positions)
  #TODO: should we add a threshold in a property database to decide if gain
  #      has changed?
  #TODO: Where do we keep the bad pixel mask stored?
  
  
  #TODO: where is delta calculation to base done?
  
  """
  Returns the simple gain table for the spectrograph CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    stagePosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the lamp stage for each ramp.
    backgroundDark : !!! difference in data cube size for slow and fast camera mode
        slow (#NDRs, 2048, 2048) ndarray, float32
        fast (#NDRs-1, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldGain: (2048, 2048) ndarray, float32
        the prevously determined gain table
    linThreshold : float, default=0 for slow mode
        threshold for the quadratic fit. Everything above will be excluded.
    mode : string, default="SLOW"
        defines the readout mode of the camera
    threshold : float32/uint64
        change beyond which change flag is set
    simulateChange: boolean
        test flag
    debug : boolean, default False
        printing debug messages 
    logPath : string, default None
        if None, printing to standard output, if path provided write to 
        calGain3.log
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
    newGain: (2048, 2048) ndarray, float32
        new gain table for multiplicative gain correction
    changeFlag: boolean
        indication whether gain has changed

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

    >>> 
   """
  
  
  if not logPath is None:
#    file = open(logPath+"ddOne_"+datetime.now().strftime("%Y-%m-%d_%H:%M:%S")+".log", 'w')
    file = open(logPath+"calGain3.log", 'w')
 
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
  
  # in this case we use the NDR number as our 'time' axis
  dataTime = np.arange(data.shape[1])
  
  if mode is not "SLOW":
    # fast mode does not use the first frame
    data = data[:,1:,:,:]
    # time vector is shorter but must maintain the values
    dataTime = dataTime[1:]
  
  linearityCorrected = np.zeros(data.shape,dtype=np.float32)
  # need to loop to use cnPolyfit
  
  # Check for order of correction
  if len(dataTime) == 2:
    order = 1
  elif len(dataTime) > 2:
    order = 2
  else:
    raise ValueError("sequence to short to apply polyfit")
    
  for i in np.arange(data.shape[0]):
    # Do quadratic fit, use data up to threshold
    coef = cnPolyfit(np.squeeze(data[i,:,:,:]), order, mode, linThreshold)
    if order == 2:
      nonLinear = np.multiply(coef[0,:,:],dataTime[:,None,None]**2.)
      linearityCorrected[i,:,:,:] = data[i,:,:,:] - nonLinear
    else:
      linearityCorrected[i,:,:,:] = data[i,:,:,:]
  
    
  ################# 4. subtract the background ################################  
  backgroundSubtracted = linearityCorrected-backgroundDark
  
  
  ################# 5. make the gain table ####################################
  temp = np.mean(backgroundSubtracted,axis=0)
  gainTable = np.mean(temp[-1,:,:])/temp[-1,:,:]
  
  if debug:
    try:
      file.write("mean variance is "+str(np.mean(np.var(backgroundSubtracted,axis=0)))+"\n")
      file.write("std of variance is "+str(np.std(np.var(backgroundSubtracted,axis=0)))+"\n")
    except:
      print("mean variance is "+str(np.mean(np.var(backgroundSubtracted,axis=0))))
      print("std of variance is "+str(np.std(np.var(backgroundSubtracted,axis=0))))
      
  # for test purposes lets write the files to fits
  if writeToFile:
    if fileFormat == "fits":
      hdu = fits.PrimaryHDU(gainTable)
      hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(0)+'.fits',overwrite=True)
    elif fileFormat == "raw":
      gainTable.tofile(filePath+sequenceName+'{0:03d}'.format(0)+'.raw',sep="")
    elif fileFormat == "both":
      hdu = fits.PrimaryHDU(gainTable)
      hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(0)+'.fits',overwrite=True)
      gainTable.tofile(filePath+sequenceName+'.{0:03d}'.format(0)+'.raw',sep="")
    else:
      raise ValueError("Specify valid file type.")
  
  if simulateChange:
    newGain = oldGain*1.1 # simulate change by 10%
    changeFlag = True
  else:
    newGain = gainTable
    changeFlag = False
  
  return newGain, changeFlag

### reading the data
## cssStyle needs ramps and ndr information
#a=cnH2rgRamps("data/coronalObs-sensitivity/spMasterBackgroundDark",
#              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#              ramps=1, ndr=2)
#dark = np.squeeze(a.read("fits",dtype=np.uint16))
#
#b=cnH2rgRamps("data/coronalObs-sensitivity/spGain3",
#              "fits",readMode="SLOW",subArray=None,verbose=True,cssStyle=True,
#              ramps=3, ndr=2)
#data=b.read("fits",dtype=np.uint16)
#
#
#
#linThreshold = 0
#mode = "SLOW"
#stagePosition = np.arange(3)
#badPixels = np.zeros((2048,2048),dtype="uint8")
#oldGain = np.ones((2048,2048),dtype="float32")

#c= calGain3(data,
#            stagePosition,
#            dark,
#            badPixels,
#            oldGain,
#            linThreshold,
#            mode,
#            changeThreshold=None,
#            simulateChange=False,
#            debug=True,
#            logPath=None,
#            writeToFile=True,
#            filePath="data/coronalObs-sensitivity/",
#            sequenceName="spMasterGain3",
#            fileFormat="both")
#  
#fig, ax=plt.subplots()
#plt.imshow(c[0],vmin=0.9, vmax=1.1)
#plt.show()
