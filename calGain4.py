#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:46:37 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
from astropy.io import fits
from cnPipeline import *
from ddLinearity import *

def calGain4(data,
             gratingPosition,
             backgroundDark,
             badPixels,
             oldGain,
             threshold=None,
             simulateChange=False,
             mode=None,
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
  Returns the dither gain table for the spectrograph CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    gratingPosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the grating stage for each ramp.
    dark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldGain: (2048, 2048) ndarray, float32
        the prevously determined gain table
    threshold : float32/uint64
        change beyond which change flag is set
    simulateChange: bolean
        test flag
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
    newGain: (2048, 2048) ndarray, float32
        new gain table for multiplicative gain correction
    changeFlag: bolean
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
    file = open(logPath+"calGain4.log", 'w')
 
  ################# 1. make some checks and conversions #######################
  
  
  # turn into float for what is to come
  data = np.float32(data)
  
  ################# 2. reference pixel correction #############################
   #TODO: reference pixel correction should be done prior to averaging. Need to check when and if to do it.
  # needed for slow mode only?
  if mode == "SLOW":
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
  
  ################# 5. make gain table         ################################
  # simple median for now
  temp = np.median(backgroundSubtracted,axis=0)
  # mitigate divide by 0 issue
  gainTable = (temp[-1,:,:]/np.median(temp[-1,:,:]))
  # mitigate divide by 0 issue
  gainTable[gainTable==0] = 1e-15
  gainTable = 1/gainTable
  
  #TODO: fix for some very large gains... should be taken care by bad pixels?
  gainTable = np.where(np.abs(gainTable) > 10, 0, gainTable)
  
  if debug:
    try:
      file.write("mean is "+str(np.mean((newGain)))+"\n")
      file.write("std is "+str(np.std((newGain)))+"\n")
    except:
      print("mean is "+str(np.mean((newGain))))
      print("std is "+str(np.std(newGain)))
  
  if simulateChange:
    newGain = oldGain*1.1 # simulate change by 10%
    changeFlag = True
  else:
    newGain = gainTable
    changeFlag = False
  
  #TODO: Don't know how to compare yet for change. Set change flat always True
  changeFlag = True
  
  return newGain, changeFlag

# ### reading the data
# ## cssStyle needs ramps and ndr information
# linThreshold = 66000
# mode = "FAST"
# a=cnH2rgRamps("data/calibrations/spFast5-masterBackgroundDark",
#               "fits",readMode=mode,subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=5)
# backgroundDark = np.squeeze(a.read("fits",dtype=np.uint16))

# b=cnH2rgRamps("data/calibrations/spGain4",
#               "fits",readMode=mode,subArray=None,verbose=True,cssStyle=True,
#               ramps=10, ndr=5)
# data=b.read("fits",dtype=np.uint16)



# gratingPosition = np.arange(10)
# badPixels = np.zeros((2048,2048),dtype="uint8")
# oldGain = np.ones((2048,2048),dtype="float32")

# c,flag= calGain4(data,
#                  gratingPosition,
#                  backgroundDark,
#                  badPixels,
#                  oldGain,
#                  threshold=None,
#                  simulateChange=False,
#                  mode=None,
#                  debug=False,
#                  logPath=None,
#                  writeToFile=False,
#                  filePath=None,
#                  sequenceName=None,
#                  fileFormat='fits')
  
# fig, ax=plt.subplots()
# plt.imshow(c,vmin=0.9, vmax=1.1)
# plt.show()
