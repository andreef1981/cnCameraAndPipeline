#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
15 August 2018:
    Removed instrument dark subtractiona as we only use background dark stand
    alone and thus the instrument dark has to be contained. Andre
28 August 2018:
    Changed file and function name to match other calibration functions.
24 September 2018:
    - added debuging and logging
    - now use loop and cnPolyfit for linear correction
    - added mode, linThreshold inputs
    
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *
from helperFunctions import *


def calBackgroundDark(data,
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

    >>> data = np.zeros((3,5,10,10),dtype='uint16')+6
    >>> InstrumentDark = np.zeros((5,10,10),dtype='float32')+3
    >>> BackgroundDark = calBackgroundDark(data, InstrumentDark,
                                              writeToFits=True,
                                              path='data/backgroundDark/',
                                              sequenceName='masterBackgroundDark')
   """
  
  if not logPath is None:
#    file = open(logPath+"ddOne_"+datetime.now().strftime("%Y-%m-%d_%H:%M:%S")+".log", 'w')
    file = open(logPath+"calBackgroundDark.log", 'w')
 
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
    
    # nonLinear = np.multiply(coef[0,:,:],dataTime[:,None,None]**2.)
    # linearityCorrected[i,:,:,:] = data[i,:,:,:] - nonLinear
    # this returns the total linearize flux
    if order == 1:
      linearityCorrected[i,:,:,:] = np.multiply(coef[0,:,:],dataTime[:,None,None])
    else:
      linearityCorrected[i,:,:,:] = np.multiply(coef[1,:,:],dataTime[:,None,None])
    
      
      
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



# b=cnH2rgRamps("data/coronalObs-sensitivity/spBackgroundDark",
#               "fits",readMode="SLOW",subArray=None,verbose=True,cssStyle=True,
#               ramps=3, ndr=2)
# data=b.read("fits",dtype=np.uint16)
# linThreshold = 0
# mode = "SLOW"
# c= calBackgroundDark(data,
#                       linThreshold,
#                       mode,
#                       debug=True,
#                       logPath=None,
#                       writeToFile=True,
#                       filePath='data/coronalObs-sensitivity/',
#                       sequenceName='spMasterBackgroundDark',
#                       fileFormat="both")

# #%%
# print(np.min(c),np.max(c))
# fig, ax=plt.subplots()
# ax.imshow(c[1])
