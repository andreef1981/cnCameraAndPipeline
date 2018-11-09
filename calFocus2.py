#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:45:58 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------

    
"""
import numpy as np
from helperFunctions import *
from astropy.io import fits
from cnPipeline import *
from ddLinearity import *
from calBackgroundDark import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def calFocus2(data,
              stagePosition,
              backgroundDark,
              gainTable,
              badPixels,
              oldFocus,
              mode,
              linThreshold,
              changeThreshold,
              camera,
              debug=False,
              logPath=None,
              simulateChange=False,
              writeToFile=False,
              path=None,
              sequenceName=None,
              fileFormat='fits'):
  
  #TODO: 
  
  #TODO: should we return the focus position in user or raw units
  #      (float32/uint64 for positions)
  #TODO: should we add a threshold in a property database to decide if focus
  #      has changed?
  #TODO: Where do we keep the bad pixel mask stored?
  
  
  #TODO: where is delta calculation to base done?
  
  """
  Returns the best focus position for the external warm focus FM2.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that will be used to determine the best focus position.
    stagePosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the focus stage for each ramp.
    backgroundDark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gainTable: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldFocus: float32/int64
        the prevously determined focus position
    mode : string, default="SLOW"
        defines the readout mode of the camera
    linThreshold : float
        threshold for the quadratic fit. Everything above/below will be excluded.
    changeThreshold : float32
        change beyond which change flag is set
    camera : string
        string describing which instrument is used "SP" or "CI"
    simulateChange: bolean
        test flag
    writeToFile : bolean, optional, default=False
        writing to files for testing
    path: string
        directory to which optinal fits files will be written to
    sequenceName: string
        name of sequence to which optinal fits files will be written to
    fileFormat: string ('fits' | 'arr')
        writing fits files or little endian binary arrays
        
    Returns
    -------
    newFocus : float32/int64
        the new focus position for SM5
    changeFlag: bolean
        indication whether focus has changed

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
  if not logPath is None:
    file = open(logPath+"wavecal_"+datetime.now().strftime("%Y-%m-%d_%H:%M:%S")+".log", 'w')
    # file = open(logPath+"wavecal.log", 'w')
 
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
  
  if debug:
    fig, ax=plt.subplots(num=2)
    im=plt.imshow(linearityCorrected[4,4])#, norm=LogNorm(vmin=100, vmax=10000))
    fig.colorbar(im)
    plt.show()
  ################# 4. subtract the background ################################  
  backgroundSubtracted = linearityCorrected-backgroundDark
  
  if debug:
    fig, ax=plt.subplots(num=3)
    im=plt.imshow(backgroundSubtracted[4,4])#, norm=LogNorm(vmin=100, vmax=10000))
    fig.colorbar(im)
    plt.show()
  
  ################# 5. flat fielding ##########################################
  gainCorrected = backgroundSubtracted*gainTable
  if debug:
    fig, ax=plt.subplots(num=4)
    im=plt.imshow(gainCorrected[4,4])#, norm=LogNorm(vmin=100, vmax=10000))
    fig.colorbar(im)
    plt.show()
  # if len(data.shape)==4:
  #   gainCorrected = np.average(gainCorrected, axis=0)
  # for slow mode signal is inverted
  if mode is "SLOW":
    result = gainCorrected[:,0,:,:] - gainCorrected[:,-1,:,:]
  else:
    result = gainCorrected[:,-1,:,:] - gainCorrected[:,0,:,:]
  ################# 6.beam mapping ############################################
  
  #TODO: beam mapping
  
  ################# 7. find spatial & spectral profiles #######################
  #TODO: pinhole mask implementation
  #TODO: spatial focus metric
  if debug:
    fig, ax=plt.subplots(num=5)
    im=plt.imshow(result[4])#, norm=LogNorm(vmin=100, vmax=10000))
    fig.colorbar(im)
    plt.show()
    
  
  # find the spatial peaks
  spatialPeaksLeft = []
  spatialWidthsLeft = []
  for i in range(result.shape[0]):
    ind, widths, prominences, heights, profile = cnFindSpatialPeaks(result[i,:,:],
                                                                    [3,50],
                                                                    [10,980],
                                                                    None,
                                                                    [100, 20000],
                                                                    [0,20000],
                                                                    invert=False)
                                                                    # [100, 0.9*np.max(-1*result[i,:1024,:])-np.mean(-1*result[i,:1024,:])],
                                                                    # [np.mean(-1*result[i,:1024,:]),0],
                                                                    # invert=True)
    
    # There should be only one peak from GOS pinhole
    if debug:
      print(ind, 'widths')
      print(widths, 'widths')
    spatialPeaksLeft.append(int(ind))
    spatialWidthsLeft.append(int(widths))
  spatialPeaksRight = spatialPeaksLeft
  spatialWidthsRight = spatialWidthsLeft
  
  avLeft = np.zeros(result.shape[0],dtype=np.float16())
  avRight = np.zeros(result.shape[0],dtype=np.float16())
  
  for i in range(result.shape[0]):
    if debug:
      print(spatialPeaksLeft[i],spatialWidthsLeft[i],'ind,width')
    # left side
    spectralPeaks, leftProfile = cnFindSpectralPeaks(result[i,:,:],
                                      [2,100],
                                      spatialPeaksLeft[i],
                                      spatialWidthsLeft[i],
                                      [0,1024],
                                      [100,20000],
                                      [-20000,-50],
                                      invert=True)
    
    
    avLeft[i] = np.mean(spectralPeaks[0][1])
    # right side
    spectralPeaks, rightProfile = cnFindSpectralPeaks(result[i,:,:],
                                      [2,100],
                                      spatialPeaksRight[i],
                                      spatialWidthsRight[i],
                                      [1024,2048],
                                      [100,20000],
                                      [-20000,-50],
                                      invert=True)
    avRight[i] = np.mean(spectralPeaks[0][1])
  if debug:
    print(avLeft,avRight)
  
  #TODO: use both sides to determine best focus
  # fit quadratic function to line width values
  parLeft = np.polyfit(stagePosition,np.float64(avLeft),2)
  hr = np.arange(np.min(stagePosition),np.max(stagePosition),0.001)
  fitterLeft = np.polyval(parLeft,hr)
  indLeft = np.argmin(fitterLeft)
  
  newFocus = hr[indLeft]  

    
  if simulateChange:
    newFocus = oldFocus+5.0 # assuming user units simulate change by 5mm
    
    
  ################# 8. check if focus changed #################################  
  if np.abs(np.median(newFocus-oldFocus)) > changeThreshold:
    changeFlag = True
  else:
    changeFlag = False

  
  # for test purposes lets write the files to fits
  if writeToFile:
    if fileFormat == 'fits':
      hdu = fits.PrimaryHDU(newFocus)
      hdu.writeto(path+sequenceName+'.fits',overwrite=True)
    else:
      newFocus.tofile(path+sequenceName+'.arr',sep="")
  
  return newFocus, changeFlag


# reading the data
# cssStyle needs ramps and ndr information
linThreshold = 0
mode = "SLOW"

a=cnH2rgRamps("data/coronalObs-sensitivity/spFocus2-GOSpinhole",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=9, ndr=5)
data = np.squeeze(a.read("fits",dtype=np.uint16))
# reading the background data
b=cnH2rgRamps("data/coronalObs-sensitivity/spFocus-masterBackground",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=1, ndr=5)

backgroundDark=np.squeeze(b.read("fits",dtype=np.float32))

#%%
gainTable = in_im = fits.open("data/coronalObs-sensitivity/spMasterGain3.000.fits")[0].data.astype(np.float32)

changeThreshold = 0.5
#exact focus
oldFocus = -10.193
#within threshold
#oldFocus = -9.8
#outside threshold
#oldFocus = -9.
camera = "SP"
oldWavecal = fits.open("data/coronalObs-sensitivity/spMasterWavecal-plus10.000.fits")[0].data.astype(np.float32)

badPixels = fits.open("data/coronalObs-sensitivity/spBadPixels.000.fits")[0].data.astype(np.uint8)

c=cnH2rgRamps("data/coronalObs-sensitivity/spBeamMapping",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=1, ndr=2)
beamMapping = np.squeeze(c.read("fits",dtype=np.float32))

#%%
# fig, ax=plt.subplots(num=1)
# im=plt.imshow(data[0,4])#, norm=LogNorm(vmin=100, vmax=10000))
# fig.colorbar(im)
# plt.show()
#%%

stagePositions = np.array([-9.4, -9.6, -9.8, -10., -10.2, -10.4, -10.6, -10.8, -11.])

newFocus, changeFlag = calFocus2(data,
                                    stagePositions,
                                    backgroundDark,
                                    gainTable,
                                    badPixels,
                                    oldFocus,
                                    mode,
                                    linThreshold,
                                    changeThreshold,
                                    camera,
                                    debug=True,
                                    logPath=None,
                                    simulateChange=False,
                                    writeToFile=False,
                                    path=None,
                                    sequenceName=None,
                                    fileFormat='fits')