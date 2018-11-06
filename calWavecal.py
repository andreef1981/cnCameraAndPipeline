#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
    
"""
import numpy as np
from helperFunctions import *

# from astropy.io import fits
# from cnPipeline import *
# from datetime import datetime
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from calBackgroundDark import *

def calWavecal(data,
               centralWavelength,
               gratingPosition,
               slit,
               expectedDispersion,
               backgroundDark,
               gainTable,
               badPixels,
               oldWavecal,
               mode,
               linThreshold,
               changeThreshold,
               debug=False,
               logPath=None,
               simulateChange=False,
               writeToFile=False,
               path=None,
               sequenceName=None,
               fileFormat='fits'):
  
  #???: should we add a threshold in a property database to decide if wavecal
  #      position has changed?
  #???: Where do we keep the bad pixel mask stored?
  #???: does the dark need to be a "lamp off thermal dark"?
  #TODO: Implement algorithm for 52 micron and pinhole masks
  #TODO: Implement algorithm for dealing with multiple grating positions
  """
  Returns wavlength calibration table for the CryoNIRSP H2RG by using a reference
  ThAr spectrum to cross correlate against.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that for wavelength calibration.
    centralWavelength : float32
        the centralWavelength in nanometers
    gratingPosition : float32
        the position of the grating stage in degrees (user units)
    slit : string, "175", "52", "PinHM"
        identifying which slit mask is used 
    expectedDispersion : float32
        the expected dispersion for this grating angle and wavelength in 
        nanometers per pixel
    backgroundDark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gainTable: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldWavecal: (2048, 2048) ndarray, float32
        previous wavelength calibration from calibration store
    mode : string
        defines the readout mode of the camera
    linThreshold : float
        threshold for the quadratic fit. Everything above will be excluded.
    changeThreshold : float32
        change beyond which change flag is set
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
    newWavecal : (2048, 2048) ndarray, float32
        wavelength calibration table
    changeFlag: bolean
        indication whether wavecal has changed

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
  
  
  ################# 4. subtract the background ################################  
  backgroundSubtracted = linearityCorrected-backgroundDark
  
  
  ################# 5. flat fielding ##########################################
  gainCorrected = backgroundSubtracted*gainTable
  
  if len(data.shape)==4:
    gainCorrected = np.average(gainCorrected, axis=0)
  # for slow mode signal is inverted
  if mode is "SLOW":
    result = gainCorrected[0,:,:] - gainCorrected[-1,:,:]
  else:
    result = gainCorrected[-1,:,:] - gainCorrected[0,:,:]
  ################# 6.beam mapping ############################################
  
  #TODO: beam mapping
  
  ################# 7. find spatial profiles###################################
  
  # coronal slit does not require to find center of spectrum
  if slit =="175um":
    spatialPeaks = np.uint16([1024])
    spatialWidths = np.uint16([1024])
    # left side
    spectralPeaks, leftProfile = cnFindSpectralPeaks(result,
                                                    [2,20],
                                                    spatialPeaks,
                                                    spatialWidths,
                                                    [0,1024],
                                                    [np.mean(result),np.max(result)],
                                                    [0,np.max(result)+0.1*np.max(result)],
                                                    invert=False)
    # right side
    spectralPeaks, rightProfile = cnFindSpectralPeaks(result,
                                                    [2,20],
                                                    spatialPeaks,
                                                    spatialWidths,
                                                    [1024,2048],
                                                    [np.mean(result),np.max(result)],
                                                    [0,np.max(result)+0.1*np.max(result)],
                                                    invert=False)
    
    # cross correlate to reference spectrum
    refWl,refSpec = cnRefSpectrum(centralWavelength,"175um",expectedDispersion)

    left = np.squeeze(leftProfile)
    right = np.squeeze(rightProfile)
    
    crossLeft=np.correlate(refSpec,left,mode="same")
    crossRight=np.correlate(refSpec,right,mode="same")
    
    leftInd = np.where(crossLeft==np.max(crossLeft))[0][0]
    rightInd = np.where(crossRight==np.max(crossRight))[0][0]
    
    
    leftWavecal = refWl[leftInd-512:leftInd+512]
    rightWavecal = refWl[rightInd-512:rightInd+512]
    # leftRefspec = refSpec[leftInd-512:leftInd+512]
    # rightRefspec = refSpec[rightInd-512:rightInd+512]
    
    newWavecal = np.tile(np.float32(np.concatenate((leftWavecal,rightWavecal))),(2048,1))
  
  # TODO: remapping the beams would greatly help with detection of edges an fitting
  if debug:
    from matplotlib.ticker import FormatStrFormatter
    # fig, ax=plt.subplots()
    # im=plt.imshow(data[0,0]-data[0,-1], norm=LogNorm(vmin=0.1, vmax=1000))
    # fig.colorbar(im)
    fig, ax=plt.subplots()
    im=plt.imshow(result, norm=LogNorm(vmin=0.1, vmax=1000))
    wo = np.arange(0,1023,128)
    wowo = np.arange(0,2047,128)
    vec = np.concatenate((leftWavecal[wo],rightWavecal[wo]))
    labels = list(map('{:6.2f}'.format,vec))
    plt.xticks(wowo,labels,rotation=-90)
    fig.colorbar(im)
    plt.show()
    
  if simulateChange:
    newWavecal = oldWavecal+10.0 # assuming user units simulate change by 10nm
    
  if np.abs(np.median(newWavecal-oldWavecal)) > changeThreshold:
    changeFlag = True
  else:
    changeFlag = False

  
  # for test purposes lets write the files to fits
  if writeToFile:
    if fileFormat == 'fits':
      hdu = fits.PrimaryHDU(newWavecal)
      hdu.writeto(path+sequenceName+'.fits',overwrite=True)
    else:
      newWavecal.tofile(path+sequenceName+'.arr',sep="")
  
  return newWavecal, changeFlag


# # reading the data
# # cssStyle needs ramps and ndr information
# linThreshold = 0
# mode = "SLOW"

# a=cnH2rgRamps("data/coronalObs-sensitivity/spWavecal175um",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=3, ndr=2)
# data = np.squeeze(a.read("fits",dtype=np.uint16))
# # reading the background data
# b=cnH2rgRamps("data/coronalObs-sensitivity/spBackgroundDark",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=3, ndr=2)

# dark=b.read("fits",dtype=np.uint16)
# backgroundDark= calBackgroundDark(dark,
#                                     linThreshold,
#                                     mode,
#                                     debug=False,
#                                     logPath=None,
#                                     writeToFile=False)

# gainTable = in_im = fits.open("data/coronalObs-sensitivity/spMasterGain3.000.fits")[0].data.astype(np.float32)

# oldWavecal = fits.open("data/coronalObs-sensitivity/spMasterWavecal-equal.000.fits")[0].data.astype(np.float32)
# # oldWavecal = fits.open("data/coronalObs-sensitivity/spMasterWavecal-plus10.000.fits")[0].data.astype(np.float32)

# badPixels = fits.open("data/coronalObs-sensitivity/spBadPixels.000.fits")[0].data.astype(np.uint8)

# c=cnH2rgRamps("data/coronalObs-sensitivity/spBeamMapping",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=2)
# beamMapping = np.squeeze(c.read("fits",dtype=np.float32))

# fileFormat = "both"
# filePath = "data/coronalObs-sensitivity/"
# sequenceName = "spObserve-ddOneProcessed"
# gratingPosition = 60.96132140227516
# expectedDispersion=0.0174347
# slit="175um"
# centralWavelength = 3934.3
# changeThreshold = 5.

# newWavecal, changeFlag = calWavecal(data,
#                                     centralWavelength,
#                                     gratingPosition,
#                                     slit,
#                                     expectedDispersion,
#                                     backgroundDark,
#                                     gainTable,
#                                     badPixels,
#                                     oldWavecal,
#                                     mode,
#                                     linThreshold,
#                                     changeThreshold,
#                                     debug=False,
#                                     logPath=None,
#                                     simulateChange=False,
#                                     writeToFile=False,
#                                     path=None,
#                                     sequenceName=None,
#                                     fileFormat='fits')


#%%
# fig, ax=plt.subplots()
# ax.plot(a,b,)#a,c,'r')

# from matplotlib.ticker import FormatStrFormatter
# # fig, ax=plt.subplots()
# # im=plt.imshow(data[0,0]-data[0,-1], norm=LogNorm(vmin=0.1, vmax=1000))
# # fig.colorbar(im)
# fig, ax=plt.subplots()
# im=plt.imshow(newWavecal, norm=LogNorm(vmin=0.1, vmax=1000))
# wo = np.arange(0,1023,128)
# wowo = np.arange(0,2047,128)
# vec = np.concatenate((a[wo],d[wo]))
# labels = list(map('{:6.2f}'.format,vec))
# print(labels)
# plt.xticks(wowo,labels,rotation=-90)
# fig.colorbar(im)
# plt.show()