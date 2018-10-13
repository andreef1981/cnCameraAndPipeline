#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
    
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *
from helperFunctions import *
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def calWavecal(data,
               centralWavelength,
               gratingPosition,
               slit,
               backgroundDark,
               gainTable,
               badPixels,
               oldWavecal,
               refSpectrum,
               mode,
               linThreshold,
               debug=False,
               logPath=None,
               simulateChange=False,
               writeToFile=False,
               path=None,
               sequenceName=None,
               fileFormat='fits'):
  #???: should we return the align position in user or raw units
  #      (float32/uint64 for positions)
  #???: should we add a threshold in a property database to decide if align
  #      position has changed?
  #???: Where do we keep the bad pixel mask stored?
  #???: where is delta calculation to base done?
  #???: does the dark need to be a "lamp off thermal dark"?
  """
  Returns wavlength calibration table for the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that for wavelength calibration.
    slit : string
        identifying which slit mask is used 
    gratingPosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the grating stage for each ramp.
    backgroundDark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gainTable: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldWavecal: (2048, 2048) ndarray, float32
        previous wavelength calibration from calibration store
    refSpectrum: (2, x) ndarray, float32
        reference spectrum for given filter from calibration store
    linThreshold : (2,3) ndarray, float32/uint64
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
    gainTable : (2048, 2048) ndarray, float32
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
  
  ################# 6.beam mapping ############################################
  
  #TODO: beam mapping
  
  ################# 7. find spatial profiles###################################
  if len(data.shape)==4:
    gainCorrected = np.average(gainCorrected, axis=0)
  # for slow mode signal is inverted
  if mode is "SLOW":
    result = gainCorrected[0,:,:] - gainCorrected[-1,:,:]
  else:
    result = gainCorrected[-1,:,:] - gainCorrected[0,:,:]
    
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
                                                    [30,2000],
                                                    [30,2000],
                                                    invert=False)
    
    #!!! problem, why are there NaNs in the result (image?)
    
    # cross correlate to reference spectrum
    a = np.squeeze(refSpectrum[1,:])
    b = np.squeeze(leftProfile)
    cross=np.correlate(a,b,mode="same")
    myInd = np.where(cross==np.max(cross))[0][0]
    print("this",myInd)
  newWavecal = result
  # Fifth loop through sequences to find lines from ThAr spectrum
  # TODO: helper function to find emission lines in slit or pinhole spectra
  # TODO: helper function to identify spectral lines from ThAr spectrum
  # TODO: remapping the beams would greatly help with detection of edges an fitting
  
  # Sixth determine if wavecal has changed
  # TODO: what is the threshold for a change
  
  if simulateChange:
    newWavecal = oldWavecal+10.0 # assuming user units simulate change by 10nm
    changeFlag = True
  else:
    # newWavecal = oldWavecal
    changeFlag = False
  
  # for test purposes lets write the files to fits
  if writeToFile:
    if fileFormat == 'fits':
      hdu = fits.PrimaryHDU(newWavecal)
      hdu.writeto(path+sequenceName+'.fits',overwrite=True)
    else:
      newWavecal.tofile(path+sequenceName+'.arr',sep="")
  
  return newWavecal, changeFlag


# reading the data
# cssStyle needs ramps and ndr information
a=cnH2rgRamps("data/coronalObs-sensitivity/spWavecal",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=3, ndr=2)
data = np.squeeze(a.read("fits",dtype=np.uint16))
# reading the background data
b=cnH2rgRamps("data/coronalObs-sensitivity/spMasterBackgroundDark",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=1, ndr=2)
backgroundDark = np.squeeze(b.read("fits",dtype=np.float32))

gainTable = in_im = fits.open("data/coronalObs-sensitivity/spMasterGain3.000.fits")[0].data.astype(np.float32)

oldWavecal = fits.open("data/coronalObs-sensitivity/spWavecal.000.fits")[0].data.astype(np.float32)

badPixels = fits.open("data/coronalObs-sensitivity/spBadPixels.000.fits")[0].data.astype(np.uint8)

c=cnH2rgRamps("data/coronalObs-sensitivity/spBeamMapping",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=1, ndr=2)
beamMapping = np.squeeze(c.read("fits",dtype=np.float32))
#data = np.squeeze(data[0,:,:,:])


######### run this after first part with first part commented
linThreshold = 0
mode = "SLOW"
fileFormat = "both"
filePath = "data/coronalObs-sensitivity/"
sequenceName = "spObserve-ddOneProcessed"
gratingPosition = 60.96132140227516
slit="175um"
centralWavelength = 3934.3
refSpectrum = np.load("data/spectra/SiIX-ThAr-reference-spectrum-175um.npy")

newWavecal, changeFlag = calWavecal(data,
                                    centralWavelength,
                                    gratingPosition,
                                    slit,
                                    backgroundDark,
                                    gainTable,
                                    badPixels,
                                    oldWavecal,
                                    refSpectrum,
                                    mode,
                                    linThreshold,
                                    debug=False,
                                    logPath=None,
                                    simulateChange=False,
                                    writeToFile=False,
                                    path=None,
                                    sequenceName=None,
                                    fileFormat='fits')
#  if fileFormat == "fits":
#    hdu = fits.PrimaryHDU(result)
#    hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(ii)+'.fits',overwrite=True)
#  elif fileFormat == "raw":
#    result.tofile(filePath+sequenceName+'{0:03d}'.format(ii)+'.raw',sep="")
#  elif fileFormat == "both":
#    hdu = fits.PrimaryHDU(result)
#    hdu.writeto(filePath+sequenceName+'.{0:03d}'.format(ii)+'.fits',overwrite=True)
#    result.tofile(filePath+sequenceName+'.{0:03d}'.format(ii)+'.raw',sep="")
#  else:
#    raise ValueError("Specify valid file type.")

fig, ax=plt.subplots()
im=plt.imshow(data[0,0]-data[0,-1], norm=LogNorm(vmin=0.1, vmax=1000))
fig.colorbar(im)
fig, ax=plt.subplots()
im=plt.imshow(newWavecal, norm=LogNorm(vmin=0.1, vmax=1000))
fig.colorbar(im)