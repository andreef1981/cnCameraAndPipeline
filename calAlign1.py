#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:46:37 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
  06 Feb 2016: - added calculation to figure out alignment, step size and angle
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *
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
  #!!!: need to add slit information
  #!!!: need to add stage scale information or becky to use it to convert
  #     arcsec information into mu for stage
  #???: How do we ensure the center position in mu is updated in porperty database 
  #     from calibration store when IC is brought up?
  #???: What is the propper scale in the property data base arcsec/mu
  
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
    pixelScaleFM1A : float32
        Number describing the spatial coverage of a single pixel.
    pixelScaleFM1A : float32
        Number describing the spatial coverage of a single pixel.
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
      - GOS has two pinholse 0.2 and 1 arcsec
      - Assumes a raster scanning pattern (not optimized)
    
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
  
  
  ################# 5. flat fielding ##########################################
  gainCorrected = backgroundSubtracted*gainTable
  
  # if len(data.shape)==4:
  #   gainCorrected = np.average(gainCorrected, axis=0)
  # for slow mode signal is inverted
  if mode == "SLOW":
    result = gainCorrected[:,0,:,:] - gainCorrected[:,-1,:,:]
  else:
    result = gainCorrected[:,-1,:,:] - gainCorrected[:,0,:,:]
  
  ################ 6. finding spatial peaks #################################
  xloc = np.zeros([stagePositionFM1A.size])
  yloc = np.zeros([stagePositionFM1B.size])
  
  fig, ax=plt.subplots(num=2)
  ax.imshow((result[0,:,:]))  
  for i in range(len(stagePositionFM1A)):
    if camera == "SP":
      # only use left side of spectrum
      im = result[i,:,:1024]
      ysize = 2048
      xsize = 1024
    else:
      im = result[i,:,:]
      ysize = 2048
      xsize = 2048
    # find peaks along x axis (0)
    xInd, xWidths, xProminences, xHeights, xProfile = cnFind2DSpatialPeak(im,
                                                                          0,
                                                                          [6,30],
                                                                          [0,ysize-1],
                                                                          None,
                                                                          [20.,65000.],
                                                                          [10,65000.],
                                                                          invert=False)
    # find peaks along y axis (1)
    yInd, yWidths, yProminences, yHeights, yProfile = cnFind2DSpatialPeak(im,
                                                                          1,
                                                                          [6,30],
                                                                          [0,xsize-1],
                                                                          None,
                                                                          [20.,65000.],
                                                                          [10,65000.],
                                                                          invert=False)
    
    # raise warning if one is bigger than size==1
    assert(xInd.size==1),"More than one peak in x direction found"
    assert(yInd.size==1),"More than one peak in y direction found"
    xloc[i] = xInd
    yloc[i] = yInd
  
  ################## 7. Find align position, step size and rotation #################
  # for a raster scanning the y axis is the slow one
  givenY = np.unique(stagePositionFM1B)  
  givenX= np.unique(stagePositionFM1A) 
  # reshap to this 
  xloc = np.reshape(xloc,[len(givenY),int(len(stagePositionFM1B)/len(givenY))])
  yloc = np.reshape(yloc,[len(givenY),int(len(stagePositionFM1B)/len(givenY))])
  
  xCenterPixel = np.zeros(givenY.size)
  xStep = np.zeros(givenY.size)
  yCenterPixel = np.zeros(givenX.size)
  yStep = np.zeros(givenX.size)
  # do linear fit to find x scale and center at each y
  # p[0] is in arcsec/pixel
  for j in range(givenY.size):
    p = np.polyfit(xloc[j,:],givenX,1)
    xCenterPixel[j] = -1*p[1]/p[0]
    xStep[j] = 1/p[0] # pixel/arcsec
  for i in range(givenX.size):
    p = np.polyfit(yloc[:,i],givenY,1)
    yCenterPixel[i] = -1*p[1]/p[0]
    yStep[i] = 1/p[0] # pixel/arcsec
  
  # need to find the arcsec value that puts spot on center pixel
  xCenterOfY = np.zeros(givenY.size)
  xSlopes = np.zeros(givenY.size)
  for j in range(givenY.size):
    p = np.polyfit(xCenterPixel,yloc[:,j],1)
    xCenterOfY[j] = (ysize/2-p[1])/p[0]
    xSlopes[j] = p[0]
  xShiftArcsec = (np.mean(xCenterOfY)-xsize/2)/np.mean(xStep)
  xAngle = np.mean(np.rad2deg(np.arctan(xSlopes)))

  yCenterOfX = np.zeros(givenY.size)
  ySlopes = np.zeros(givenY.size)
  for j in range(givenY.size):
    p = np.polyfit(xloc[j,:],yCenterPixel,1)
    yCenterOfX[j] = (xsize/2*p[0]) + p[1]
    ySlopes[j] = p[0]
  yShiftArcsec = (np.mean(yCenterOfX)-ysize/2)/np.mean(yStep)
  yAngle = np.mean(np.rad2deg(np.arctan(ySlopes)))

  newAlignFM1A = xShiftArcsec
  newAlignFM1B = yShiftArcsec
  newStepSizeFM1A = np.mean(xStep)
  newStepSizeFM1B = np.mean(yStep)
  newRotationFM1A = xAngle
  newRotationFM1B = yAngle
  
  
  if simulateChange[0,0]:
    newAlignFM1A = oldAlignFM1A+0.05 # assuming user units simulate change by 0.05 degrees
    changeFlag[0,0] = True
  
  if simulateChange[1,0]:
    newAlignFM1B = oldAlignFM1B-0.05 # assuming user units simulate change by 0.05 degrees
    changeFlag[1,0] = True
  
  if simulateChange[0,1]:
    newStepSizeFM1A = oldStepSizeFM1A+3.0 # simulate change by 3 mu per pixel
    changeFlag[0,1] = True
  
  if simulateChange[1,1]:
    newStepSizeFM1B = oldStepSizeFM1B-3.0 # simulate change by 3 mu per pixel
    changeFlag[1,1] = True
  
  if simulateChange[0,2]:
    newRotationFM1A = oldRotationFM1A+0.8 # simulate change by 0.8 degrees
    changeFlag[0,2] = True
  
  if simulateChange[1,2]:
    newRotationFM1B = oldRotationFM1B-0.8 # simulate change by 0.8 degrees
    changeFlag[1,2] = True
  
  
  ################# 8. Check if values have changed ############################
  if np.abs(newAlignFM1A - oldAlignFM1A) > changeThreshold[0,0]:
    changeFlag[0,0] = True
  else:
    changeFlag[0,0] = False 
  
  if np.abs(newAlignFM1B - oldAlignFM1B) > changeThreshold[1,0]:
    changeFlag[1,0] = True
  else:
    changeFlag[1,0] = False 
  
  if np.abs(newStepSizeFM1A - oldStepSizeFM1A) > changeThreshold[0,1]:
    changeFlag[0,1] = True
  else:
    changeFlag[0,1] = False 
    
  if np.abs(newStepSizeFM1B - oldStepSizeFM1B) > changeThreshold[1,1]:
    changeFlag[1,1] = True
  else:
    changeFlag[1,1] = False 
    
  if np.abs(newRotationFM1A - oldRotationFM1A) > changeThreshold[0,2]:
    changeFlag[0,2] = True
  else:
    changeFlag[0,2] = False 
    
  if np.abs(newRotationFM1B - oldRotationFM1B) > changeThreshold[1,2]:
    changeFlag[1,2] = True
  else:
    changeFlag[1,2] = False 
  
  
  return newAlignFM1A, newAlignFM1B, newStepSizeFM1A, newStepSizeFM1B,\
         newRotationFM1A, newRotationFM1B, changeFlag

#%%
# spectrograph 0 order
xscale = 8.3333 #pixel/arcsec
yscale = 8.3333 #pixel/arcsec
givenX = np.arange(-19, 20, 5)
givenY = np.arange(-110, 120, 50)

# raster type scanning
stagePositionFM1A = np.zeros([givenX.size*givenY.size])
stagePositionFM1B = np.zeros([givenX.size*givenY.size])
for i in range(len(givenY)):
  for j in range(len(givenX)):

    stagePositionFM1A[i*len(givenX)+j] = givenX[j]
    stagePositionFM1B[i*len(givenX)+j] = givenY[i]
#%%    
# reading the data
# cssStyle needs ramps and ndr information


a=cnH2rgRamps("data/coronalObs-sensitivity/spAlign1.*",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=40, ndr=5)
data = np.squeeze(a.read("fits",dtype=np.uint16))
#%%
# reading the background data
b=cnH2rgRamps("data/coronalObs-sensitivity/spAlign1-masterBackground",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=1, ndr=5)
backgroundDark=np.squeeze(b.read("fits",dtype=np.float32))

#%%
gainTable = fits.open("data/coronalObs-sensitivity/spMasterGain3.000.fits")[0].data.astype(np.float32)
badPixels = fits.open("data/coronalObs-sensitivity/spBadPixels.000.fits")[0].data.astype(np.uint8)

c=cnH2rgRamps("data/coronalObs-sensitivity/spBeamMapping",
              "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
              ramps=1, ndr=2)
beamMapping = np.squeeze(c.read("fits",dtype=np.float32))

#%%
linThreshold = 0
mode = "SLOW"
camera = "SP"
# camera = "CI"
oldAlignFM1A = 0.  # arcsec
oldAlignFM1B = 0.  # arcsec
oldStepSizeFM1A = 10. # pixel/arcsec
oldStepSizeFM1B = 10.
oldRotationFM1A = 0. # degrees
oldRotationFM1B = 0.
changeFlag = np.zeros([2,3])
changeThreshold = np.zeros([2,3])

newAlignFM1A, newAlignFM1B, newStepSizeFM1A,\
newStepSizeFM1B, newRotationFM1A,\
newRotationFM1B, changeFlag = calAlign1(data[:,:,:,:],
                                        stagePositionFM1A[:],
                                        stagePositionFM1B[:],
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
                                        changeThreshold=changeThreshold,
                                        simulateChange=np.zeros((2,3),dtype=bool),
                                        debug=False,
                                        logPath=None,
                                        writeToFile=False,
                                        filePath=None,
                                        sequenceName=None,
                                        fileFormat='fits')



# #%%
# # context imager
# xscale = 20 #pixel/arcsec
# yscale = 20 #pixel/arcsec
# givenX = np.arange(-45, 45, 10)
# givenY = np.arange(-45, 45, 20)

# # raster type scanning
# stagePositionFM1A = np.zeros([givenX.size*givenY.size])
# stagePositionFM1B = np.zeros([givenX.size*givenY.size])
# for i in range(len(givenY)):
#   for j in range(len(givenX)):

#     stagePositionFM1A[i*len(givenX)+j] = givenX[j]
#     stagePositionFM1B[i*len(givenX)+j] = givenY[i]
# #%%    
# # reading the data
# # cssStyle needs ramps and ndr information
# a=cnH2rgRamps("data/coronalObs-sensitivity/ciAlign1.*",
#             "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#             ramps=45, ndr=5)
# data = np.squeeze(a.read("fits",dtype=np.uint16))
# #%%
# # reading the background data
# b=cnH2rgRamps("data/coronalObs-sensitivity/ciAlign1-masterBackground",
#               "fits",readMode="SLOW",subArray=None,verbose=True, cssStyle=True,
#               ramps=1, ndr=5)
# backgroundDark=np.squeeze(b.read("fits",dtype=np.float32))

# #%%
# gainTable = fits.open("data/coronalObs-sensitivity/ciMasterGain1.000.fits")[0].data.astype(np.float32)

# badPixels = fits.open("data/coronalObs-sensitivity/ciBadPixels.000.fits")[0].data.astype(np.uint8)

# beamMapping = None

# #%%
# linThreshold = 0
# mode = "SLOW"
# camera = "CI"
# oldAlignFM1A = 0.  # arcsec
# oldAlignFM1B = 0.  # arcsec
# oldStepSizeFM1A = 10. # pixel/arcsec
# oldStepSizeFM1B = 10.
# oldRotationFM1A = 0. # degrees
# oldRotationFM1B = 0.
# changeFlag = np.zeros([2,3])
# changeThreshold = np.zeros([2,3])

# newAlignFM1A, newAlignFM1B, newStepSizeFM1A,\
# newStepSizeFM1B, newRotationFM1A,\
# newRotationFM1B, changeFlag = calAlign1(data[:,:,:,:],
#                                         stagePositionFM1A[:],
#                                         stagePositionFM1B[:],
#                                         camera,
#                                         backgroundDark,
#                                         gainTable,
#                                         badPixels,
#                                         oldAlignFM1A,
#                                         oldAlignFM1B,
#                                         oldStepSizeFM1A,
#                                         oldStepSizeFM1B,
#                                         oldRotationFM1A,
#                                         oldRotationFM1B,
#                                         linThreshold,
#                                         mode,
#                                         changeThreshold=changeThreshold,
#                                         simulateChange=np.zeros((2,3),dtype=bool),
#                                         debug=False,
#                                         logPath=None,
#                                         writeToFile=False,
#                                         filePath=None,
#                                         sequenceName=None,
#                                         fileFormat='fits')

