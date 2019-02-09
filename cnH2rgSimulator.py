#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 15:43:59 2018

@author: andreef
"""
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from helperFunctions import *

class cnH2rgFrame():
  def __init__(self):
    self.xDim = 2048
    self.yDim = 2048
    self.nrChannels = 32
    self.channelWidth = int(self.xDim/self.nrChannels)
    self.refRows = 4
    self.refCol = 4
    self.adu2e = 4.
    self.frame = np.zeros((self.xDim,self.yDim), dtype=np.uint16)
    self.refDarkCurrent = 0.027 # e-/s
    self.refDarkCurrentTemperature = 37. # k
    # bias Level offsets is random number in [-0.5,0.5]
    self.biasLevelOffsets = np.array([-0.31433745,  0.02008064, -0.40256452,  
                                      0.28055018, -0.01820799,-0.02614954, 
                                      -0.30583004,  0.42315239, -0.06008208, 
                                      0.13673453, -0.09015975,  0.32653709, 
                                      -0.01020625, -0.18490512,  0.19872614,
                                      0.34519787,  0.14878068, -0.02831683, 
                                      0.45009987, -0.4508446 , 0.25490326,  
                                      0.11697454, -0.29125675, -0.42248867,  
                                      0.089145  , 0.23986628, -0.11452058, 
                                      -0.45732304,  0.27197376,  0.26960433,
                                      -0.31100639,  0.44790468])
    
class cnH2rgRamp():
  def __init__(self, mode, ndr, frameTime, frameDealy, biasLevel, biasLevelOffsetScaling, readNoise, arrayTemperature,
               sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, addInstrumentDark=False,
               addInstrumentDarkNoise=False,addThermalDark=False, addThermalDarkNoise=False,
               addFlatQuadraticSignal=False,quadraticCoeff=[2.,1.,0], addFlatQuadraticNoise=False,
               gainVariation=None, 
               spectrum=None, spectrumType=None, spatialProfile=None,
               alignPattern=None, 
               contextImage = None, fileFormat=None):
#  def __init__(self, mode, frameRate, integrationTime, biasChannelDifference,
#               biasChannelVariations, channelReadNoise, channelReadNoiseVariation,
#               gainVariation, deadPixels, negPixels, darkCurrent, ):
    self.mode = mode
    self.ndr = ndr
    self.biasLevel = biasLevel
    self.biasLevelOffsetScaling = biasLevelOffsetScaling
    self.readNoise = readNoise
    self.nrSequences = nrSequences
    self.sequenceName = sequenceName
    self.arrayTemperature = arrayTemperature
    self.frameTime = frameTime
    self.frameDelay = frameDelay
    emptyFrame = cnH2rgFrame()
    if self.mode is "slow":
      self.addingSign = -1.
    else:
      self.addingSign = 1.
    self.spectrum = spectrum  
    # Bias level
    # scale the random channel bias offest [-0.5,0.5] with bias level offset scaling
    # to a percentage of bias level
    bias = cnH2rgFrame() # all zeros
    if addChannelBiases:
      tmp = bias.biasLevelOffsets *self.biasLevelOffsetScaling*self.biasLevel+self.biasLevel
      for i in range(bias.nrChannels):
        bias.frame[:,i*bias.channelWidth:(i+1)*bias.channelWidth] = tmp[i]
    
    # read noise 
    # standard_normal returns a sample with mean 0 and stdev=1
    readNoiseFrame = np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))*self.readNoise
    
    for j in range(self.nrSequences):
      if spectrum is not None:
        if spectrumType is "modulated":
          spectrumPrep = np.tile(self.spectrum[j,:],(2048,1))
        elif spectrumType is "defocused":
          spectrumPrep = np.tile(self.spectrum[j,:],(2048,1))
        else:
          spectrumPrep = np.tile(np.concatenate((self.spectrum,self.spectrum)),(2048,1))
      if contextImage is not None:
        conIm = contextImage[j,:,:]
      if alignPattern is not None:
        # reform into 3D cube
        aSize = alignPattern.shape
        # 'F' to change first index fastest
        alignPatternPrep = np.reshape(alignPattern,[aSize[0],aSize[1],aSize[2]*aSize[3]], order='F')
        if aSize[1] < 2048:
          # sp data needs tiling
          alignPatternPrep = np.concatenate((alignPatternPrep[:,:,j],alignPatternPrep[:,:,j]),axis=1)
        else:
          alignPatternPrep = alignPatternPrep[:,:,j]
      for i in range(self.ndr):
        print(i)
        if addReadnoise:
          readNoiseFrame = np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))*self.readNoise
#          print('mean and std read noise',np.mean(readNoiseFrame),np.std(readNoiseFrame))
        else:
          readNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        #TODO: add gradient in the first frame
        
        # dark current and noise
        # slow mode does pixel by pixel reset, fast has global reset.
        # doubling interval 8K taken from CI dark current measurements.
        instrumentDark = cnH2rgFrame() # all zeros
        instrumentDarkNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        if addInstrumentDark:
          self.instrumentDarkRate = emptyFrame.refDarkCurrent*2.**((self.arrayTemperature-emptyFrame.refDarkCurrentTemperature)/8.)
          tempDarkCurrent = i*self.frameTime*self.instrumentDarkRate
#          print('instrument dark current',tempDarkCurrent)
          # create a dark and in the middle of the instrument dark frame
          bandwidth = 800 # needs to be even
          instrumentDark.frame = instrumentDark.frame+tempDarkCurrent
          instrumentDark.frame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)] = (instrumentDark.frame[:,1024-int(bandwidth/2):1024+
                                        int(bandwidth/2.)]+i*self.frameTime*20.)
          if addInstrumentDarkNoise:
            # use shot noise to simulate dark noise
            instrumentDarkNoiseFrame = np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))*np.sqrt(tempDarkCurrent)
            instrumentDarkNoiseFrame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)] = np.sqrt(instrumentDarkNoiseFrame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)]**2.+
                                                instrumentDark.frame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)])*np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))
            
        thermalDark = cnH2rgFrame() # all zeros
        thermalDarkNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        if addThermalDark:
          # create a dark band in the two beams of the thermal dark frame
          bandwidth = 980 # needs to be even
          thermalScale= 0.1
          thermalDark.frame[:,0:int(bandwidth)] = (thermalDark.frame[:,0:int(bandwidth)]+
                           (i*self.frameTime*quadraticCoeff[1]*thermalScale)+
                          (i*self.frameTime)**2.*quadraticCoeff[0]*thermalScale)
          thermalDark.frame[:,-int(bandwidth)+1:] = (thermalDark.frame[:,-int(bandwidth)+1:]+
                           (i*self.frameTime*quadraticCoeff[1]*thermalScale)+
                          (i*self.frameTime)**2.*quadraticCoeff[0]*thermalScale)
          if addThermalDarkNoise:
            # use shot noise to simulate dark noise
            thermalDarkNoiseFrame = np.sqrt(thermalDark.frame)*np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))
        
        
        #add a flat signal with quadratic increase
        flatSignal = cnH2rgFrame() # all zeros
        flatSignalNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        if addFlatQuadraticSignal:
          # create an illuminated band in the two beams of the  frame
          bandwidth = 990 # needs to be even
          
           # do the gain variation row wise to clearly identify it in the data
          varVector = np.ones(emptyFrame.yDim)
          if gainVariation is not None:
            varVector = gainVariation
          
          flatSignal.frame[:,0:int(bandwidth)] = (flatSignal.frame[:,0:int(bandwidth)]+
                          (i*self.frameTime*quadraticCoeff[1]*varVector)[:,None]+
                          (i*self.frameTime)**2.*quadraticCoeff[0])
          flatSignal.frame[:,-int(bandwidth)+1:] = (flatSignal.frame[:,-int(bandwidth)+1:]+
                          (i*self.frameTime*quadraticCoeff[1]*varVector)[:,None]+
                          (i*self.frameTime)**2.*quadraticCoeff[0])
         
          if spectrum is not None:
#            if len(spectrum.shape) > 2:
#              flatSignal.frame = flatSignal.frame * spectrum[j,:,:]
#            else:
            flatSignal.frame = flatSignal.frame * spectrumPrep
          # TODO: Add slit rotation
          # TODO: Add slit curvature
          # TODO: Add pinhole spectra
          
          if alignPattern is not None:
            flatSignal.frame = flatSignal.frame * alignPatternPrep
          
          if spatialProfile is not None:
            flatSignal.frame = flatSignal.frame * spatialProfile
            
          if contextImage is not None:
            flatSignal.frame = np.multiply(flatSignal.frame,conIm)
            
          if addFlatQuadraticNoise:
            # use normal distribution to simulate dark noise
            flatSignalNoiseFrame = np.sqrt(flatSignal.frame)*np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))
          
        # put it all together and save data
        finalFrame = (bias.frame +
                      readNoiseFrame +
                      self.addingSign*instrumentDark.frame + instrumentDarkNoiseFrame +
                      self.addingSign*thermalDark.frame + thermalDarkNoiseFrame +
                      self.addingSign*flatSignal.frame + flatSignalNoiseFrame)
        if self.mode is "rstrdrd":
          #TODO: conversion to 32 bit might need to happen earlier to get non overflowing calculations.
          finalFrame = np.uint32(finalFrame)
        else:
          finalFrame = np.uint16(finalFrame)
        
        if fileFormat == 'fits':
          #TODO: add the right fits header to the data
          #TODO: change file names to naming convention given in MAN-0008
#          self.writeFits(finalFrame,'data/', self.sequenceName+'-{0:05d}-{1:05d}'.format(j,i)+'.fits',overwrite=True)
          self.writeFits(finalFrame,'data/', self.sequenceName+'.{0:03d}'.format(j*self.ndr+i)+'.fits',overwrite=True)

        elif fileFormat == 'raw':
#          self.writeBinary(finalFrame,'data/', self.sequenceName+'-{0:05d}-{1:05d}'.format(j,i)+'.raw')
          self.writeBinary(finalFrame,'data/', self.sequenceName+'.{0:03d}'.format(j*self.ndr+i)+'.raw')
        elif fileFormat == "both":
          self.writeFits(finalFrame,'data/', self.sequenceName+'.{0:03d}'.format(j*self.ndr+i)+'.fits',overwrite=True)
          self.writeBinary(finalFrame,'data/', self.sequenceName+'.{0:03d}'.format(j*self.ndr+i)+'.raw')
          
  
  
    
    
  def writeBinary(self,data,filepath,filename):
    data.tofile(filepath+filename,sep="")
    
  def writeFits(self,n,filepath,filename,overwrite=False):
    hdu = fits.PrimaryHDU(n)
    #hdul = fits.HDUList([hdu])
    hdu.writeto(filepath+filename,overwrite=overwrite)

# create spatial profiles
# total of 42 groups with an extra center pinhole
spacing = 14
center = np.arange(36,2048,47)
left = np.concatenate((center[0:21],center[22:]))+ spacing
right = np.concatenate((center[0:21],center[22:]))- spacing
y = np.arange(2048)

allY = np.concatenate((left,center,right))
profilePinhole = np.zeros(2048)
profileDiskSlit = np.zeros(2048)
profileDiskSlit[512:1536] = profileDiskSlit[512:1536]+1
profileCoronaSlit = np.zeros(2048)
profileCoronaSlit[5:2042] = profileCoronaSlit[5:2042]+1
kernel = cnGauss(np.arange(100)-50,1,0,20,0)
spatialDisk = np.convolve(profileDiskSlit,kernel,mode="same")/sum(kernel)
spatialCorona = np.convolve(profileCoronaSlit,kernel,mode="same")/sum(kernel)
spatialFwhmPinhole = 3 # pixels
imDisk = np.tile(spatialDisk, (2048,1)).T
imCorona = np.tile(spatialCorona, (2048,1)).T
for ii in range(len(allY)):
  profilePinhole = profilePinhole + cnGauss(y,0.1,allY[ii],spatialFwhmPinhole/2.35,0)
imPinhole = np.tile(profilePinhole, (2048,1)).T  
imSinglePinhole = np.tile(cnGauss(y,1.,1024,10/2.35,0), (2048,1)).T
mode = "slow"

arrayTemperature = 130. 
biasLevelOffsetScaling =  0. # realistic value is 0.001
# bias, readnoise values given in ADU
if mode is "slow":
  biasLevel = 52000.
  readNoise = 5. # 20e-
elif mode is "fast":
  biasLevel = 25000.
  readNoise = 20. # 80e-


heSpectrum = fits.open("data/spectra/he_spectrum_combined-32bitfloat.fits")[0].data.astype(np.float32)
gainVariation = fits.open("data/gain/2048row_gainVariation.fits")[0].data.astype(np.float32)
siSpectrum = np.load("data/spectra/modulated-8-SiIX.npy")
siSpectrumDefocused = np.load("data/spectra/defocused-9-SiIX.npy")
tharSpectrum = np.load("data/spectra/SiIX-ThAr-spectrum.npy")
spAlignPattern = np.load("data/align/spZero.npy")
ciAlignPattern = np.load("data/align/ciIm.npy")
ciModulated = np.load("data/spectra/modulated-8-9pos-contextImager.npy")
frameTime = np.load("frameTime.npy")
print(frameTime)
frameTime = 0.5
frameDelay = 0.
#sequenceName = "/instrumentDark/simInstrumentDark"
#sequenceName = "/backgroundDark/simBackgroundDark"
#sequenceName = "/flatSignal/simFlatSignal"
#sequenceName = "/gain/simGainVariation"
#sequenceName = "/spectra/simSpectrum"

# ########## ci dark
# ndr = 4
# nrSequences = 3
# sequenceName = "/coronalObs-sensitivity/ciBackgroundDark"
# a=cnH2rgRamp(mode, ndr, frameTime, frameDelay, biasLevel, biasLevelOffsetScaling, readNoise,
#               arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
#               addInstrumentDark=True,addInstrumentDarkNoise=False,
#               addThermalDark=True, addThermalDarkNoise=False,
#               addFlatQuadraticSignal=None,quadraticCoeff=[1000,5000.,0], addFlatQuadraticNoise=False,
#               gainVariation=None,
#               spectrum=None, spectrumType=None,
#               spatialProfile=None, alignPattern=None,fileFormat='both') 

########## ci align dark
# ndr = 5
# nrSequences = 3
# sequenceName = "/coronalObs-sensitivity/ciAlign1-backgroundDark"
# a=cnH2rgRamp(mode, ndr, frameTime, frameDelay, biasLevel, biasLevelOffsetScaling, readNoise,
#               arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
#               addInstrumentDark=True,addInstrumentDarkNoise=False,
#               addThermalDark=True, addThermalDarkNoise=False,
#               addFlatQuadraticSignal=None,quadraticCoeff=[1000,5000.,0], addFlatQuadraticNoise=False,
#               gainVariation=None,
#               spectrum=None, spectrumType=None,
#               spatialProfile=None, alignPattern=None,fileFormat='both') 
 
# sequenceName = "/coronalObs-sensitivity/spBackgroundDark"
#sequenceName = "/coronalObs-sensitivity/spGain3"

# ############ ci gain1
# sequenceName = "/coronalObs-sensitivity/ciGain1"
# ndr = 4
# nrSequences = 5
# a=cnH2rgRamp(mode, ndr, frameTime, frameDelay, biasLevel, biasLevelOffsetScaling, readNoise,
#               arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
#               addInstrumentDark=True,addInstrumentDarkNoise=False,
#               addThermalDark=True, addThermalDarkNoise=False,
#               addFlatQuadraticSignal=True,quadraticCoeff=[1000,5000.,0], addFlatQuadraticNoise=False,
#               gainVariation=gainVariation,
#               spectrum=None, spectrumType=None,
#               spatialProfile=None, alignPattern=None,fileFormat='both') 

############ ci gain1
sequenceName = "/coronalObs-sensitivity/ciObserve"
ndr = 4
nrSequences = 72
a=cnH2rgRamp(mode, ndr, frameTime, frameDelay, biasLevel, biasLevelOffsetScaling, readNoise,
              arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
              addInstrumentDark=True,addInstrumentDarkNoise=False,
              addThermalDark=True, addThermalDarkNoise=False,
              addFlatQuadraticSignal=True,quadraticCoeff=[1000,5000.,0], addFlatQuadraticNoise=False,
              gainVariation=gainVariation,
              spectrum=None, spectrumType=None,
              spatialProfile=None, alignPattern=None,
              contextImage=ciModulated, fileFormat='both') 


#sequenceName = "/coronalObs-sensitivity/spObserve"
#sequenceName = "/coronalObs-sensitivity/ciObserve"
# sequenceName = "/coronalObs-sensitivity/spWavecal175um"
# sequenceName = "/coronalObs-sensitivity/spFocus1-175um"
# sequenceName = "/coronalObs-sensitivity/spFocus-background"
# sequenceName = "/coronalObs-sensitivity/spFocus2-GOSpinhole"
# sequenceName = "/coronalObs-sensitivity/spAlign1"
# sequenceName = "/coronalObs-sensitivity/ciAlign1"
#sequenceName = "/generic/observe"

# sequenceName = "/test/test"

# create fixe gain variation
#varVector = np.random.normal(loc=1.0, scale=0.1, size=2048)
#varVector[:,:1024]=varVector[:,:1024]*0.95
#hdu = fits.PrimaryHDU(varVector)
#hdu.writeto('data/gain/2048row_gainVariation.fits',overwrite=True)


# the pattern form is 'raster'
# nrSequences = spAlignPattern.shape[2]*spAlignPattern.shape[3]
# nrSequences = ciAlignPattern.shape[2]*ciAlignPattern.shape[3]
#plt.imshow(siSpectrum)
#TODO: instrument dark noise has issue with non matching array sizes
# a=cnH2rgRamp(mode, ndr, frameTime, frameDelay, biasLevel, biasLevelOffsetScaling, readNoise,
#               arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
#               addInstrumentDark=True,addInstrumentDarkNoise=False,
#               addThermalDark=True, addThermalDarkNoise=False,
#               addFlatQuadraticSignal=None,quadraticCoeff=[1000,5000.,0], addFlatQuadraticNoise=False,
#               gainVariation=None,
#               spectrum=None, spectrumType=None,
#               spatialProfile=None, alignPattern=None,fileFormat='both')    

# Generic data set
#a=cnH2rgRamp(mode, ndr, frameTime/1e9, frameDelay/1e9, biasLevel, biasLevelOffsetScaling, readNoise,
#             arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
#             addInstrumentDark=True,addInstrumentDarkNoise=False,
#             addThermalDark=True, addThermalDarkNoise=False,
#             addFlatQuadraticSignal=True,quadraticCoeff=[100,500.,0], addFlatQuadraticNoise=False,
#             gainVariation=gainVariation,
#             spectrum=siSpectrum, fileFormat='both')    


