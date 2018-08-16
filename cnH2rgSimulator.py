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
               addInstrumentDarkNoise=False,addThermalDark=False, addThermalDarkNoise=False, fileFormat=None):
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
      for i in range(self.ndr):
        print(i)
        if addReadnoise:
          readNoiseFrame = np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))*self.readNoise
          print('mean and std read noise',np.mean(readNoiseFrame),np.std(readNoiseFrame))
        else:
          readNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        
        # dark current and noise
        # slow mode does pixel by pixel reset, fast has global reset.
        # doubling interval 8K taken from CI dark current measurements.
        instrumentDark = cnH2rgFrame() # all zeros
        instrumentDarkNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        if addInstrumentDark:
          self.instrumentDarkRate = emptyFrame.refDarkCurrent*2.**((self.arrayTemperature-emptyFrame.refDarkCurrentTemperature)/8.)
          tempDarkCurrent = i*self.frameTime*self.instrumentDarkRate
          print('instrument dark current',tempDarkCurrent)
          # create a dark and in the middle of the instrument dark frame
          bandwidth = 800 # needs to be even
          instrumentDark.frame = instrumentDark.frame+tempDarkCurrent
          instrumentDark.frame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)] = (instrumentDark.frame[:,1024-int(bandwidth/2):1024+
                                        int(bandwidth/2.)]+i*self.frameTime*20.)
          if addInstrumentDarkNoise:
            # use normal distribution to simulate dark noise
            instrumentDarkNoiseFrame = np.random.standard_normal(size=(emptyFrame.xDim,emptyFrame.yDim))*np.sqrt(tempDarkCurrent)
            instrumentDarkNoiseFrame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)] = np.sqrt(instrumentDarkNoiseFrame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)]**2.+
                                                instrumentDark.frame[:,1024-int(bandwidth/2):1024+int(bandwidth/2.)])
            
        thermalDark = cnH2rgFrame() # all zeros
        thermalDarkNoiseFrame = np.zeros((emptyFrame.xDim,emptyFrame.yDim))
        if addThermalDark:
          # create a dark band in the two beams of the thermal dark frame
          bandwidth = 980 # needs to be even
          
          thermalDark.frame[:,0:int(bandwidth)] = (thermalDark.frame[:,0:int(bandwidth)]+i*self.frameTime*40.)
          thermalDark.frame[:,-int(bandwidth)+1:] = (thermalDark.frame[:,-int(bandwidth)+1:]+i*self.frameTime*40.)
          if addThermalDarkNoise:
            # use normal distribution to simulate dark noise
            thermalDarkNoiseFrame = np.sqrt(thermalDark.frame)
        
        # put it all together and save data
        finalFrame = (bias.frame +
                      readNoiseFrame +
                      self.addingSign*instrumentDark.frame + instrumentDarkNoiseFrame +
                      self.addingSign*thermalDark.frame + thermalDarkNoiseFrame )
        finalFrame = np.uint16(finalFrame)#+readNoiseFrame)
        if fileFormat == 'fits':
          self.writeFits(finalFrame,'data/', self.sequenceName+'-{0:05d}-{1:05d}'.format(j,i)+'.fits',overwrite=True)
        elif fileFormat == 'arr':
          self.writeBinary(finalFrame,'data/', self.sequenceName+'-{0:05d}-{1:05d}'.format(j,i)+'.arr')

          
  
  
    
    
  def writeBinary(self,data,filepath,filename):
    data.tofile(filepath+filename,sep="")
    
  def writeFits(self,n,filepath,filename,overwrite=False):
    hdu = fits.PrimaryHDU(n)
    #hdul = fits.HDUList([hdu])
    hdu.writeto(filepath+filename,overwrite=overwrite)

#%%   
class cnPb1Calculator():
  def __init__(self,mode,ndr,coadd, rows, yCoord, pixelClkOscPulses, 
               resetLength, betweenFrameDelay, rstrdrdIntegration, 
               interspersedDelay, resetSettlePause, exposureTime=None):
    # internally all times are stored in ns
    # exposure time has to be entered in seconds
    self.mode = mode
    self.clockPulse = 20 # 20ns clock pulses
    self.ndr = ndr
    self.coadd = coadd
    self.channelWidth = 64     # in pixel
    self.rows = rows
    self.yCoord = yCoord
    assert (self.rows <= 2048 - self.yCoord), "Number of rows is too big for chosen yCoordinate."
    self.pixelClkOscPulses = pixelClkOscPulses
    self.resetLength = resetLength
    self.betweenFrameDelay = betweenFrameDelay     # number of 20 ns clock pulses
    self.rstrdrdIntegration = rstrdrdIntegration   # number of 20 ns clock pulses
    self.interspersedDelay = interspersedDelay     # 1us pulses
    self.resetSettlePause = resetSettlePause       # number of 20 ns clock pulses
    
    # slow mode calculations ##################################################
    if mode is 'slow':
      tFpgaReset = 2000      # time in ns
      tFsynch = 620          # time in ns
      tVclk = 420            # time in ns
      tFirstDwell = None     # time in ns
      tLastLow = None        # time in ns
      self.lineTime = (tVclk+(self.channelWidth+2)*(self.pixelClkOscPulses+1)*self.clockPulse)
      # tVclk*yCoord is time to address (move to) first row of subarray
      self.frameTime = (tFsynch+self.yCoord*tVclk+self.rows*self.lineTime)
      if exposureTime is not None:
        self.exposureTime = exposureTime*1e9  # time is in ms
        self.betweenFrameDelayTime = (self.exposureTime/(self.ndr-1) - self.interspersedDelay*2000. - self.frameTime)
        self.betweenFrameDelay = self.betweenFrameDelayTime/self.clockPulse
        assert (self.betweenFrameDelayTime >= 0.), "Timing ERROR: Exposur time too short."
      else:
        self.betweenFrameDelayTime = self.betweenFrameDelay*self.clockPulse
        self.exposureTime = (self.ndr-1)*(self.betweenFrameDelayTime+ self.interspersedDelay*2000. + self.frameTime)
      #!!! check the in between multiplicator should be ndr-1
      self.totalbetweenFrameDelayTime = self.betweenFrameDelayTime+ self.interspersedDelay*2000.
      self.rampTime = (tFpgaReset+self.coadd*(self.ndr+1)*(self.frameTime)+
                       self.coadd*((self.ndr-1)*self.betweenFrameDelay*self.clockPulse+self.interspersedDelay*2000))
      
      self.firstPixelRead = tFpgaReset + tFsynch + tVclk*self.yCoord + tVclk
    
    # fast mode calculations ##################################################
    elif mode is 'fast':
      tFpgaReset = 2000      # time in ns
      tFsynch = 320          # time in ns
      tVclk = 220            # time in ns
      tFirstDwell = 2040     # time in ns
      tLastLow = 120         # time in ns
      self.lineTime = (tVclk+tFirstDwell+tLastLow+(self.channelWidth+1.5)*(self.pixelClkOscPulses+1)*self.clockPulse)
      self.frameTime = (tFsynch+self.yCoord*tVclk+(self.rows+1)*self.lineTime)
      if exposureTime is not None:
        self.exposureTime = exposureTime*1e9  # time is in ms
        self.betweenFrameDelayTime = (self.exposureTime/(self.ndr-1) - self.interspersedDelay*1000. - self.frameTime)
        self.betweenFrameDelay = self.betweenFrameDelayTime/self.clockPulse
        assert (self.betweenFrameDelayTime >= 0.), "Timing ERROR: Exposur time too short."
      else:
        self.betweenFrameDelayTime = self.betweenFrameDelay*self.clockPulse
        self.exposureTime = (self.ndr-1)*(self.betweenFrameDelayTime+ self.interspersedDelay*1000. + self.frameTime)
      # !!! check the in between multiplicator should be ndr-1 
      self.totalbetweenFrameDelayTime = self.betweenFrameDelayTime+ self.interspersedDelay*2000.
      self.rampTime = (tFpgaReset + self.coadd*(self.ndr*self.frameTime +
                         (self.resetLength+self.resetSettlePause)*self.clockPulse) +
                         self.coadd*((self.ndr-1)*self.betweenFrameDelay*self.clockPulse + self.interspersedDelay*1000))
      self.firstPixelRead = -99.
    
    # line mode calculations ##################################################
    elif mode is 'rstrdrd':
      tFpgaReset = 2000      # time in ns
      tFsynch = 320          # time in ns
      tVclk = 220            # time in ns
      tFirstDwell = 2040     # time in ns
      tLastLow = 120         # time in ns
      self.firstPixelRead = -99.
    
    
    self.frameRate = 1/(self.frameTime/1e9)
    
    
#%%  

mode = "slow"
ndr = 10
coadd =1

nrSequences = 5
frameRate = 12.
arrayTemperature = 130. 
rows = 2048
ycoord = 0
t = None

biasLevelOffsetScaling =  0. # realistic value is 0.001
# bias, readnoise values given in ADU
if mode is "slow":
  biasLevel = 52000.
  readNoise = 5. # 20e-
elif mode is "fast":
  biasLevel = 25000.
  readNoise = 20. # 80e-


if mode is 'slow':
  pixelClkOscPulses = 150.  # number of 20 ns clock pulses
  resetLength = None        # number of 20 ns clock pulses
  betweenFrameDelay = 0     # number of 20 ns clock pulses
  rstrdrdIntegration = None # number of 20 ns clock pulses
  interspersedDelay = 1000  # 1us pulses
  resetSettlePause = 0      # number of 20 ns clock pulses
    
elif mode is 'fast':
  pixelClkOscPulses = 28.   # number of 20 ns clock pulses
  resetLength = 1000        # number of 20 ns clock pulses
  betweenFrameDelay = 0     # number of 20 ns clock pulses
  rstrdrdIntegration = None # number of 20 ns clock pulses
  interspersedDelay = 1000  # 1us pulses
  resetSettlePause = 0      # number of 20 ns clock pulses
  
elif mode is 'rstrdrd':
  pixelClkOscPulses = 28.   # number of 20 ns clock pulses
  resetLength = 1000        # number of 20 ns clock pulses
  betweenFrameDelay = 0     # number of 20 ns clock pulses
  rstrdrdIntegration = 0    # number of 20 ns clock pulses
  interspersedDelay = 1000  # 1us pulses
  resetSettlePause = 0      # number of 20 ns clock pulses   
  coadd = None   

# specifying exposureTime overwrites betweenFrameDelay
# provide exposureTime in [seconds]
c = cnPb1Calculator(mode, ndr, coadd, rows, ycoord, pixelClkOscPulses, 
               resetLength, betweenFrameDelay, rstrdrdIntegration, 
               interspersedDelay, resetSettlePause, exposureTime=t)
frameTime = c.frameTime
frameDelay = c.totalbetweenFrameDelayTime
print('line time [ms]',c.lineTime/1e6)  
print('frame time [ms]',c.frameTime/1e6)
print('frame rate [Hz]',c.frameRate) 
print('ramp time [ms]',c.rampTime/1e6)
print('first pixel read after [us]',c.firstPixelRead*1000.)
print('in between frame dealy [ms]',c.betweenFrameDelayTime/1e6)
print('exposure time [ms]',c.exposureTime/1e6)

#sequenceName = "/instrumentDark/simInstrumentDark"
sequenceName = "/backgroundDark/simBackgroundDark"
a=cnH2rgRamp(mode, ndr, frameTime/1e9, frameDelay/1e9, biasLevel, biasLevelOffsetScaling, readNoise,
             arrayTemperature, sequenceName, nrSequences, addChannelBiases=True, addReadnoise=True, 
             addInstrumentDark=True,addInstrumentDarkNoise=True,
              addThermalDark=True, addThermalDarkNoise=True,fileFormat='arr')    

# flux tester
#f1= 100
#f2 = 200
#a = -5e-4
#b = 4  
#c = 0     
#t = np.arange(0,20.1, 0.1)
#t1 = f1*t
#t2 = f2*t
##adu = a*t**2 + b*t +c
#adu1 = a*t1**2. + b*t1 + c
#adu2 = a*t2**2. + b*t2 + c
#fig, ax=plt.subplots()
#ax.plot(t,adu1, t, adu2, 'r')
##ax.plot(t,adu)
#plt.show()
      