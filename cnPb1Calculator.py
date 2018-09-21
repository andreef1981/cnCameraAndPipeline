#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:04:21 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision History
----------------
  - 20 Sep 2018 AFe Implement final changes from Erics spread sheet
"""

import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

class cnPb1Calculator():
  def __init__(self,mode,ndr,coadd, rows, yCoord, pixelClkOscPulses, 
               betweenFrameDelay, rstrdrdIntegration, 
               interspersedDelay, exposureTime=None):
    """
    if exposure time is specified then the inter frame delay will be overwritten
    if ndr==None then maximizing but requires exposure time
    """
    
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
    #self.resetLength = resetLength
    self.betweenFrameDelay = betweenFrameDelay     # number of 20 ns clock pulses
    self.rstrdrdIntegration = rstrdrdIntegration   # number of 20 ns clock pulses
    self.interspersedDelay = interspersedDelay     # 1us pulses
#    self.resetSettlePause = resetSettlePause       # number of 20 ns clock pulses
    
    # slow mode calculations ##################################################
    if mode is 'slow':
      self.tFpgaReset = 2000      # time in ns
      self.tFsynch = 620          # time in ns
      self.tVclk = 420            # time in ns
      self.tFirstDwell = None     # time in ns
      self.tLastLow = None        # time in ns
      self.tFullFrameReset = None # time in ns
      # line time in ns
      self.lineTime = (self.tVclk+(self.channelWidth+2)*(self.pixelClkOscPulses+1)*self.clockPulse)
      # tVclk*yCoord is time to address (move to) first row of subarray
      self.frameTime = (self.tFsynch+self.yCoord*self.tVclk+self.rows*self.lineTime)
      if exposureTime is not None:
        self.exposureTime = exposureTime*1e9  # all times are in ns
        assert (self.exposureTime >= self.frameTime), "Timing ERROR: Exposur time too cannot be shorter than frame time."
        if self.ndr is None:
          # maximize number of NDRs
          self.ndr = np.floor(self.exposureTime/self.frameTime)+1
        # calculate what delay is needed in between reads
        self.betweenFrameDelayTime = (self.exposureTime-(self.ndr-1)*self.frameTime)/(self.ndr-1)
        assert (self.betweenFrameDelayTime >= 0.), "Timing ERROR: Exposur time too short for given #NDR."
        # betweenFrameDelay should be an integer to match clocking
        # force rounding to integer
        self.betweenFrameDelay = np.round(self.betweenFrameDelayTime/self.clockPulse)
        self.betweenFrameDelayTime = self.betweenFrameDelay*self.clockPulse
        # adjust exposure time accordingly
        self.requestedExposureTime = self.exposureTime
        self.exposureTime = ((self.ndr-1)*self.betweenFrameDelayTime) + (self.ndr-1)*self.frameTime
        if not self.exposureTime==self.requestedExposureTime:
          print("Warning: Exposure time had to be adjusted to meet 20ns clock cycles")
      else:
        assert (self.ndr is not None), "If ndr==None, a exposure time has to be specified to maximize #NDRs"
        self.betweenFrameDelayTime = self.betweenFrameDelay*self.clockPulse
        self.exposureTime = ((self.ndr-1)*self.betweenFrameDelayTime) + (self.ndr-1)*self.frameTime
      
      # calculate time for a ramp
      self.rampTime = (self.tFpgaReset+self.coadd*(self.ndr+1)*(self.frameTime)+
                       self.coadd*((self.ndr-1)*self.betweenFrameDelayTime+
                                   self.interspersedDelay*2000))
      
      self.firstFrameStart = self.tFpgaReset + self.frameTime
      self.firstFrameEnd = self.tFpgaReset + 2*self.frameTime
      self.firstPixelRead = self.tFpgaReset + self.frameTime + self.tFsynch + self.tVclk*self.yCoord + self.tVclk + (self.pixelClkOscPulses+1)*self.clockPulse
      self.lastPixelRead = self.tFpgaReset + 2*self.frameTime - 2*(self.pixelClkOscPulses+1)*self.clockPulse
      # frame time includes Fsynch, Vclk
      self.timeVectorStart = (np.arange(self.ndr)*(self.frameTime+self.betweenFrameDelayTime) +
                              self.firstPixelRead) 
      self.timeVectorEnd = (np.arange(self.ndr)*(self.frameTime+self.betweenFrameDelayTime)+
                            self.lastPixelRead)
    
    # fast mode calculations ##################################################
    elif mode is 'fast':
      self.tFpgaReset = 2000         # time in ns
      self.tFsynch = 180             # time in ns
      self.tVclk = 120               # time in ns
      self.tFirstDwell = 4000        # time in ns
      self.tLastLow = 60             # time in ns
      self.tFullFrameReset = 2090000 # time in ns
      
      # line time in ns
      self.lineTime = (self.tVclk+self.tFirstDwell+self.tLastLow+(self.channelWidth+1.5)*(self.pixelClkOscPulses+1)*self.clockPulse)
      # tVclk*yCoord is time to address (move to) first row of subarray
      self.frameTime = (self.tFsynch+self.yCoord*self.tVclk+(self.rows+1)*self.lineTime)
      if exposureTime is not None:
        self.exposureTime = exposureTime*1e9  # all times are in ns
        assert (self.exposureTime >= self.frameTime), "Timing ERROR: Exposur time too cannot be shorter than frame time."
        if self.ndr is None:
          # maximize number of NDRs fast mode transmits unusable first frame
          self.ndr = np.floor(self.exposureTime/self.frameTime)+2
        # calculate what delay is needed in between reads
        self.betweenFrameDelayTime = (self.exposureTime-(self.ndr-2)*self.frameTime)/(self.ndr-2)
        assert (self.betweenFrameDelayTime >= 0.), "Timing ERROR: Exposur time too short for given #NDR."
        # betweenFrameDelay should be an integer to match clocking
        # force rounding to integer
        self.betweenFrameDelay = np.round(self.betweenFrameDelayTime/self.clockPulse)
        self.betweenFrameDelayTime = self.betweenFrameDelay*self.clockPulse
        # adjust exposure time accordingly
        self.requestedExposureTime = self.exposureTime
        self.exposureTime = ((self.ndr-2)*self.betweenFrameDelayTime) + (self.ndr-2)*self.frameTime
        if not self.exposureTime==self.requestedExposureTime:
          print("Warning: Exposure time had to be adjusted to meet 20ns clock cycles")
      else:
        assert (self.ndr is not None), "If ndr==None, a exposure time has to be specified to maximize #NDRs"
        self.betweenFrameDelayTime = self.betweenFrameDelay*self.clockPulse
#        self.requestedExposureTime = exposureTime
        self.exposureTime = ((self.ndr-2)*self.betweenFrameDelayTime) + (self.ndr-2)*self.frameTime
      
      # calculate time for a ramp
      self.rampTime = self.tFullFrameReset+(self.tFpgaReset+self.coadd*(self.ndr)*(self.frameTime)+
                       self.coadd*((self.ndr-1)*self.betweenFrameDelayTime+self.interspersedDelay*1000))
      
      self.firstFrameStart = self.tFullFrameReset + self.tFpgaReset 
      self.firstFrameEnd = self.tFullFrameReset + self.tFpgaReset + self.frameTime
      self.firstPixelRead = self.tFullFrameReset + self.tFpgaReset + self.tFsynch + self.tVclk*self.yCoord + self.tFirstDwell + self.tVclk + (self.pixelClkOscPulses+1)*self.clockPulse
      self.lastPixelRead = self.tFullFrameReset + self.tFpgaReset + self.frameTime - self.tLastLow - 1.5*(self.pixelClkOscPulses+1)*self.clockPulse
      # frame time includes Fsynch, Vclk
      self.timeVectorStart = (np.arange(self.ndr)*(self.frameTime+self.betweenFrameDelayTime) +
                              self.firstPixelRead) 
      self.timeVectorEnd = (np.arange(self.ndr)*(self.frameTime+self.betweenFrameDelayTime)+
                            self.lastPixelRead)
#      self.timeVectorEnd[-1] = self.timeVectorEnd[-1] - self.betweenFrameDelayTime
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

mode = "fast"
ndr = None
coadd =1 
rows = 2048
ycoord = 0
t = 0.1 # exposure time given in seconds

if mode is 'slow':
  # pcopf value
  pixelClkOscPulses = 184.  # number of 20 ns clock pulses
  #resetLength = None        # number of 20 ns clock pulses
  # -t value
  betweenFrameDelay = 0     # number of 20 ns clock pulses
  rstrdrdIntegration = None # number of 20 ns clock pulses
  interspersedDelay = 0     # 2us pulses
#  resetSettlePause = 0      # number of 20 ns clock pulses
    
elif mode is 'fast':
  # pcopf value
  pixelClkOscPulses = 28.   # number of 20 ns clock pulses
  #resetLength = None        # number of 20 ns clock pulses
  # -t value
  betweenFrameDelay = 0     # number of 20 ns clock pulses
  rstrdrdIntegration = None # number of 20 ns clock pulses
  interspersedDelay = 1000  # 1us pulses
#  resetSettlePause = 0      # number of 20 ns clock pulses
  
elif mode is 'rstrdrd':
  # pcopf value
  pixelClkOscPulses = 28.   # number of 20 ns clock pulses
  #resetLength = 1000        # number of 20 ns clock pulses
  betweenFrameDelay = None     # number of 20 ns clock pulses
  # -t value
  rstrdrdIntegration = 0    # number of 20 ns clock pulses
  interspersedDelay = 1000  # 1us pulses
#  resetSettlePause = 0      # number of 20 ns clock pulses
  # make sure no coadding is implemented   
  coadd = None   

# specifying exposureTime overwrites betweenFrameDelay
# provide exposureTime in [seconds]
c = cnPb1Calculator(mode, ndr, coadd, rows, ycoord, pixelClkOscPulses, 
               betweenFrameDelay, rstrdrdIntegration, 
               interspersedDelay, exposureTime=t)
frameTime = c.frameTime
#frameDelay = c.totalbetweenFrameDelayTime
print('line time [ms]',c.lineTime/1e6)  
print('frame time [ms]',c.frameTime/1e6)
print('frame rate [Hz]',c.frameRate) 
print('ramp time [ms]',c.rampTime/1e6)
print('#NDRs',c.ndr) 
print('first frame start [ms]',c.firstFrameStart/1e6)
print('first pixel read after [ms]',c.firstPixelRead/1e6)
print('last pixel read after [ms]',c.lastPixelRead/1e6)
print('first frame end [ms]',c.firstFrameEnd/1e6)

print('in between frame dealy [ms]',c.betweenFrameDelayTime/1e6)
try:
  print('requested exposure time [ms]',c.requestedExposureTime/1e6)
except:
  print("No specific exposure time was requested")
print('effective exposure time [ms]',c.exposureTime/1e6)
print(c.timeVectorStart/1e6)
print(c.timeVectorEnd/1e6)
print((c.timeVectorEnd-c.timeVectorStart)/1e6)