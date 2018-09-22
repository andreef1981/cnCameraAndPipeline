#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)
helper library
"""

import numpy as np  # Fast numeric arrays/matrices (FITS files imported as numpy matrices)

def inversematrices(data):
    # From Kathleen Tatem
    N = len(data[0,:,0,0])
    dim = (N, 6)
    
    #initial inverse array that is N by 6
    Cinv = np.zeros(dim)
    
    times = np.arange(N)
    for t in times:
        #create inverse matrix depending on the length of the ramp 
        t0 = np.ones(t+1)
        t1 = np.arange(t+1)
        t2 = t1**2
        
        #the Cij matrix is a symmetric matrix, so there are only 6 elements we need to calculate
        c00 = t0.sum(axis=0)
        c01 = t1.sum(axis=0)
        c02 = t2.sum(axis=0)
        c11 = np.dot(t1, t1)
        c12 = np.dot(t1, t2)
        c22 = np.dot(t2, t2)
        
        #fill the Cij matrix
        Cij = np.zeros((3, 3))
        Cij[0, 0] = c00
        Cij[0, 1] = c01
        Cij[1, 0] = c01
        Cij[0, 2] = c02
        Cij[2, 0] = c02
        Cij[1, 1] = c11
        Cij[1, 2] = c12
        Cij[2, 1] = c12
        Cij[2, 2] = c22        
        
        #take the inverse of the Cij matrix
        #creates a 3x3 array that is the inverse of Cij (when t = 2, 1, or 0, the matrix is not invertible)
        c = np.zeros(6)
        if t > 2:
            invCji = np.linalg.inv(Cij)
            c[0] = invCji[0,0]
            c[1] = invCji[1,1]
            c[2] = invCji[2,2]
            c[3] = invCji[0,1]
            c[4] = invCji[0,2]
            c[5] = invCji[1,2]
        
            Cinv[t] = c
                
    return Cinv    

def matrixQuadfit(data,threshold=66000, mode="SLOW",ignoreRef=False):  
    # From Kathleen Tatem
    #don't use the reference pixels in calculations
    
    #TODO: figure out how to treat data with only two frames up the ramp
    
    #TODO: sequnce with only 3 valid frames should work but does not
    if ignoreRef:
      data = data[:, :, 4:-4, 4:-4]

#    threshold = 66000#45000 #for slow mode
    
    #put a zero everywhere in the 4D data cube that the intensity is above or equal to threshold 
    #doesn't account for when data goes below threshold again after saturating
    if mode == "SLOW":
      unsatdata = np.where(data >= threshold, data, 0)
    else:
      unsatdata = np.where(data <= threshold, data, 0) 
    #make a mask data set, where data above threshold = 0, and data below threshold = 1
    mask = np.where(unsatdata == 0, unsatdata, 1) 
    #1-D array of the frame times from 0 to N 
    t = np.arange(len(data[0,:,0,0])) 
    #read in the inverse matrix elements given the number of frames in the data set
    Cinv = inversematrices(data)
   
    #create a for loop over the frame numbers, starting from frame 1
    bij = np.zeros((len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    cij = np.zeros((len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:]))) 
    lastUnsatFrame = np.zeros((len(data[:,0,0,0]), len(data[0,:,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    invCijperpix = np.zeros((6, len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    for i in t[1:]:
        #makes sure that pixels that cross the threshold say out of our calculations even if they go under threshold again
        mask[:,i] = mask[:,i]*mask[:, i-1]    
        #FIXED: next line ...*mask[:,i] not i-1
        unsatdata[:,i] = unsatdata[:,i]*mask[:, i]
        #dot the data with time vectors of order 1, and 2
        bij = (unsatdata[:,i]*i) + bij
        cij = (unsatdata[:,i]*(i**2)) + cij
        #put a one at the frame time when the pixel was last below threshold
        lastUnsatFrame[:,i-1] = np.abs(mask[:,i] - mask[:, i-1]) #make sure it is full of 1s and 0s not -1s and 0s
        #if a pixel never saturated, put a 1 in the in the last frame of the lasUnsatFrame
        lastUnsatFrame[:, t[-1]] = mask[:, t[-1]]

    #find the correct inverse matrix elements to use for each pixel (depending on time that the pixel was last below saturation)
    indMatElements = np.arange(6)
    for i in t:
        for j in indMatElements:
            invCijperpix[j] = (lastUnsatFrame[:,i]*(Cinv[i, j])) + invCijperpix[j]

    #dot the data with the time vector of order 0
    aij = unsatdata.sum(axis=1)

    #create an array of the polynomial fit coefficients for each pixel
    coefs = np.zeros((3, len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    coefs[0] = (invCijperpix[0]*aij) + (invCijperpix[3]*bij) + (invCijperpix[4]*cij)
    coefs[1] = (invCijperpix[3]*aij) + (invCijperpix[1]*bij) + (invCijperpix[5]*cij)
    coefs[2] = (invCijperpix[4]*aij) + (invCijperpix[5]*bij) + (invCijperpix[2]*cij)
        
    return coefs[0], coefs[1], coefs[2]
"""
Notes 
-----
  timing for a 10 NDR ramp of 2048x2048 the timing is as follows
    -up to where fit is done 0.69 seconds (mostly due to reshaping)
    -up to where fited data points ar calculated 2.5 seconds
    -end 2.5 seconds
  NDR = 4 total time is 1.1 sec
    -if fitted data points are done/reshaping adds 0.6 seconds
  in comparison Quadfit function (no fitted points calculated) takes 2.14 sec
  but does not check for how many points are in series
"""

#import matplotlib.pyplot as plt
#import time

def cnPolyfit(ramp, order, mode, threshold):
  """
  returns coefficient for quadratic fit [3,2048,2048], where dim 0 is quadratic
  term
  """
  # currently only works for up to quadratic
#  start = time.time()

  # define x vector for fit
  x = np.arange(ramp.shape[0])
  
  # reshape the data
  y = ramp.reshape(ramp.shape[0],ramp.shape[1]*ramp.shape[2])
  
  # find data beyond threshold
  if mode == "SLOW":
    nanY = np.where(y <= threshold, np.nan, y) 
  else:
    nanY = np.where(y >= threshold, np.nan, y)
    
  # now find first occurance of nan along the x axis
  #BUG: This only works if there actually is a nan otherwise the first entry
  # might be the biggest entry <AFE 2018-09-21 a:AFE>
  #firstNan = nanY.argmax(axis=0)
  firstNan = np.sum(np.isfinite(nanY),axis=0)
  
  # are there any series where number of data points is not sufficient for 
  # order of polynom?
  if order == 2:
    quadratic = np.where(firstNan>order)[0]
    linear = np.where(firstNan==order)[0]
    justOne = np.where(firstNan<order)[0]
  if order == 1:
    quadratic = np.array([])
    linear = np.where(firstNan>order)[0]
    justOne = np.where(firstNan<=order)[0]
  
  # set everything after first NaN to NaN
  nanY = np.where(np.repeat(np.arange(x.shape[0])[:,None],y.shape[1],axis=1)-
                  firstNan<=0, nanY, np.nan)
  
  # for only one data point set whole series to 0
  nanY[:,justOne] = 0 
  
  # create masked array that ignores everything with NaNs
  idx = ~np.isfinite(nanY)
  maskedY = np.ma.masked_array(nanY,mask=idx)
  
  # initialize coef array
  coef = np.zeros((order+1, y.shape[1]))
  
#  one = time.time()
#  print(one - start)
  
  # now use masked polyfit
  if order == 2:
    coef[:,quadratic] = np.ma.polyfit(x, maskedY[:,quadratic], 2)
    if not linear.shape[0] == 0:
      coef[1:,linear] = np.ma.polyfit(x, maskedY[:,linear], 1)
  else:
    coef[:,linear] = np.ma.polyfit(x, maskedY[:,linear], 1)
    
#  two = time.time()
#  print(two - start)
    
  # create fited data points
#  fit = np.polynomial.polynomial.polyval(x,np.flipud(coef), tensor=True).T
  
  # make coefficients and fits for invalid series Nans
  coef[:,justOne] = np.nan
#  fit[:,justOne] = np.nan
  
  # if for order 2 only linear fit was done make quadratic coef NaN
  if order ==2:
    coef[0,linear] = np.nan
    
  # reshape the fitted values
#  fit = fit.reshape(ramp.shape[0],ramp.shape[1],ramp.shape[2])
  coef = coef.reshape(order+1,ramp.shape[1],ramp.shape[2])
#  end = time.time()
#  print(end - start)
  return coef

#test = np.float32(np.arange(2048*2048*10).reshape((10,2048,2048)))
#ndrs = 4
#test = np.ones((ndrs,2048,2048))
#test = test + test*np.arange(ndrs)[:,None,None]+ test*(0.5*np.arange(ndrs)[:,None,None]**2)
##test[5,0,0]=400
#test[2,0,0]=300
#test[3,0,0]=90
##test[4,0,0]=20
#test[1,2,0]=300
##y = test.reshape(test.shape[0],test.shape[1]*test.shape[2])
#
##start = time.time()
#
#res = cnPolyfit(test, 2, "FAST", 200.)
#end = time.time()
#print(end - start)
#ndrs = 4
#test = np.ones((ndrs,2048,2048))
#test = test + test*np.arange(ndrs)[:,None,None]+ test*(0.5*np.arange(ndrs)[:,None,None]**2)
##test[5,0,0]=400
#test[2,0,0]=300
#test[3,0,0]=90
##test[4,0,0]=20
#test[1,2,0]=300
#start = time.time()
#res = matrixQuadfit(test[None,:,:,:],threshold=200, mode="FAST",ignoreRef=False)
#end = time.time()
#print(end - start)