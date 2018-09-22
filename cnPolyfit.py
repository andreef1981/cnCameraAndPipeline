#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 14:22:22 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

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

import numpy as np
import matplotlib.pyplot as plt
import time

def cnPolyfit(ramp, order, mode, threshold):
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
ndrs = 4
test = np.ones((ndrs,2048,2048))
test = test + test*np.arange(ndrs)[:,None,None]+ test*(0.5*np.arange(ndrs)[:,None,None]**2)
#test[5,0,0]=400
test[2,0,0]=300
test[3,0,0]=90
#test[4,0,0]=20
test[1,2,0]=300
#y = test.reshape(test.shape[0],test.shape[1]*test.shape[2])

start = time.time()

res = cnPolyfit(test, 2, "SLOW", 0.)
end = time.time()
#print(end - start)
