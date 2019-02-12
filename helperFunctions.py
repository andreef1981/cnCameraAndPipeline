#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)
helper library
"""

import numpy as np  # Fast numeric arrays/matrices (FITS files imported as numpy matrices)
from scipy.signal import find_peaks
from scipy import exp
# import matplotlib.pyplot as plt

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


#import matplotlib.pyplot as plt
#import time

def cnRefSpectrum(cw,
                  slit,
                  desiredDispersion):
  
  if cw >= 3900. and cw <= 3968.:
    print ("here")
    #  3934 nm
    wavenumber = np.array([2530.4232, 2531.7493, 2542.4886, 2542.8614, 2549.9993, 2552.3161])
    waveWidth = np.array([0.014, 0.021, 0.014, 0.025, 0.016, 0.025]) #fwhm = 2.355*sigma
    intensity = np.array([2,20,9,53,3,32])
    names = ["Th I", "Ar I", "Th II", "Ar I", "Th I", "Ar I"]
  elif cw >= 1063. and cw <= 1103.:
    # 1083.0 nm
    #Engleman 2003
    # wavenumber = np.array([9136.3890, 9136,9082, 9151.2409, 9160.3878, 9162.2777,
    #                       9170.7994, 9179.9819])
    # waveWidth = np.array([0.020,0.019,0.018,0.021,0.022,0.026,0.028])
    # intensity = np.array([4215,296,111,98,80,410,464])
    # names = ["Th II", "Th II", "Th II", "Th I", "Th I", "Th I", "Th I"]
    # Hinkle 2001
    scale = 1.#np.mean(np.array([80/0.0066,98/0.0191,410/0.0326,464/0.0243]))
    wavenumber = np.array([9152.115,9160.388,9162.278,
                          9170.799,9173.53,9175.288,9176.116,9178.231,9179.982,
                          9183.69,9185.371,9186.487,9187.864,9189.64,9197.142,
                          9198.903,9199.360,9199.695,9201.688,9203.462,9204.055,
                          9204.668,9207.824,9209.04,9211.130,9217.941,9218.820,
                          9221.51,9222.53,9224.63,9226.583,9227.152,9227.527,
                          9229.46,9231.543,9235.95,9236.23,9237.31,9239.545,
                          9241.423,9242.606,9245.257,9245.683,9246.209,9250.433,
                          9250.707,9252.848,9253.365,9254.553,9256.682,9258.28,
                          9258.75,9260.11,9261.815,9268.425,9268.815,9269.478,
                          9274.675,9276.494,9279.606,9282.150])
    waveWidth = 2.3e-6*wavenumber #np.ones(wavenumber.shape)*0.6#/(2.35)
    intensity = scale*np.array([0.0610,0.0191,0.0066,
                          0.0326,0.0015,0.0013,0.0023,0.2303,0.0243,0.0101,0.0024,
                          0.0077,1.5511,0.0010,0.0415,0.0063,0.0156,0.0053,0.0064,
                          0.0230,0.0017,0.1426,0.0065,0.0017,0.0196,0.0072,0.0250,
                          0.0016,0.0026,0.0013,0.0039,0.0039,0.0111,0.0017,0.0122,
                          0.0024,0.0026,0.0097,0.0786,0.0078,0.0256,0.4545,0.1455,
                          0.0626,0.0192,0.0472,0.0034,0.1019,0.0038,0.1442,0.0013,
                          0.0023,0.0055,0.0111,0.0121,0.0557,0.0079,0.0080,0.0526,
                          0.4566,0.0736])
    names = ["Ar II","Th I","Th I",
            "Th I","Ar II","Ar I","Th I","Ar I","Th I","Ar I","Th I","Th I",
            "Ar I","Undef","Th I","Ar II","Ar II","Th I","Th I","Th I","Th II",
            "Ar I","Th I","Undef","Th I","Ar I","Th I","Undef","Ar II","Undef",
            "Th blend","Th II","Th I","Ar I","Ar II","Ar I","Undef","Ar I",
            "Th II","Ar II","Ar II","Th I","Ar II","Th I","Th I","Ar I","Th I",
            "Th I","Th I","Th I","Undef","Undef","Ar I","Th II","Th I","Th I",
            "Ar II","Th I","Th I","Ar I","Ar I"]
  
  # convert to wavelengths  
  wavelength = cnWnum2ang(wavenumber)/10.
  fwhmWavelength = waveWidth/wavenumber*wavelength
  
  # use high resolution to create reference spectrum
  artRes = 0.001
  if slit == "175um":
    slitProfile = np.int(np.ceil(6*desiredDispersion/artRes))
  else:
    slitProfile = np.int(np.ceil(2*desiredDispersion/artRes))
  
  x = np.arange(-3*slitProfile,3*slitProfile,1)
  kernel = cnGauss(x,1,0,slitProfile/(2.355),0)
  
  # refsepctrum goes 1% above beyond max of lookup values
  wl = np.arange(np.min(wavelength)-0.01*np.min(wavelength),
                 np.max(wavelength)+0.01*np.max(wavelength),artRes)
  spec = np.zeros(wl.shape)
  
  # create spectrum
  for ii in range(len(wavelength)):
    spec = spec + cnGauss(wl,intensity[ii],wavelength[ii],fwhmWavelength[ii]/(2.355),0)
    
  
  # convolve with instrument profile
  conv = np.convolve(spec,kernel,mode="same")/sum(kernel)
  # normalize to the max value
  conv = conv/np.max(conv)
  # interpolate to desired dispersion values
  allx = np.arange(np.min(wl), np.max(wl), desiredDispersion)
  allfinal = np.interp(allx,wl,conv)
  return allx, allfinal


def cnFindSpatialPeaks(im,
                       peakWidth,
                       xLoc,
                       averageFwhm,
                       prominenceLimits,
                       heightLimits,
                       invert=False):
  
  # try to estimate a maximum height from the data
  #TODO: use bad pixel mask to refine this
  if invert:
    im = -1*im
  # profile is calculated from certain x location
  if isinstance(xLoc,int):
    # averaging over several columns
    if averageFwhm>0:
      profile=np.median(im[:,xLoc-averageFwhm:xLoc+averageFwhm-1], axis=1)
    # only use a single column (could be affected by bad pixels)
    else:
      profile=im[:,xLoc]
  # use given xLoc interval for median profile
  else:
    profile=np.median(im[:,xLoc[0]:xLoc[1]], axis=1)
  peaks = find_peaks(profile,width=peakWidth,prominence=prominenceLimits,
                     height=heightLimits)
  ind = peaks[0] 
  widths = peaks[1]['widths']
  prominences = peaks[1]['prominences']
  heights = peaks[1]['peak_heights']
  
  return ind, widths, prominences, heights, profile

def cnFind2DSpatialPeak(im,
                         myAxis,
                         peakWidth,
                         xLoc,
                         averageFwhm,
                         prominenceLimits,
                         heightLimits,
                         invert=False):
  # TODO: is identical to 1D version exept the myAxis distinction before find_peaks
  # try to estimate a maximum height from the data
  #TODO: use bad pixel mask to refine this
  if invert:
    im = -1*im
  # profile is calculated from certain x location
  if isinstance(xLoc,int):
    # averaging over several columns
    if averageFwhm>0:
      profile=np.mean(im[:,xLoc-averageFwhm:xLoc+averageFwhm-1], axis=myAxis)
    # only use a single column (could be affected by bad pixels)
    else:
      profile=im[:,xLoc]
  # use given xLoc interval for median profile
  else:
    if myAxis==0:
      profile=np.mean(im[xLoc[0]:xLoc[1],:], axis=myAxis)
    elif myAxis==1:
      profile=np.mean(im[:,xLoc[0]:xLoc[1]], axis=myAxis)
  peaks = find_peaks(profile,width=peakWidth,prominence=prominenceLimits,
                     height=heightLimits)
  ind = peaks[0] 
  widths = peaks[1]['widths']
  prominences = peaks[1]['prominences']
  heights = peaks[1]['peak_heights']
  
  return ind, widths, prominences, heights, profile

def cnFindSpectralPeaks(im,
                        peakWidth,
                        spatialPeaks,
                        spatialWidths,
                        xLoc,
                        prominenceLimits,
                        heightLimits,
                        invert=False):
  #TODO: implement subpixel width support
  # assert(spatialWidths.dtype == np.int16()),\
  # "Spatial widths need to be integers for indexing"
  
  if invert:
    im = -1*im
  
  if len(np.shape(im)) > 1:  
    #special treatment if there is only one spatial profile
    if np.size(spatialPeaks) == 1:
      result = []
      # a= xLoc[0]
      # b=xLoc[1]
      spatialPeaks = spatialPeaks[0]
      spatialWidths = spatialWidths[0]
      # c = im[spatialPeaks-spatialWidths:spatialPeaks+spatialWidths-1,xLoc[0]:xLoc[1]]
      profile = np.median(im[spatialPeaks-spatialWidths:spatialPeaks+spatialWidths-1,xLoc[0]:xLoc[1]],axis=0)
      
      peaks = find_peaks(profile,width=peakWidth,
                           prominence=prominenceLimits,
                           height=heightLimits)
      ind = peaks[0] 
      widths = peaks[1]['widths']
      prominences = peaks[1]['prominences']
      heights = peaks[1]['peak_heights']
      result.append((ind,widths,prominences,heights))
    else:
      nrSpectra = len(spatialPeaks)
      profile = np.zeros((nrSpectra,xLoc[1]-xLoc[0]))
      result = []
      
      for ii in range(nrSpectra):
        profile[ii,:] = np.median(im[spatialPeaks[ii]-spatialWidths[ii]:
          spatialPeaks[ii]+spatialWidths[ii]-1,xLoc[0]:xLoc[1]],axis=0)
        peaks = find_peaks(profile[ii,:],width=peakWidth,
                           prominence=prominenceLimits,
                           height=heightLimits)
        ind = peaks[0] 
        widths = peaks[1]['widths']
        prominences = peaks[1]['prominences']
        heights = peaks[1]['peak_heights']
        result.append((ind,widths,prominences,heights))
      
    return result, profile
  else:
    result = []
    peaks = find_peaks(im,width=peakWidth,
                         prominence=prominenceLimits,
                         height=heightLimits)
    ind = peaks[0] 
    widths = peaks[1]['widths']
    prominences = peaks[1]['prominences']
    heights = peaks[1]['peak_heights']
    result.append((ind,widths,prominences,heights))
    return result


def cnGauss(x,
            a,
            x0,
            sigma,
            c):
  """
  Simple Gauss function.
  """
  return c+a*exp(-(x-x0)**2/(2*sigma**2))
  

def cnNonLinearityCorrection(data,
                             mode,
                             linThreshold,
                             multiRamp=False):  
  """
  Revision history
  ------------------
   11 Feb 2019:
     - changing fast mode to return same amount of NDRs as input
  """
  if multiRamp:
    # in this case we use the NDR number as our 'time' axis
    dataTime = np.arange(data.shape[1])
    
    linearityCorrected = np.zeros(data.shape,dtype=np.float32)
    
    if mode == "FAST":
      # fast mode does not use the first frame for fit
      data = data[:,1:,:,:]
      # time vector is shorter but must maintain the values
      orgDataTime = dataTime
      dataTime = dataTime[1:]
    else:
      orgDataTime = dataTime
    
    
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
      
      if order == 2:
        linearityCorrected[i,:,:,:] = np.multiply(coef[1,:,:],orgDataTime[:,None,None])
      else:
        linearityCorrected[i,:,:,:] = np.multiply(coef[0,:,:],orgDataTime[:,None,None])
  else:
    # in this case we use the NDR number as our 'time' axis
    dataTime = np.arange(data.shape[0])
    
    if mode is not "SLOW":
      # fast mode does not use the first frame
      data = data[1:,:,:]
      # time vector is shorter but must maintain the values
      orgDataTime = dataTime
      dataTime = dataTime[1:]
    else:
      orgDataTime = dataTime
    
    # Check for order of correction
    if len(dataTime) == 2:
      order = 1
    elif len(dataTime) > 2:
      order = 2
    else:
      raise ValueError("sequence to short to apply polyfit")
      
    # Do quadratic fit, use data up to threshold
    coef = cnPolyfit(data, order, mode, linThreshold)
    if order == 2:
      linearityCorrected = np.multiply(coef[1,:,:],orgDataTime[:,None,None])
    else:
      linearityCorrected = np.multiply(coef[0,:,:],orgDataTime[:,None,None])
      
  return linearityCorrected
  

def cnPolyfit(ramp, order, mode, threshold):
  """
  returns coefficient for quadratic fit [3,2048,2048], where dim 0 is quadratic
  term
  

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
  # currently only works for up to quadratic
#  start = time.time()

  # define x vector for fit
  x = np.arange(ramp.shape[0])
  if mode == "FAST":
    x = x+1
  
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
  
  #BUG: This does not consider cases where NaNs and values are alternating
  # <AFE 2018-10-15 a:AFE>
  #firstNan = np.sum(np.isfinite(nanY),axis=0)
  firstNan = nanY.argmax(axis=0)
  # but if there is no Nan
  tester = np.max(nanY, axis=0)
  length = y.shape[0]
  firstNan = np.where(np.isfinite(tester),length,firstNan)
  
  # set everything after first NaN to NaN
  nanY = np.where(np.repeat(np.arange(x.shape[0])[:,None],y.shape[1],axis=1)-
                  firstNan<=0, nanY, np.nan)
  
  
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
  

def cnWnum2ang(wavenum, help=False):
 
  if help is True:
    print('Convert input wavenumber in cm^-1 (vacuum) to air Angstroms')
    print('angstrom=wnum2ang(wavenum)')
    print('wavenum - floating point wavenumber in Kaysers')
    print('(cm^-1 in vacuum)')
    print('angstrom - wavelength in Angstroms computed for air')
    angstrom=0.
    return angstrom
	

  #Compute index of refraction from NBS formulation
  a=1.+6432.8e-8
  b=2949810.0
  c=146.0e8
  d=25540.0
  e=41.0e8
  n=a+(b/(c-wavenum**2.))+(d/(e-wavenum**2.))
  
  angstrom=1./(n*wavenum)*1.e8
  return angstrom


def test():
  return 1

#test= cnGauss(np.arange(100),1,50,10,0)
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