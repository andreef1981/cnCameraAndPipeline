#QuadraticFit_SlowMode 
#Calculates the 3 coefficients of a quadratic fit function for each pixel. The linear coefficient is the flux. This codes uses the saturation threshold determined for slow mode. 
#Kathleen V. Tatem
#May 23, 2018

#Description:
        #The linear coffecient of polynomial fit function for the intensity over time can be used to determine the flux at a pixel. Instead of doing a time/memory-consuming loop over every pixel using a numpy polyfit function, this code does a quadratic fit for each pixel by only looping over the time variable (number of frames in the ramp). Each image is 2048x2048 pixels, but the time of the ramp is usually on the order of 20 frames or less. Since we are interested linear term, that is the output. We only use the pixels that are below saturation level, and never use them again even if they go below saturation level later. We also don't use the reference pixels. You should use this code of data that has already been calibrated (corrected based on the reference pixels, and with the first two frames corrected based on dark images.)
    
#Dependencies:
    #Depends on the InverseMatrices.py routine, which produces an array of the 6 unique elements of the inverse fit matrix for each possible time value. 
    #Read in correctd data.
    
#Limitations:
    #There is no inverse matrix for t = 0, 1, or 2, so those inverse matrix elements are set to 0.
    #This code does not include a mask of bad pixels or account for them in any way.
    #Once a pixel crosses the threshold (saturates), it is not included in calculations. 
    
#Input: 
    #Use data that has already been calibrated. The data will be read in as a 4D array, with the 0th axis being the number of ramps (which can be just 1), the 1st axis being the number of frames per ramp, and the last two axes being the x-y coordinates of each pixel in the array. 
    
#Output:
    #This code produces a 2D image array of flux values (the linear coefficient) for each pixel. 

import numpy as np
import InverseMatrices as m

def quadfit(data):    
    #don't use the reference pixels in calculations
    data = data[:, :, 4:-4, 4:-4]

    threshold = 66000#45000 #for slow mode
    
    #put a zero everywhere in the 4D data cube that the intensity is above or equal to threshold 
    #doesn't account for when data goes below threshold again after saturating
    unsatdata = np.where(data <= threshold, data, 0) 
    #make a mask data set, where data above threshold = 0, and data below threshold = 1
    mask = np.where(unsatdata == 0, unsatdata, 1) 
    #1-D array of the frame times from 0 to N 
    t = np.arange(len(data[0,:,0,0])) 
    #read in the inverse matrix elements given the number of frames in the data set
    Cinv = m.inversematrices(data)
   
    #create a for loop over the frame numbers, starting from frame 1
    bij = np.zeros((len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    cij = np.zeros((len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:]))) 
    lastUnsatFrame = np.zeros((len(data[:,0,0,0]), len(data[0,:,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    invCijperpix = np.zeros((6, len(data[:,0,0,0]), len(data[0,0,:,0]), len(data[0,0,0,:])))
    for i in t[1:]:
        #makes sure that pixels that cross the threshold say out of our calculations even if they go under threshold again
        mask[:,i] = mask[:,i]*mask[:, i-1]        
        unsatdata[:,i] = unsatdata[:,i]*mask[:, i-1]
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
        
    return [coefs[0], coefs[1]]
