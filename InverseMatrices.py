#Inverse Matrices
#Calculates the inverse matrices needed to do a quadratic fit.
#Kathleen V. Tatem
#May 23, 2018

#Description:
    #Each pixel might saturate at the different frame along the ramp. The time matrix is a 3 x 3 symmetric matrix. The elements of the time matrix are the dot products of the time vectors to 0th order(1s), 1st order(t), and 2nd order(t**2). Since it is symmetric, it has 6 unique elements, t0*t0, t0*t1 = t1*t0, t1*t1, t1*t2 = t2*t1, t0*t2 = t2*t0, t2*t2. This code calculates the 6 unique elements of the inverse of the time matrix at all possible time values in the ramp. 
    
#Dependencies:
    #Does not depend on other routines. 
    
#Limitations:
    #There is no inverse matrix for t = 0, 1, or 2, so those inverse matrix elements are set to 0. 
    
#Input: 
    #In general, the data will be read in as a 4D array, with the 0th axis being the number of ramps, the 1st axis being the number of frames per ramp, and the last two axes being the x-y coordinates of each pixel in the array. 

#Output:
    #This code outputs a 2D array called Cinv, which has dimensios (number of frames per ramp x 6). There are 6 unique inverse matrix elements associated with each frame number (time value) of the ramp. 

import numpy as np

def inversematrices(data):
    
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
