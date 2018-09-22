#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:37:04 2018

@author: andreef
"""
import numpy as np  # Fast numeric arrays/matrices (FITS files imported as numpy matrices)
from astropy.io import fits  # Reading/writing FITS data
import glob
import matplotlib.pyplot as plt
#import InverseMatrices as m

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

    
class cnH2rgRamps():
  def __init__(self, 
                   fileString,
                   fileType,
                   readMode="SLOW", 
                   subArray=None, 
                   verbose=True):
    
    self.mode = readMode
    self.subArray = subArray
    self.verbose = verbose
    self.sequenceName = fileString
    
    # wild card search for Linux style folder structure
    if fileType is "fits":
      frameFiles = glob.glob(fileString+'*.fits')
    else:
      frameFiles = glob.glob(fileString+'*.arr')
    
    # sort file into sequences and make sure all are the same length
    self.fileList = self.sortFileList(frameFiles)   
    # creat list of input files and determind nrNDR, nSeq
    self.setFileList(self.fileList)
    
    if verbose:
      print(self.files)
      print(self.nSeq, self.nNdr)
      
    # Determine the resolution of intenral frames 
    if self.subArray != None: 
      self.xDim = self.subArray[0][1]-self.subArray[0][0]
      self.yDim = self.subArray[1][1]-self.subArray[1][0]
    else:
      self.xDim = 2048
      self.yDim = 2048
      
#  def masterInstrumentDark(self,data):
#  
#    # reference pixel correction should be done prior to feeding into this function
##    v = np.var(data,axis=0)
##    print(v.shape)
##    wo=9
##    print(np.max(v[wo,:,:]),np.min(v[wo,:,:]),np.mean(v[wo,:,:]),np.std(v[wo,:,:]))
##    print('mean variance and std of variance', np.mean(np.var(data,axis=0)),np.std(np.var(data,axis=0)))
#    averagedDark = np.average(data, axis=0)
#    
#    return averagedDark
      
  def inversematrices(self,data):
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
  
  def polyfit(self,data, order=1, time=None):
    
    assert (len(data.shape)==3), 'function can only deal with 3D data cubes'
    
    # re-arrange data to 2D
    new = data.reshape(-1, data.shape[-1])
    
    
    a=1
    
  def quadfitKathleen(self,data,threshold=66000,ignoreRef=False):  
    # From Kathleen Tatem
    #don't use the reference pixels in calculations
    if ignoreRef:
      data = data[:, :, 4:-4, 4:-4]

#    threshold = 66000#45000 #for slow mode
    
    #put a zero everywhere in the 4D data cube that the intensity is above or equal to threshold 
    #doesn't account for when data goes below threshold again after saturating
    unsatdata = np.where(data <= threshold, data, 0) 
    #make a mask data set, where data above threshold = 0, and data below threshold = 1
    mask = np.where(unsatdata == 0, unsatdata, 1) 
    #1-D array of the frame times from 0 to N 
    t = np.arange(len(data[0,:,0,0])) 
    #read in the inverse matrix elements given the number of frames in the data set
    Cinv = self.inversematrices(data)
   
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
        
    return [coefs[0], coefs[1], coefs[2]]
    
  def read(self, fileType, nSeqArr=None, dtype=np.uint16):
#    dtype = np.uint16
    if nSeqArr is None:
      nSeq = self.nSeq
      nSeqArr = range(nSeq)
    else:
      nSeq = len(nSeqArr)
    all_ims = np.zeros([nSeq, self.nNdr, self.xDim, self.yDim], dtype=dtype)
    
    for i_frame in range(self.nNdr):
      seqCounter=0
      for i_seq in nSeqArr:
        cur_file = self.files[i_seq][i_frame]
        if self.verbose:
          print('Reading', cur_file)
        if fileType is "fits":
          in_im = fits.open( cur_file )[0].data.astype(dtype)
        elif fileType is "arr":
          in_im = np.fromfile(cur_file,dtype=dtype)
          in_im = np.reshape(in_im,(self.xDim, self.yDim))
        else:
          raise ValueError("Unknown file type")
        
        # Trim image to subARray specification
        if self.subArray != None:
          in_im = in_im[ self.subArray[1][0] : self.subArray[1][1],
                           self.subArray[0][0] : self.subArray[0][1] ] 
        all_ims[seqCounter, i_frame, :, :] = in_im
        seqCounter += 1
    return all_ims
    
    
  def setFileList(self, files=[]):
    if len(files) == 0 or len(files[0]) == 0:
      raise ValueError("Expecting sequence(s) of fits files in list of list format")
    
    # Number of frame sequences
    self.nSeq = len(files)
    
    # Check to make sure all file sequences are the same length
    seq_lens = np.array( [ len( seq ) for seq in files ] )
    if np.all(seq_lens == seq_lens[0]):
      self.nNdr = int(seq_lens[0])
      if self.verbose:
        print("%d sequences of %d frames found..." % (self.nSeq, self.nNdr))
    else:
      raise ValueError("Not all sequences the same length")

#    # Change the resulting sequence length to account for refernece frame start
#    self.nNdr = self.nNdr - 1 - self.refFrame

    # Complain if not enough frames
    if self.nNdr < 1:
      raise ValueError(str("Not enough images found to compute anything (reference frame set to %d)" % (self.refFrame)) )
    
    self.files = files
    
  def sortFileList(self,frame_files):
    # - - - Sort files into sequences, make sure all are the same length
    # Assuming all fileneames are the same format 
    # [path][prefix]-[seq nu (2+ digits)-[frame number (4 digits)].fits
    # (08-29 Found out not true) The sequence number should always be [-12:-10] in the string
    # Have to use a strategy of splitting on dashes to get the sequence number
    
    # TODO: add some error checking to make sure correct format
    
    # Strip all the sequence numbers from the files
    # Store them in a set ( keeps only the unique ones )
    # 2017-07-28: Turn them back into a list and sort them so we get in order
    # 2017-08-29: Sequences aren't always 2 digits, split on dashes
    uniq_seq_nums = list(set( [ x.split("-")[-2] for x in frame_files ] ))
    uniq_seq_nums.sort()
    
   # nSeq = len(uniq_seq_nums)
    
    files = []
    
    for i_seq, seq_num in enumerate(uniq_seq_nums):
      files.append([])
  
      for cur_file in frame_files:
        if cur_file.split("-")[-2] == seq_num:
          files[i_seq].append(cur_file)
      files[i_seq].sort()
    
    return files

  def subtractSignal(self, mode='linear'):
    print('yes')        
    
  def writeFits(self,n,filepath,filename,overwrite=False):
    hdu = fits.PrimaryHDU(n)
    #hdul = fits.HDUList([hdu])
    hdu.writeto(filepath+filename,overwrite=overwrite)
    
  def writeBinary(self,data,filepath,filename):
    data.tofile(filepath+filename,sep="")
        
#a=cnH2rgRamps("data/slowdark*",readMode="SLOW",subArray=None,verbose=True)
#a=cnH2rgRamps("data/simdark*",readMode="SLOW",subArray=None,verbose=True)
#b=a.read()
#c=a.masterInstrumentDark(b)
#c = a.quadfit(np.float32(b))
#
#fig, ax=plt.subplots(1,3)
#plt.subplot(1,3,1,title=str(0))
#plt.hist(np.ndarray.flatten(c[0]),bins=1000, range=[45000.,55000.])
#plt.subplot(1,3,2,title=str(1))
#plt.hist(np.ndarray.flatten(c[1]),bins=1000, range=[-250.,0.])
#plt.subplot(1,3,3,title=str(2))
#plt.hist(np.ndarray.flatten(c[2]),bins=1000, range=[-5.,5.])
#
#t=b[0,:,:,:4]
#refCol = np.concatenate((b[:,:,:,:4],b[:,:,:,-4:]),axis=3)
#refRow = np.concatenate((b[:,:,:4,:],b[:,:,-4:,:]),axis=2)
#
#refColSubtracted= -1.*(np.float32(refCol)-np.expand_dims(np.float32(refCol[:,0,:,:]),axis=1))
#
## investigate reference pixels
## look at up the ramp raw reference columns
#ramp = 0
#fig, ax=plt.subplots(2,4)
#fig.suptitle("raw up the ramp ref columns (ramp "+str(ramp)+")", fontsize=16)
#for i in range(8):
#  plt.subplot(2,4,i+1,title=str(i))
#  plt.plot(refCol[ramp,:,:,i].T)
#
#fig, ax=plt.subplots(2,4)
#fig.suptitle("raw up the ramp ref columns (no first row)(ramp "+str(ramp)+")", fontsize=16)
#for i in range(8):
#  plt.subplot(2,4,i+1,title=str(i))
#  plt.plot(refCol[ramp,:,1:,i].T)
#
#
## look at up the ramp reference columns bias subtracted
#fig, ax=plt.subplots(2,4)
#fig.suptitle("up the ramp ref columnsbias subtracted(ramp "+str(ramp)+")", fontsize=16)
#for i in range(8):
#  plt.subplot(2,4,i+1,title=str(i))
#  plt.plot(refColSubtracted[ramp,:,:,i].T)
#
#fig, ax=plt.subplots(2,4)
#fig.suptitle("up the ramp ref columnsbias subtracted (no first row)(ramp "+str(ramp)+")", fontsize=16)
#for i in range(8):
#  plt.subplot(2,4,i+1,title=str(i))
#  plt.plot(refColSubtracted[ramp,:,1:,i].T)
#
#refCol = np.concatenate((b[:,:,:,:4],b[:,:,:,-4:]),axis=3)
#refColSubtracted= -1.*(np.float32(refCol)-np.float32(refCol[0,:,:,:]))
#
## look at raw reference columns accross ramps
#frame = 1
#fig, ax=plt.subplots(2,4)
#fig.suptitle("raw across ramps, 1st frame (subtracted first ramp)(frame "+str(frame)+")", fontsize=16)
#for i in range(8):
#  plt.subplot(2,4,i+1,title=str(i))
#  plt.text(1,0,str(np.std(refColSubtracted[:,0,1,i])))
#  plt.plot(refColSubtracted[:,frame,:,i].T)

