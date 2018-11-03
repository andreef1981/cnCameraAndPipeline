# -*- coding: utf-8 -*-
"""
Created on Mon Jan 8 08:17:02 2017

@author: Andre Fehlmann
"""

import numpy as np


#-------------------------------------------------------
# solve spectrograph equation to calculate alpha and beta
#-------------------------------------------------------
def alpha_beta_angle(k, n, lam, theta):
    # input:  - diffraction order k
    #         - grating ruling density n [groves/mm]
    #         - wavelength lam [um]
    #         - Littorw angle theta [degrees]
    #
    # output: - angle of incidence alpha [degrees]
    #         - angle of diffraction beta [degrees]

    beta = np.degrees( np.arcsin( (1e-3*k*n*lam) / 
                             (2. * np.cos(-1.*np.radians(theta))) ) ) - theta
    alpha = np.degrees( np.arcsin((1e-3*k*n*lam) -
                              (np.sin(np.radians(beta))) ) )

    return alpha, beta



#-------------------------------------------------------
# calculate the angular dispersion of the grating
#-------------------------------------------------------
def angular_dispersion(beta, k, n):
    # input:  - angle of diffraction beta [degrees]
    #         - diffraction order k
    #         - grating ruling density n [groves/mm]
    #
    # output: - angular dispersion dbdlam [rad/um]
    
    dbdlam = (1e-3*k*n) / np.cos(np.radians(beta))
    
    return dbdlam

#-------------------------------------------------------
# solve spectrograph equation to calculate beta
#-------------------------------------------------------
def beta_angle(k, n, lam, alpha):
    # input:  - diffraction order k
    #         - grating ruling density n [groves/mm]
    #         - wavelength lam [um]
    #         - angle of incidence alpha [degrees]
    # output: - angle of diffraction beta [degrees]
    
    bet = np.degrees( np.arcsin( (1e-3*k*n*lam) - np.sin(np.radians(alpha)) ))
    return bet


#-------------------------------------------------------
# calculate theoretical spectrograph blaze order
#-------------------------------------------------------
def blaze_order(lam, n, grating_blaze, blaze_function=False, theta=-5.5):
    # input:  - grating ruling density n [groves/mm]
    #         - wavelength lam [um]
    #         - blaze angle t of the grating [degrees]
    #         - optional blaze_function keyword
    #         - optional littrow angle theta [degrees]
    #
    # output: - blaze order
    #         - estimated max efficiency of order

    blorder = np.array([np.around(np.sin(np.radians(grating_blaze))*2. /
                                   (lam*n*1e-3)), 1.0])
  
    if blaze_function == True:
        res = np.zeros((2,np.int(blorder[0]+3)))
        for k in range(0,np.int(blorder[0]+3)):
            angles = alpha_beta_angle(k, n, lam, theta)
            res[0,k] = k
            res[1,k] = gratingEfficiencyTheoretical(lam, n, k, grating_blaze,
               angles[0], angles[1])
        #ss =  np.nanmax(res[1,:])
        blorder=res[:,res[1,:] == np.nanmax(res[1,:])]
        
        blorder = np.squeeze(blorder)
    
    return blorder
  
#-------------------------------------------------------
# calculate the linear dispersion
#-------------------------------------------------------
def linear_dispersion(k, n, beta, Lb, inclination_angle=None,
  lh=None, angular_dispersion=None):
    # input:  - grating ruling density n [groves/mm]
    #         - diffraction order k
    #         - angle of diffraction beta [degrees]
    #         - camera focal length Lb [mm]
    #   optional  - inclination angle in [degrees]
    #             - Perpendicular distance from the spectral plane to grating
    #               Lh [mm]
    #             - angular dispersion dbdlam [rad/um]
    # output: - linear dispersion dlamdx [um/mm]

    if inclination_angle is not None:
        dlamdx = (1e3*np.cos(np.radians(beta))*
                  (np.cos(np.radians(inclination_angle))**2.)) / (k*n*lh)
    else:
        if angular_dispersion is not None:
            dlamdx = 1. /(angular_dispersion) / Lb
        else:
            dlamdx = (1e3*np.cos(np.radians(beta))) / (k*n*Lb)
    return dlamdx

def gratingEfficiencyTheoretical(lam,n,k,gratingBlaze,alpha,beta):
      # !function is only working for scalars 
      
      # input:  - grating ruling density n [groves/mm]
      #         - wavelength lam [um]
      #         - grating order k
      #         - blaze angle t of the grating [degrees]
      #         - incidence angle alpha [degrees]
      #
      #
      # output: - estimated efficiency of order
      # d = 1.d /n # in mm
       # Don Mickey
      gam = ((n*np.pi*np.cos(np.deg2rad(gratingBlaze))/lam) *
             (np.sin(np.deg2rad(beta-gratingBlaze)) + 
              np.sin(np.deg2rad(alpha-gratingBlaze))))
      res = np.square(np.sin(gam)/gam)
      
      #        if alpha < gratingBlaze:
      #            rho = np.cos(np.radians(gratingBlaze))
      #        else:
      #            rho = (np.cos(np.radians(gratingBlaze)) / 
      #                   np.cos(np.radians((alpha -gratingBlaze))))
      #      
      #      
      #        s = ((np.cos(np.radians(beta)) * 
      #                      np.cos(np.radians((alpha-gratingBlaze)))) / 
      #          (np.cos(np.radians(alpha))*np.cos(np.radians(
      #                  (beta-gratingBlaze)))))
      #        
      #        fak = np.nanmin([s,1.])
      #        
      #        if alpha + beta == 0 or alpha + beta == np.nan:
      #            res = 0
      #        else:
      #            h = ((np.pi*k*rho)*((np.cos(np.radians(gratingBlaze))) - 
      #                 (np.sin(np.radians(gratingBlaze)) /
      #                  np.tan(np.radians(alpha+beta))/2.) ) ) 
      #          
      #            res = fak * np.square(np.sin(h)/h)

      return res
# some fixed design values for spectrograph
pixel_size = 18.            # pixel size of H2RG [um]
f_number_feed = 18.         # f number of feed to spectrograph dimensionless
d_primary = 4000.           # telescope primary mirror diameter in [mm]
f_collimator = 2096.        # collimator mirror focal length [mm]
f_camera = 932.             # camera mirror focal length in [mm]
n = 31.6                    # grove density in [groves/mm]
grating_blaze = 63.9        # blaze angle in [deg] 63.9
littrow = -5.5              # littrow angle in [deg]
grating_width = 408.        # grating width in mm
grating_height = 153.       # grating heigth in mm
spectral_pixels = 1024.     # number of pixels in spectral direction

# the slit width defines the spatial sampling
slit_width = 175.           # (52 or 175 [um]
slit_height = 81.           # slit height in mm (42 or 81)

# select wavelength
# lam = 3.9343
chosenOrder = 48
lam = 1.07725
# blaze order is an integer
bl_order, eff = blaze_order (lam, n, grating_blaze,blaze_function=True,theta=littrow)

print(bl_order, eff)

# angles are in degrees

alpha, beta = alpha_beta_angle(bl_order,n,lam,littrow)
alphaC, betaC = alpha_beta_angle(chosenOrder,n,lam,littrow)
chosenEfficiency = gratingEfficiencyTheoretical(lam,n,chosenOrder,grating_blaze,alphaC,betaC)
print(alpha,beta)
print(alpha-(alpha-beta)/2.)
print(alphaC,betaC)
print(alphaC-(alphaC-betaC)/2.)
print(alphaC-60.9612,betaC-60.9612)
print('Effciency of chosen order is: ',chosenEfficiency)
# linear dispersion is in nm/um 
dlin = linear_dispersion(bl_order,n,beta, f_camera)
dlinC = linear_dispersion(chosenOrder,n,betaC, f_camera)
print(dlin)
print(dlinC,dlinC*18)
print('range on array', dlin*pixel_size*spectral_pixels, ' nm')
print('range on array with custom order', dlinC*pixel_size*spectral_pixels, ' nm')

vec = (np.arange(1024)-512)*dlinC*18.+lam*1000.
np.save("siWavelength.npy",vec)




