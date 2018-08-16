#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 12:54:56 2018

@author: andreef
"""

 """
  Returns the instrument and thermal dark subtracted ramp of the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be dark corrected.
    instrumentDark : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that contains the instrument dark ramp.
    thermalDark : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that contains the thermal dark ramp.

    Returns
    -------
    darkSubtracted : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that is dark corrected.

    Other Parameters
    ----------------
    

    Raises
    ------
    AssertationError
        If the shapes of the input arrays do not match.

    See Also
    --------

    Notes
    -----
    

    References
    ----------
    

    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> data = np.zeros((5,10,10),dtype='uint16')+6
    >>> instrumentDark = np.zeros((5,10,10),dtype='float32')+1
    >>> thermalDark = np.zeros((5,10,10),dtype='float32')+2
    >>> darkSubtracted = dark(data,instrumentDark,thermalDark)
   """