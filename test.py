#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 12:33:51 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

from skimage import viewer,data, io, filters


image = data.coins()

new_viewer = viewer.ImageViewer(image) 

new_viewer.show() 