#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 13:51:37 2018

@author: andreef
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools

def main():
    x, data = generate_data(256, 6)
    model = [invert(x, y) for y in data.T]
    sigma, mu, height = [np.array(item) for item in zip(*model)]
    prediction = gaussian(x, sigma, mu, height)

    plot(x, data, linestyle='none', marker='o')
    plot(x, prediction, linestyle='-')
    plt.show()

def parallel_main():
    import multiprocessing
    p = multiprocessing.Pool()
    x, data = generate_data(256, 262)
    args = zip(itertools.repeat(x), data.T)
    model = p.imap(parallel_func, args, chunksize=500)
    sigma, mu, height = [np.array(item) for item in zip(*model)]
    prediction = gaussian(x, sigma, mu, height)
    
    plot(x, data, linestyle='none', marker='o')
    plot(x, prediction, linestyle='-')
    plt.show()

def parallel_func(args):
    return invert(*args)

def invert(x, y):
    # Use only data within the "peak" (20% of the max value...)
    key_points = y > (0.2 * y.max())
    x = x[key_points]
    y = y[key_points]

    # Fit a 2nd order polynomial to the log of the observed values
    A, B, C = np.polyfit(x, np.log(y), 2)

    # Solve for the desired parameters...
    sigma = np.sqrt(-1 / (2.0 * A))
    mu = B * sigma**2
    height = np.exp(C + 0.5 * mu**2 / sigma**2)
    return sigma, mu, height

def generate_data(numpoints, numcurves):
    np.random.seed(3)
    x = np.linspace(0, 500, numpoints)

    height = 100 * np.random.random(numcurves)
    mu = 200 * np.random.random(numcurves) + 200
    sigma = 100 * np.random.random(numcurves) + 0.1
    data = gaussian(x, sigma, mu, height)

    noise = 5 * (np.random.random(data.shape) - 0.5)
    return x, data + noise

def gaussian(x, sigma, mu, height):
    data = -np.subtract.outer(x, mu)**2 / (2 * sigma**2)
    return height * np.exp(data)

def plot(x, ydata, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    colors_ = ['r','b','g','m','k']
    colorcycle = itertools.cycle(colors_)
    for y, color in zip(ydata.T, colorcycle):
        ax.plot(x, y, color=color, **kwargs)

parallel_main()