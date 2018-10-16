# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 20:53:03 2018

@author: Andre
"""

t = np.arange(10)*0.36
a=1000.
b=5000.
c=0
y = a*t**2 + b*t + c
l = b*t + 0

fig, ax=plt.subplots()
ax.plot(t,y,'b',t,l,'r')
