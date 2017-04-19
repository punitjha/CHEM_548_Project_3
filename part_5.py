#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 13:45:33 2017

@author: punit
"""

#========================
#3D surface (solid color)
#========================
#
#Demonstrates a very basic plot of a 3D surface using a solid color.
#'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
print(x)
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
#theta=np.linspace(-4*np.pi, 4*np.pi, 100)
#z=np.linspace(-2,2,100)
#x= np.sin(theta)
#y=np.cos(theta)
# Plot the surface
ax.plot_surface(x, y, z, color='b')

plt.show()