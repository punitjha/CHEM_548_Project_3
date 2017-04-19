#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 13:45:33 2017

@author: punit
"""



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
k=np.loadtxt('H2D.txt')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x=k[:,0]
y=k[:,1]
z=k[:,2]
z2=k[:,3]
ax.plot(x, y, z, color='r')
ax.plot(x, y, z2, color='g')
plt.zlim(-5,5)
plt.show()