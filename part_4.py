#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:27:16 2017

@author: punit
"""

 

import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('carbon1.txt')
#eigenvalues1.sort(axis=0)
ee=eigenvalues1[::]
x1=np.linspace(-10,10,len(ee) )
#ee=eigenvalues1[::10]
#ee1=eigenvalues1[1::10]
#ee2=eigenvalues1[2::10]
#ee3=eigenvalues1[3::10]
#ee4=eigenvalues1[4::10]
#ee5=eigenvalues1[5::10]
#ee6=eigenvalues1[6::10]
#ee7=eigenvalues1[7::10]
#ee8=eigenvalues1[8::10]
#ee9=eigenvalues1[9::10]
#x1=np.linspace(-10,10,len(ee) )
#print(len(ee))
#print(len(ee1))
#print(len(ee2))
#print(len(ee3))
plt.plot(x1,ee)
#plt.plot(x1,ee1)
#plt.plot(x1,ee2)
#plt.plot(x1,ee3)
#plt.plot(x1,ee4)
#plt.plot(x1,ee5)
#plt.plot(x1,ee6)
#plt.plot(x1,ee7)
#plt.plot(x1,ee8)
#plt.plot(x1,ee9)
plt.xlabel('k -->')
plt.xlim(-10,10)
plt.show()