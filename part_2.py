# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:07:19 2017

@author: mintuser
"""

import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('bandgap1.txt')
ee=[]
ee1=[]
print(len(eigenvalues1))
x=np.linspace(-10,10,20)
for row in range(40):
        if(row%2 == 0):
            ee.append(eigenvalues1[row])
        else:
            ee1.append(eigenvalues1[row])   
print(len(ee))
plt.plot(x,ee)
plt.plot(x,ee1)
plt.xlabel('k -->')
plt.xlim(-10,10)
plt.show()