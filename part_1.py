# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('bandgap.txt')
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
plt.show()