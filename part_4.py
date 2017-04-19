#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:27:16 2017

@author: punit
"""

 

import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('carbon1.txt')
ee=eigenvalues1[::]
x1=np.linspace(-3.14,3.14,len(ee) )
plt.plot(x1,ee)
plt.xlabel('k -->')
plt.xlim(-3.14,3.14)
plt.show()