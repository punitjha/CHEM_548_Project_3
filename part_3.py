 

import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('carbon.txt')
ee=eigenvalues1[::]
x1=np.linspace(-10,10,len(ee) )
plt.plot(x1,ee)
plt.xlabel('k -->')
plt.xlim(-10,10)
plt.show()