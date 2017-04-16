 

import matplotlib.pyplot as plt
import numpy as np
eigenvalues1=np.loadtxt('carbon.txt')
ee=eigenvalues1[::4]
ee1=eigenvalues1[1::4]
ee2=eigenvalues1[2::4]
ee3=eigenvalues1[3::4]
x1=np.linspace(-10,10,len(ee) )
print(len(ee))
print(len(ee1))
print(len(ee2))
print(len(ee3))
plt.plot(x1,ee)
plt.plot(x1,ee1)
plt.plot(x1,ee2)
plt.plot(x1,ee3)
plt.xlabel('k -->')
plt.xlim(-10,10)
plt.show()