 

import matplotlib.pyplot as plt
import numpy as np






f = plt.figure(1)
eigenvalues1=np.loadtxt('carbon.txt')
ee=eigenvalues1[::]
x1=np.linspace(-10,10,len(ee) )
plt.plot(x1,ee)
plt.xlabel('k -->')
plt.xlim(-10,10)
f.show()

g = plt.figure(2)
ben=np.loadtxt('ben_sym.txt')
plt.plot(ben[:,0],ben[:,1],color='r', label='Benzene')
plt.plot(ben[:,0],ben[:,2], color='g',label='Benzene Cation')
plt.xlim(-0.25,0.25)
plt.xlabel('Delta -->')
plt.ylabel('Energies -->')
plt.legend()
g.show()

l = plt.figure(3)
poly=np.loadtxt('poly.txt')
plt.plot(poly[:,0],poly[:,1],color='r', label='Total Energy')
#plt.plot(ben[:,0],ben[:,2], color='g',label='Benzene Cation')
#plt.xlim(-0.25,0.25)
plt.xlabel('Delta -->')
plt.ylabel('Energies -->')
plt.legend()
l.show()