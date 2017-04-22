 

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import Axes3D

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
poly=np.loadtxt('poly1.txt')
plt.plot(poly[:,0],poly[:,1],color='r', label='Total Energy')
plt.plot(poly[:,0],poly[:,2],color='g', label='Energy')
plt.xlabel('Delta -->')
plt.ylabel('Energies -->')
plt.legend()
l.show()

gg = plt.figure(4)
eigenvalues11=np.loadtxt('carbon1.txt')
ee=eigenvalues11[::]
x1=np.linspace(-3.14,3.14,len(ee) )
plt.plot(x1,ee)
plt.xlabel('k -->')
plt.xlim(-3.14,3.14)
plt.show()
gg.show()



gg1 = plt.figure(5)
k=np.loadtxt('H2D.txt')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x=k[:,0]
y=k[:,1]
z=k[:,2]
z2=k[:,3]
ax.plot(x, y, z, color='r')
ax.plot(x, y, z2, color='g')
gg1.show()












