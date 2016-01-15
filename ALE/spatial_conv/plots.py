import numpy as np
import matplotlib.pyplot as plt

file='test_incompressible_nx32_dt1.00E-07_tend2.00E-06'
#file = 'testing.out'
f = open(file,'r')
nx = f.readline()
f.close()
nx=int(nx)+1

data = np.genfromtxt(file,skip_header=1)
print data
nsteps = (data.shape[0])/(nx*nx+2)
print nsteps

U = [0]*nsteps
V = [0]*nsteps
P = [0]*nsteps
vel = [0]*nsteps

for i in range(nsteps):
    P[i] = np.zeros((nx,nx))
    U[i] = np.zeros((nx,nx))
    V[i] = np.zeros((nx,nx))
    vel[i] = np.zeros((nx,nx))
    for j in range(nx):
        P[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),0]
        U[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),1]
        V[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),2]        
    vel[i] = np.sqrt(U[i]**2 + V[i]**2)
        
x = np.linspace(0,1,nx)
y = x
xx, yy = np.meshgrid(x,y)

unit=0

plt.figure()
plot2 = plt.contourf(xx,yy,U[unit],50)
plt.colorbar(plot2)
plt.axis([0, 1, 0, 1])
plt.savefig('contour-U.png')

plt.figure()
plot2 = plt.contourf(xx,yy,V[unit],50)
plt.colorbar(plot2)
plt.axis([0, 1, 0, 1])
plt.savefig('contour-V.png')

plt.figure()
plot2 = plt.contourf(xx,yy,P[unit],50)
plt.colorbar(plot2)
plt.axis([0, 1, 0, 1])
plt.savefig('contour-P.png')

plt.figure()
plot1 = plt.streamplot(xx,yy,U[unit],V[unit],color=vel[unit])
plt.colorbar(plot2)
plt.axis([0, 1, 0, 1])
plt.savefig('streamplot.png')

