import numpy as np
import matplotlib.pyplot as plt

file='testing'
unit=4
time='-test'

f = open(file,'r')
nx = f.readline()
ny = f.readline()
nx=int(nx)+1
ny = int(ny)+1
print ny/float(nx)

data = np.genfromtxt(file,skip_header=2)
print data
nsteps = (data.shape[0])/(nx*ny+2)
print nsteps

U = [0]*nsteps
V = [0]*nsteps
P = [0]*nsteps
T = [0]*nsteps
rho = [0]*nsteps
vel = [0]*nsteps

for i in range(nsteps):
    P[i] = np.zeros((nx,ny))
    U[i] = np.zeros((nx,ny))
    V[i] = np.zeros((nx,ny))
    T[i] = np.zeros((nx,ny))
    rho[i] = np.zeros((nx,ny))
    vel[i] = np.zeros((nx,ny))
    for j in range(ny):
        P[i][:,j]=data[i*(nx*ny+2)+2+nx*j:i*(nx*ny+2)+2+nx*(j+1),0]
        U[i][:,j]=data[i*(nx*ny+2)+2+nx*j:i*(nx*ny+2)+2+nx*(j+1),1]
        V[i][:,j]=data[i*(nx*ny+2)+2+nx*j:i*(nx*ny+2)+2+nx*(j+1),2]
        T[i][:,j]=data[i*(nx*ny+2)+2+nx*j:i*(nx*ny+2)+2+nx*(j+1),3]
        rho[i][:,j]=data[i*(nx*ny+2)+2+nx*j:i*(nx*ny+2)+2+nx*(j+1),4]
    vel[i] = np.sqrt(U[i]**2 + V[i]**2)
        
x = np.linspace(0,1,nx)
y = np.linspace(0,ny/float(nx),ny)
xx, yy = np.meshgrid(x,y)


plt.figure()
plot2 = plt.contourf(xx,yy,np.transpose(U[unit]))#,np.arange(-1,1,0.01), extend='both')
plt.colorbar(plot2)
#plt.axis([0, 1, 0, ny/float(nx)])
plt.savefig('contour-U' + time + '.png')

plt.figure()
plot2 = plt.contourf(xx,yy,np.transpose(V[unit]),np.arange(-0.5,0.5,0.01), extend='both')
plt.colorbar(plot2)
plt.axis([0, 1, 0, ny/float(nx)])
plt.savefig('contour-V' + time + '.png')

plt.figure()
plot2 = plt.contourf(xx,yy,np.transpose(P[unit]),200)
plt.colorbar(plot2)
plt.axis([0, 1, 0, ny/float(nx)])
plt.savefig('contour-P' + time + '.png')

plt.figure()
plot2 = plt.contourf(xx,yy,np.transpose(T[unit]),50)
plt.colorbar(plot2)
plt.axis([0, 1, 0, ny/float(nx)])
plt.savefig('contour-T' + time + '.png')

plt.figure()
plot2 = plt.contourf(xx,yy,np.transpose(rho[unit]),200)
plt.colorbar(plot2)
plt.axis([0, 1, 0, ny/float(nx)])
plt.savefig('contour-rho' + time + '.png')

plt.figure()
plot1 = plt.streamplot(xx,yy,np.transpose(U[unit]),np.transpose(V[unit]),color=np.transpose(vel[unit]))
plt.colorbar()
plt.axis([0, 1, 0, ny/float(nx)])
plt.savefig('streamplot' + time + '.png')

