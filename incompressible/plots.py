import numpy as np
import matplotlib.pyplot as plt

file='test_incompressible_nx64_dt1.00E-05_tend1.00E-02'
file = 'ss1.0'
unit=0
time='-1.0ss'

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


plt.figure()
plot2 = plt.contourf(xx,yy,U[unit],np.arange(-1,1,0.01), extend='both')
plt.colorbar(plot2)
plt.title('U', size ='20')
plt.axis([0, 1, 0, 1])
plt.savefig('contour-U' + time + '.png')

plt.figure()
plot2 = plt.contourf(xx,yy,V[unit],np.arange(-0.5,0.5,0.01), extend='both')
cb = plt.colorbar(plot2)
plt.title('V', size ='20')
plt.axis([0, 1, 0, 1])
plt.savefig('contour-V' + time + '.png')

plt.figure()
plot2 = plt.contourf(xx,yy,P[unit],50)
plt.colorbar(plot2)
plt.title('P', size ='20')
plt.axis([0, 1, 0, 1])
plt.savefig('contour-P' + time + '.png')

plt.figure()
plot1 = plt.streamplot(xx,yy,U[unit],V[unit],color=vel[unit])
plt.colorbar()
plt.title('Stream Lines', size ='20')
plt.axis([0, 1, 0, 1])
plt.savefig('streamplot' + time + '.png')

