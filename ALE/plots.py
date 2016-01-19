import numpy as np
import matplotlib.pyplot as plt

file='testing'
unit=49
time='-test'

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
T = [0]*nsteps
rho = [0]*nsteps
vel = [0]*nsteps
y_g = [0]*nsteps

for i in range(nsteps):
    P[i] = np.zeros((nx,nx))
    U[i] = np.zeros((nx,nx))
    V[i] = np.zeros((nx,nx))
    T[i] = np.zeros((nx,nx))
    rho[i] = np.zeros((nx,nx))
    vel[i] = np.zeros((nx,nx))
    y_g[i] = np.zeros((nx,nx))
    for j in range(nx):
        P[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),0]
        U[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),1]
        V[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),2]
        T[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),3]
        rho[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),4]
        y_g[i][:,j]=data[i*(nx*nx+2)+2+nx*j:i*(nx*nx+2)+2+nx*(j+1),5]
    vel[i] = np.sqrt(U[i]**2 + V[i]**2)
        
x = np.linspace(0,1,nx)
#y = np.linspace(0,1,nx)
#y = np.linspace(0,np.amax(y_g[unit]),nx)
#print np.amax(y_g[unit])
#xx, yy = np.meshgrid(x,y)


for unit in range(nsteps):
    y = np.linspace(0,np.amax(y_g[unit]),nx)
    xx, yy = np.meshgrid(x,y)

    plt.figure()
    plot2 = plt.contourf(xx,yy,U[unit],np.arange(-0.5,0.5,0.01), extend='both')
    plt.colorbar(plot2)
    plt.axis([0, 1, 0, 1.4])
    plt.savefig('contour-U' + str(unit) + '.png')
    plt.close()

    plt.figure()
    plot2 = plt.contourf(xx,yy,V[unit],np.arange(-1,1,0.01), extend='both')
    plt.colorbar(plot2)
    plt.axis([0, 1, 0, 1.4])
    plt.savefig('contour-V' + str(unit) + '.png')
    plt.close()

    plt.figure()
    plot2 = plt.contourf(xx,yy,P[unit],np.arange(0,200000,1000), extend='both')
    plt.colorbar(plot2)
    plt.axis([0, 1, 0, 1.4])
    plt.savefig('contour-P' + str(unit) + '.png')
    plt.close()
    
    plt.figure()
    plot2 = plt.contourf(xx,yy,T[unit],50)
    plt.colorbar(plot2)
    plt.axis([0, 1, 0, 1.4])
    plt.savefig('contour-T' + str(unit) + '.png')
    plt.close()
    
    print 'hey'
#    plt.figure()
#    plot2 = plt.contourf(xx,yy,rho[unit],200)
#    plt.colorbar(plot2)
#    plt.axis([0, 1, 0, 1.4])
#    plt.savefig('contour-rho' + str(unit) + '.png')
    
#plt.figure()
#plot1 = plt.streamplot(xx,yy,U[unit],V[unit],color=vel[unit])
#plt.colorbar()
#plt.axis([0, 1, 0, 1.4])
#plt.savefig('streamplot' + time + '.png')

