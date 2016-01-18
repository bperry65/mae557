import numpy as np
import matplotlib.pyplot as plt

#  Command to make gif out of many png
#  convert -delay 10 -loop 0 $(ls -v *V*png)  V.gif

 
file='testing'
unit=1
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
    P[i] = np.zeros((ny,nx))
    U[i] = np.zeros((ny,nx))
    V[i] = np.zeros((ny,nx))
    T[i] = np.zeros((ny,nx))
    rho[i] = np.zeros((ny,nx))
    vel[i] = np.zeros((ny,nx))
    for j in range(nx):
        P[i][:,j]=data[i*(nx*ny+2)+2+ny*j:i*(nx*ny+2)+2+ny*(j+1),0]
        U[i][:,j]=data[i*(nx*ny+2)+2+ny*j:i*(nx*ny+2)+2+ny*(j+1),1]
        V[i][:,j]=data[i*(nx*ny+2)+2+ny*j:i*(nx*ny+2)+2+ny*(j+1),2]
        T[i][:,j]=data[i*(nx*ny+2)+2+ny*j:i*(nx*ny+2)+2+ny*(j+1),3]
        rho[i][:,j]=data[i*(nx*ny+2)+2+ny*j:i*(nx*ny+2)+2+ny*(j+1),4]
    vel[i] = np.sqrt(U[i]**2 + V[i]**2)
        
x = np.linspace(0,1,nx)
y = np.linspace(0,ny/float(nx),ny)
xx, yy = np.meshgrid(x,y)

for unit in range(nsteps):
    # plt.figure()
    # plot2 = plt.contourf(xx,yy,U[unit],np.arange(-1,1,0.01), extend='both')
    # plt.colorbar(plot2)
    # plt.axis([0, 1, 0, ny/float(nx)])
    # plt.savefig('contour-U' + str(unit) + '.png')
    # plt.close()
    #
    plt.figure()
    plot2 = plt.contourf(xx,yy,V[unit],np.arange(-1,1,0.01), extend='both')
    plt.colorbar(plot2)
    plt.axis([0, 1, 0, ny/float(nx)])
    plt.savefig('contour-V' + str(unit) + '.png')
    plt.close()
    #
    plt.figure()
    plot2 = plt.contourf(xx,yy,P[unit],np.arange(0,2,0.01), extend='both')
    plt.colorbar(plot2)
    plt.axis([0, 1, 0, ny/float(nx)])
    plt.savefig('contour-P' + str(unit) + '.png')
    plt.close()
    #
    # plt.figure()
    # plot2 = plt.contourf(xx,yy,T[unit],50)
    # plt.colorbar(plot2)
    # plt.axis([0, 1, 0, ny/float(nx)])
    # plt.savefig('contour-T' + str(unit) + '.png')
    # plt.close()
    #
    # plt.figure()
    # plot2 = plt.contourf(xx,yy,rho[unit],200)
    # plt.colorbar(plot2)
    # plt.axis([0, 1, 0, ny/float(nx)])
    # plt.savefig('contour-rho' + str(unit) + '.png')
    # plt.close()
    #
    # plt.figure()
    # plot1 = plt.streamplot(xx,yy,U[unit],V[unit],color=vel[unit])
    # plt.colorbar()
    # plt.axis([0, 1, 0, ny/float(nx)])
    # plt.savefig('streamplot' + str(unit) + '.png')
    # plt.close()
    print unit
