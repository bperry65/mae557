import numpy as np
files = ['test_centered_nx40_dt1.60E-06_tend2.00E-01',
         'test_centered_nx40_dt8.00E-07_tend2.00E-01',
         'test_centered_nx40_dt4.00E-07_tend2.00E-01',
         'test_centered_nx40_dt2.00E-07_tend2.00E-01',
         'test_centered_nx40_dt1.00E-07_tend2.00E-01',
         'test_centered_nx40_dt5.00E-08_tend2.00E-01',
         'test_centered_nx40_dt2.50E-08_tend2.00E-01']
steps = ( 0.01600e-4,
          0.00800e-4,
          0.00400e-4,
          0.00200e-4,
          0.00100e-4,
          0.00050e-4,
          0.00025e-4)

nofiles = len(files)

x = [0] * nofiles
error = np.zeros(nofiles-1)

for i in range(0,nofiles):
    x[i] = np.genfromtxt(files[i],skip_header=3,usecols=(1,2))
    
nx = len(x[nofiles-1])


for i in range(0,nofiles-1): #nofiles-
#    print x[i]
#    for j in range(0,nx):
    j = 237
    error[i] = error[i] + abs(x[i][j,1] - x[i+1][j,1])
#    print abs(x[i][j,1] - x[i+1][j,1])
#    print error[i]
error = error/nx
logerror = np.log(error)
logsteps = np.log(steps[0:nofiles-1])
for i in range(0,nofiles-1):
    print steps[i],error[i]

#for i in range(0,nofiles-1):
#    print logsteps[i],logerror[i]
    
         
