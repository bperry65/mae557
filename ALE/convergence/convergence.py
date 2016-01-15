import numpy as np
files = ['test_compressible_nx40_dt1.02E-04_tend2.05E-02',
         'test_compressible_nx40_dt5.12E-05_tend2.05E-02',
         'test_compressible_nx40_dt2.56E-05_tend2.05E-02',
         'test_compressible_nx40_dt1.28E-05_tend2.05E-02',
         'test_compressible_nx40_dt6.40E-06_tend2.05E-02',
         'test_compressible_nx40_dt3.20E-06_tend2.05E-02',
         'test_compressible_nx40_dt1.60E-06_tend2.05E-02',
         'test_compressible_nx40_dt8.00E-07_tend2.05E-02',
         'test_compressible_nx40_dt4.00E-07_tend2.05E-02',
         'test_compressible_nx40_dt2.00E-07_tend2.05E-02',
         'test_compressible_nx40_dt1.00E-07_tend2.05E-02']
steps = ( 1.024e-4,
          0.512e-4,
          0.256e-4,
          0.128e-4,
          0.064e-4,
          0.032e-4,
          0.016e-4,
          0.008e-4,
          0.004e-4,
          0.002e-4,
          0.001e-4)

nofiles = len(files)

x = [0] * nofiles
error = np.zeros(nofiles-1)

for i in range(0,nofiles):
    x[i] = np.genfromtxt(files[i],skip_header=3,usecols=(1,2))
    
nx = len(x[nofiles-1])

for i in range(0,nofiles-1):
    for j in range(0,nx):
        error[i] = error[i] + abs(x[i][j,1] - x[i+1][j,1])
error = error/nx
logerror = np.log(error)
logsteps = np.log(steps[0:nofiles-1])
for i in range(0,nofiles-1):
    print steps[i],error[i]

#for i in range(0,nofiles-1):
#    print logsteps[i],logerror[i]
    
         
