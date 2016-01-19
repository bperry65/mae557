import numpy as np
files = ['test_central_nx40_dt2.05E-04_tend2.05E-02',
         'test_central_nx40_dt1.02E-04_tend2.05E-02',
         'test_central_nx40_dt5.12E-05_tend2.05E-02',
         'test_central_nx40_dt2.56E-05_tend2.05E-02',
         'test_central_nx40_dt1.28E-05_tend2.05E-02',
         'test_central_nx40_dt6.40E-06_tend2.05E-02',
         'test_central_nx40_dt3.20E-06_tend2.05E-02',
         'test_central_nx40_dt1.60E-06_tend2.05E-02',
         'test_central_nx40_dt8.00E-07_tend2.05E-02',
         'test_central_nx40_dt4.00E-07_tend2.05E-02',
         'test_central_nx40_dt2.00E-07_tend2.05E-02',
         'test_central_nx40_dt1.00E-07_tend2.05E-02']
steps = ( 2.048e-4,
          1.024e-4,
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
error_p = np.zeros(nofiles-1)
error_u = np.zeros(nofiles-1)
error_v = np.zeros(nofiles-1)
error_T = np.zeros(nofiles-1)

for i in range(0,nofiles):
    x[i] = np.genfromtxt(files[i],skip_header=4,usecols=(0,1,2,3))
    
nx = len(x[nofiles-1])

for i in range(0,nofiles-1):
    for j in range(0,nx):
        error_p[i] = error_p[i] + abs(x[i][j,0] - x[i+1][j,0])
        error_u[i] = error_u[i] + abs(x[i][j,1] - x[i+1][j,1])
        error_v[i] = error_v[i] + abs(x[i][j,2] - x[i+1][j,2])
        error_T[i] = error_T[i] + abs(x[i][j,3] - x[i+1][j,3])

error_p = error_p/nx
error_u = error_u/nx
error_v = error_v/nx
error_T = error_T/nx

for i in range(0,nofiles-1):
    print steps[i], error_p[i], error_u[i], error_v[i], error_T[i]

#for i in range(0,nofiles-1):
#    print logsteps[i],logerror[i]
    
         
