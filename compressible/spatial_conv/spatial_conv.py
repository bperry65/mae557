import numpy as np

files = ['test_compressible_nx8_dt5.00E-07_tend2.00E-03',
         'test_compressible_nx16_dt5.00E-07_tend2.00E-03',
         'test_compressible_nx32_dt5.00E-07_tend2.00E-03',
         'test_compressible_nx64_dt5.00E-07_tend2.00E-03',
         'test_compressible_nx128_dt5.00E-07_tend2.00E-03']

# files = ['test_compressible_nx8_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx16_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx32_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx64_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx128_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx256_dt2.50E-08_tend2.00E-03']

steps = [ 8,
          16,
          32,
          64,
          128]
steplength = [0]*len(steps)

for i in range (0,len(steps)):
    steplength[i] = 1.0/float(steps[i])

nofiles = len(files)
x = [0] * nofiles
error = np.zeros(nofiles-1)

for i in range(0,nofiles):
    x[i] = np.genfromtxt(files[i],skip_header=3,usecols=(1,2))

for i in range(0,nofiles-1):
    nx = len(x[i])
    j1 = steps[i]*6/8 * (steps[i]+1) + steps[i]/2 +1
    j2 = steps[i+1]*6/8 * (steps[i+1]+1) + steps[i+1]/2 + 1
    error[i] = error[i] + abs(x[i][j1,1] - x[i+1][j2,1])
    for k in range(1,4):
        error[i] = error[i] + abs(x[i][j1+k*steps[i]/8,1] - x[i+1][j2+k*steps[i+1]/8,1])
        error[i] = error[i] + abs(x[i][j1-k*steps[i]/8,1] - x[i+1][j2-k*steps[i+1]/8,1])
        
    print steplength[i], error[i]
