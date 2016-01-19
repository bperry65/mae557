import numpy as np

files = ['test-re10_central_nx8_dt5.00E-06_tend1.00E-01',
         'test-re10_central_nx16_dt5.00E-06_tend1.00E-01',
         'test-re10_central_nx32_dt5.00E-06_tend1.00E-01',
         'test-re10_central_nx64_dt5.00E-06_tend1.00E-01',
         'test-re10_central_nx128_dt5.00E-06_tend1.00E-01',
         'test-re10_central_nx256_dt5.00E-06_tend1.00E-01']

# files = ['test_central_nx8_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx16_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx32_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx64_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx128_dt2.50E-08_tend2.00E-03',
#          'test_compressible_nx256_dt2.50E-08_tend2.00E-03']

steps = [ 8,
          16,
          32,
          64,
          128,
          256]
steps2 = [ 16,
          30,
          59,
          117,
          232,
           462]
steplength = [0]*len(steps)

for i in range (0,len(steps)):
    steplength[i] = 1.0/float(steps[i])

nofiles = len(files)
x = [0] * nofiles
error_p = np.zeros(nofiles-1)
error_u = np.zeros(nofiles-1)
error_v = np.zeros(nofiles-1)
error_T = np.zeros(nofiles-1)

for i in range(0,nofiles):
    x[i] = np.genfromtxt(files[i],skip_header=4,usecols=(0,1,2,3,5,6))

for i in range(0,nofiles-1):
    #nx = len(x[i])
    for m in range(0,7):
        for k in range(0,7):
            j1 = steps[i]*m/8 * steps2[i] + steps[i]*k/8
            j2 = steps[i+1]*m/8 * steps2[i+1] + steps[i+1]*k/8
            #print 'x', x[i][j1,1], x[i+1][j2,1]
            #print 'y', x[i][j1,2], x[i+1][j2,2]
            error_p[i] = error_p[i] + abs(x[i][j1,0] - x[i+1][j2,0])
            error_u[i] = error_u[i] + abs(x[i][j1,1] - x[i+1][j2,1])
            error_v[i] = error_v[i] + abs(x[i][j1,2] - x[i+1][j2,2])
            error_T[i] = error_T[i] + abs(x[i][j1,3] - x[i+1][j2,3])
    print steplength[i], error_p[i], error_u[i], error_v[i], error_T[i]
