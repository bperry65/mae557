set terminal postscript eps enhanced color

set size 0.7,0.7

set output 'spaceconv.eps'
set logscale xy

set xlabel '{/Symbol D}x'
set ylabel 'error = |u_{{/Symbol D}x} - u_{{/Symbol D}x/2}|
set xrange [0.005:0.5]

set key at graph 0.9,0.25
set arrow from 0.01,0.01 to 0.03,0.03 nohead
set label at graph 0.12,0.6 '1st order'
set arrow from 0.01,0.0003 to 0.03,0.0027 nohead
set label at graph 0.12,0.3 '2nd order'

p '../../incompressible/spatial_conv/conv-incompressible' u 1:2 w l lt 1 lw 2.5 lc 1 t 'Incompressible', \
  'conv-compressible' u 1:2 w l lt 3 lw 1.5 lc 3 t 'Compressible'