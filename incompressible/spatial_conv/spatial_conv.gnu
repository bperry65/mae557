set terminal postscript eps enhanced color

set size 0.7,0.7

set output 'spaceconv.eps'
set logscale xy

set xlabel '{/Symbol D}x'
set ylabel 'error = |u_{{/Symbol D}x} - u_{{/Symbol D}x/2}|

set key at graph 0.9,0.25
set arrow from 0.0025,3e-5 to 0.02,2.4e-4 nohead
set label at graph 0.12,0.68 '1st order'
set arrow from 0.0025,2e-7 to 0.02,12.8e-6 nohead
set label at graph 0.12,0.4 '2nd order'

p 'upwind_euler' u 1:2 w l lt 1 lw 2.5 lc 1 t 'Upwind (euler)', \
  'upwind_trapezoid' u 1:2 w l lt 3 lw 1.5 lc 3 t 'Upwind (trapezoid)', \
  'central_euler' u 1:2 w l lt 1 lw 2.5 lc 2  t 'Central (euler)', \
  'central_trapezoid' u 1:2 w l lt 3 lw 1.5 lc 4 t 'Central (trapezoid)'