set terminal postscript eps enhanced color

set size 0.7,0.7

set output 'spaceconv.eps'
set logscale xy

set xlabel '{/Symbol D}x'
set ylabel 'error = |{/Symbol f}_{{/Symbol D}x} - {/Symbol f}_{{/Symbol D}x/2}|
set xrange [0.005:0.5]

set key at graph 0.9,0.25
set arrow from 0.1,0.005 to 0.3,0.015 nohead
set label at graph 0.65,0.62 '1st order'
set arrow from 0.01,0.0001 to 0.03,0.0009 nohead
set label at graph 0.12,0.33 '2nd order'

p 'results' u 1:2 w l lt 1 lw 2.5 lc 1 t 'Pressure', \
  'results' u 1:3 w l lt 1 lw 2.5 lc 3 t 'Velocity - u', \
  'results' u 1:4 w l lt 1 lw 2.5 lc 2 t 'Velocity - v', \
  'results' u 1:5 w l lt 1 lw 2.5 lc 4 t 'Temperature', \