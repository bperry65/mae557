set terminal postscript eps enhanced color

set size 0.7,0.7

set output 'timeconv.eps'
set logscale xy

set xlabel '{/Symbol D}t'
set ylabel 'error = |{/Symbol f}_{{/Symbol D}t} - {/Symbol f}_{{/Symbol D}t/2}|
set format y "%4.0e"

set key at graph 0.93,0.3
set arrow from 1e-6,5e-6 to 1e-5,5e-5 nohead
set label at graph 0.16,0.73 '1st order'
#set arrow from 2e-5,1.5e-10 to 4e-4,6e-8 nohead
#set label at graph 0.2,0.37 '2nd order'

p 'results' u 1:2 w l lt 1 lw 2.5 lc 1 t 'Pressure', \
  'results' u 1:3 w l lt 1 lw 2.5 lc 3 t 'Velocity - u', \
  'results' u 1:4 w l lt 1 lw 2.5 lc 2 t 'Velocity - v', \
  'results' u 1:5 w l lt 1 lw 2.5 lc 4 t 'Temperature', \