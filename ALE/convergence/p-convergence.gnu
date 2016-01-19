set terminal postscript eps enhanced color

set size 0.7,0.7

set output 'timeconv.eps'
set logscale xy

set xlabel '{/Symbol D}t'
set ylabel 'error = |u_{{/Symbol D}t} - u_{{/Symbol D}t/2}|
set format y "%4.0e"

set key at graph 0.93,0.2
set arrow from 1e-7,1e-10 to 1e-6,1e-9 nohead
set label at graph 0.30,0.90 '1st order'
#set arrow from 2e-5,1.5e-10 to 4e-4,6e-8 nohead
#set label at graph 0.2,0.37 '2nd order'

p 'convint' u 1:2 w l lt 1 lw 2 lc 3 t 'Compressible'
