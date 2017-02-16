#!/usr/bin/gnuplot


analytical(x) = ( ( (100/0.02) + (1E6/(2*0.5))*(0.02-x)) *x ) + 100
plot analytical(x) title "Analytical solution" w lines, "output-no-source.dat" title "FVM solution"
set xrange [0:0.02]
set terminal png  enhanced
set output "conduction01-no-source.png"
set xlabel 'X [m]'
set ylabel 'Temperature [K]'
set title 'heat conduction in a rod'
replot
