#!/usr/bin/gnuplot


analytical(x) = (800*x) + 100
plot analytical(x) title "Analytical solution" w lines, "output.dat" title "FVM solution"
set xrange [0:0.5]
set terminal png  enhanced
set output "conduction01.png"
set xlabel 'X [m]'
set ylabel 'Temperature [K]'
set title 'heat conduction in a rod'
replot
