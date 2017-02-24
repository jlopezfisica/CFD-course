splot 'output.dat' u 1:2:6 w lines
set xlabel 'x[m]'
set ylabel 'y[m]'
set zlabel 'T[K]'
set terminal png enhanced
set output 'temperature.png'
replot
