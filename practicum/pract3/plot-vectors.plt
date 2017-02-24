plot 'output.dat' w vector
set xlabel 'x[m]'
set ylabel 'y[m]'
set title 'velocity'
set terminal png enhanced
set output 'velocity.png'
replot
