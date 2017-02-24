plot 'CDS.dat', 'UPW.dat', 'HYB.dat'
set xlabel 'Number of grid points'
set ylabel 'Temperature[K]'
set terminal png enhanced
set title 'Convergence of schemes'
set output 'plot-schemes.png'
replot
