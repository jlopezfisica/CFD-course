splot 'output-CDS.dat' w line
set terminal png  enhanced
set output "plot-CDS-2x2.png"
set xlabel 'X [m]'
set ylabel 'Y [m]'
set zlabel 'T [K]'
set title '2D steady convection diffusion'
set pm3d
replot


reset
splot 'output-UPW.dat' w line
set terminal png  enhanced
set output "plot-UPW-2x2.png"
set xlabel 'X [m]'
set ylabel 'Y [m]'
set zlabel 'T [K]'
set title '2D steady convection diffusion'
set pm3d
replot


reset
splot 'output-HYB.dat' w line
set terminal png  enhanced
set output "plot-HYB-2x2.png"
set xlabel 'X [m]'
set ylabel 'Y [m]'
set zlabel 'T [K]'
set title '2D steady convection diffusion'
set pm3d
replot
