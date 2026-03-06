set title 'Histograms for L=64' font ',22'
set tics font ',18'
set xtics font ',16'
set ytics font ',16'
set key font ',18'
set lmargin 18
set bmargin 5
set xlabel 'ψ' font ',20' offset 0,-1
set ylabel 'Freq.' font ',20' offset -5,0
set grid x,y

p 'histdt1p0.dat' w errorbars title 'a=1.0', 'histdt0p8.dat' w errorbars title 'a=0.8', 'histdt0p6.dat'  w errorbars title 'a=0.6', 'histdt0p4.dat' w errorbars title 'a=0.4'
pause -1
