set title 'Histograma para L=64' font ',22'
set tics font ',18'
set xtics font ',16'
set ytics font ',16'
set key font ',18'
set lmargin 18
set bmargin 5
set xlabel 'x' font ',20' offset 0,-1
set ylabel '|ψ|^2' font ',20' offset -5,0 rotate by 0
set grid x,y

p 'hist1p0.dat' w errorbars notitle lw 2
pause -1
