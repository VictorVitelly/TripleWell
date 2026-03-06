set title 'Ajuste para L=256 (oscilador armónico)' font ',20'
set tics font ',18'
set xtics font ',16'
set ytics font ',16'
set key font ',14'
set lmargin 18
set bmargin 5
set xlabel 'x' font ',20' offset 0,-1
set ylabel '|ψ|^2' font ',20' offset -5,0 rotate by 0
set grid x,y
set yrange [0:0.7]

f(x)=a*exp(-b*x**2)

fit f(x) 'histogram.dat' u 1:2:3 via a,b
chi2=(FIT_STDFIT*FIT_STDFIT)

p 'histogram.dat' w errorbars notitle lw 2, f(x) title sprintf('k_1=%.4f(%.0f), k_2=%.4f(%.0f), χ^2/dof=%.2f',a,10000*a_err,b,10000*b_err,chi2) lw 2 lc rgb "red"
pause -1
