#set terminal qt size 1350,550
set key top center
set xlabel '|x-y|' font ',18'
set ylabel 'G(|x-y|)' font ',18' offset -2,0
set title 'Correlation Function' font ',18'
set xtics font ',14'
set ytics font ',14'
set grid x,y
set lmargin 14
set bmargin 4
set key font ',14'
set key outside 
set key right
#set key box vertical width -1 height 1
#set logscale y
#set yrange [0.001:2]
#set xrange [0:32]


nn=51
L=32
mi=2.0
mf=-5.0

array column2[nn]
array column3[nn]
array at[nn]
array mrt[nn]
array aterr[nn]
array mrterr[nn]
array ct[nn]
array chi2[nn]
array mu02[nn]

do for [i=1:nn] {
#set style line i pt 1
column2[i]=1+i
column3[i]=1+nn+i
mu02[i]=mi+(mf-mi)*(i-1)/(nn-1)
#mu02[i]=-1.265+(-1.285+1.265)*(i-1)/(nn-1)
}

set style line 1 lc rgb "#FF0000" pt 1 # Red
set style line 2 lc rgb "#00FF00" pt 1 # Green
set style line 3 lc rgb "#0000FF" pt 1 # Blue
set style line 4 lc rgb "#FFFF00" pt 1 # Yellow
set style line 5 lc rgb "#FF00FF" pt 1 # Magenta
set style line 6 lc rgb "#00FFFF" pt 1 # Cyan
set style line 7 lc rgb "#FFA500" pt 1 # Orange
set style line 8 lc rgb "#800080" pt 1 # Purple
set style line 9 lc rgb "#008000" pt 1 # Dark Green
set style line 10 lc rgb "#8B0000" pt 1 # Dark Red
set style line 11 lc rgb "#4682B4" pt 1 # Steel Blue
set style line 12 lc rgb "#D2B48C" pt 1 # Tan
set style line 13 lc rgb "#F0E68C" pt 1 # Khaki
set style line 14 lc rgb "#9ACD32" pt 1 # Yellow Green
set style line 15 lc rgb "#BA55D3" pt 1 # Medium Orchid
set style line 16 lc rgb "#FF6347" pt 1 # Tomato

f(x,a,mr,c)=a*cosh((x-L/2)/mr)+c
#f(x,a,mr)=a*(exp(-mr*x)+exp(-mr*(L-x)))/(2.0*mr)
#f(x,a,mr)=a*cosh(mr*(x-L/2))
mr=2.1
a=0.1
#c=0.0001

do for [i=1:nn] {
#mr=0.8
#a=0.8
#c=0.0001
#fit f(x,a,mr) '../data/corrfunc.dat' u 1:column2[i]:column3[i] every ::0::(L-1) via a,mr
fit f(x,a,mr,c) '../data/corrfunc.dat' u 1:column2[i]:column3[i] every ::0::(L-2) via a,mr,c
at[i]=a
aterr[i]=a_err
mrt[i]=mr
mrterr[i]=mr_err
ct[i]=c
chi2[i]=(FIT_STDFIT*FIT_STDFIT)
}

set print "results.dat"
do for [i=1:nn] {
print mu02[i], mrt[i], mrterr[i]
}
set print

#plot for [i=1:nn] '../data/corrfunc.dat' u 1:column2[i]:column3[i] w errorbars notitle linestyle i, for [i=1:nn] f(x,at[i],mrt[i]) title sprintf(' λ_0=%.2f, m_r=%.4f±%.4f, χ^2/dof=%.2f',mu02[i],mrt[i],mrterr[i],chi2[i]) linestyle i lw 2
plot for [i=1:nn] '../data/corrfunc.dat' u 1:column2[i]:column3[i] w errorbars notitle linestyle i lw 2, for [i=1:nn] f(x,at[i],mrt[i],ct[i]) linestyle i lw 2 title sprintf(' μ^2=%.2f, ξ=%.3f(%.0f), χ^2/dof=%.2f',mu02[i],mrt[i],1000*mrterr[i],chi2[i])


pause -1
