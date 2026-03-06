# Gnuplot script: normalized histogram

# Input file
datafile = "myhistogram.dat"

# Histogram parameters
binwidth = 6./200.
bin(x,width) = width*floor(x/width) + width/2.0

# Count number of data points
stats datafile using 1 nooutput
N = STATS_records

# Plot normalized histogram
set boxwidth binwidth
set style fill solid 0.7
set xlabel "x"
set ylabel "Probability density"
set title "Normalized Histogram"

plot datafile using (bin($1,binwidth)):(1.0/(N*binwidth)) smooth freq with boxes
pause -1
