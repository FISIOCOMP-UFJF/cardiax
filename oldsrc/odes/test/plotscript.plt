#
# gnuplot script to visualize the action potential data
# 
set terminal postscript
set output "ap.eps"
set xlabel "ms"
set ylabel "mV"
plot 'ap' using 1:2 title 'transmembrane voltage' with lines 
print "\nPostscript file written to ap.eps file.\n"
