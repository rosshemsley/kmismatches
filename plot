# use from prompt: gnuplot plotTTT
# This file plots columns of data from a file called "data".
# There must be 4 columns:
# 1 - index (number of games played, usually in 1000s)
# 2 - wins by x
# 3 - wins by o
# 4 - draws
# Press return to return to the command prompt (you may
# have to click on the prompt window first).
#
set title "k-Mismatches case 2"
set style data linespoints
#set yrange[0:100]
set autoscale
set ylabel "Time Taken (Seconds)"
set xlabel "Size of Pattern (Bytes)"
set key left top

# plot './data' using 1:2 t "Wins with Bolzmann 0.001", 'data' using 1:3 t "Wins by Greedy 0.1", 'data' using 1:4 t "Draws"

set term postscript enhanced color 
set output 'out.eps'
set size 0.7, 0.5 

# plot './data' using 1:4 t "Wins by X", 'data' using 1:4 t "Wins by O", 'data' using 1:4 t "Draws"

plot './stats' using 1:2 t  "Naive", './stats' using 1:3 t "Kangarooing"


