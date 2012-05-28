#!/usr/bin/env gnuplot

set term pngcairo size 1500,500;
set output "plotck.png"
plot  "plotck.txt" using 1:(sqrt($2**2+$3**2)) w l

set term pngcairo size 1500,500;
set output "plotc.png"
plot [x=0:4] "plotc.txt" using 1:2 w lp