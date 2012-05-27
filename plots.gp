#!/usr/bin/env gnuplot

set term pngcairo size 1500,500;
set output "plotck.png"
plot  "plotck.txt" using 0:(sqrt($1**2+$2**2)) w l