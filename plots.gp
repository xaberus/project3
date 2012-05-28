#!/usr/bin/env gnuplot

set term pngcairo size 1500,500;
set output "plotck.png"
plot  "plotck.txt" using 1:4 w l, "plotE.txt" u 1:(6) w i


set term pngcairo size 1500,500;
set output "plotc.png"
plot [x=0:4] "plotc.txt" using 1:2 w lp