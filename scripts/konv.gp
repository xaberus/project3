set terminal pdfcairo solid size 13cm,8cm
set output "harmosc_konv.pdf"

set grid x y2
set key center top title " "

set xlabel "Zeitschritte * 10000"
set ylabel "Abweichung"

plot  \
  "harmosc.prog.dat" using ($0+30):1 with lines,\
  "harmosc.prog.dat" using ($0+30):2 with lines,\
  "harmosc.prog.dat" using ($0+30):3 with lines,\
  0;


reset

set terminal pdfcairo solid size 13cm,8cm
set output "simple_konv.pdf"

set grid x y2
set key center top title " "

set xlabel "Zeitschritte * 10000"
set ylabel "Abweichung"

plot  \
  "simple.prog.dat" using ($0+30):1 with lines,\
  "simple.prog.dat" using ($0+30):2 with lines,\
  "simple.prog.dat" using ($0+30):3 with lines,\
  0;


reset

set terminal pdfcairo solid size 13cm,8cm
set output "square_konv.pdf"

set grid x y2
set key center top title " "

set xlabel "Zeitschritte * 10000"
set ylabel "Abweichung"

plot  \
  "square.prog.dat" using ($0+30):1 with lines,\
  "square.prog.dat" using ($0+30):2 with lines,\
  "square.prog.dat" using ($0+30):3 with lines,\
  0;

