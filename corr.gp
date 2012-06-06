bins=`cut -d: -f 1 stats`
dt=`cut -d: -f 2 stats`
rangemin=`cut -d: -f 3 stats`
rangemax=`cut -d: -f 4 stats`
steps=`cut -d: -f 5 stats`
runs=`cut -d: -f 6 stats`
Vmin=`cut -d: -f 7 stats`
Vmax=`cut -d: -f 8 stats`
apsimin=`cut -d: -f 9 stats`
apsimax=`cut -d: -f 10 stats`
tsteps=`cut -d: -f 11 stats`
dx=`cut -d: -f 12 stats`
dk=`cut -d: -f 13 stats`
dE=`cut -d: -f 14 stats`

set terminal pdfcairo solid size 13cm,8cm
set output "corr.pdf"

set grid x y2
set key center top title " "

#set title "Autokorrelation"

set xrange [0 : dt*steps*runs/5]

set xlabel "Zeit"

set datafile separator ";"

plot \
  "corr.dat" axis x1y1 with lines linestyle 3 title "c(t)"

