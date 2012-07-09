
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

Emin=`cut -d: -f 15 stats`
Emax=`cut -d: -f 16 stats`

set terminal pdfcairo enhanced solid size 13cm,8cm
set output "relen.pdf"

set xlabel "Peak"
set ylabel "Abweichung"

plot "relen.dat" using 0:(($2-$1)/$2) with lp title "(ccsen-numen)/numen";

