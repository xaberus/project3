
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

set terminal pdfcairo solid size 13cm,8cm
set output "spektrum.pdf"

set grid x y2
set key center top title " "

set title "Autokorrelation"

set xlabel "Energie"

set xrange [Emin : Emax]

plot \
  "dftcorr.dat" using 1:4 with lines title "|F[<psi|U(t)|psi>](E)|", \
  "theoenrg.dat" using 2:(-1) with impulses title "theoretische Energien";

