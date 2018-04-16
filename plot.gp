
w = 8.5
h = 5.0
set terminal epslatex size w cm,h cm color dashlength 0.3
set output 'plot.tex'
borders_width = 0.5
set termopt enhanced 

set linestyle 1 lc rgb '#000000'  lw 2   dt 1
set linestyle 2 lc rgb '#660000'  lw 3.5   dt 1
set linestyle 3 lc rgb '#990000'  lw 3   dt 1
set linestyle 4 lc rgb '#aa0000'  lw 2.5   dt 1
set linestyle 5 lc rgb '#cc0000'  lw 2   dt 1

set xrange [0.:0.4]
set xtics 0.1
set xlabel 'hole doping $p$'
set ylabel 'density of states $N(0)$' offset -1,0
#set yrange [0:7]
#set ytics 2
set tics scale 0.7
unset key

#set label 2 '\scriptsize $\tau=250$' at first 0.24,6

fac = 2*3.1416*3.1416/3

plot 'transport_vs_mu.dat' u 'p0':'dos' w l ls 1
