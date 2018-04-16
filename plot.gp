w = 10.5
h = 8
set terminal epslatex size w cm,h cm color dashlength 0.3
set output 'plot.tex'
borders_width = 0.5
set termopt enhanced 

set linestyle 1 lc rgb 'blue'	 lw 3 dt 1 ps 0.6 pt 7
set linestyle 2 lc rgb 'purple'	 lw 3 dt 1 ps 0.6 pt 7
set linestyle 3 lc rgb '#009900'	 lw 3 dt 1 ps 0.6 pt 7
set linestyle 4 lc rgb 'green' lw 3 dt 1 ps 0.6 pt 7
set linestyle 5 lc rgb 'red'	 lw 3 dt 1 ps 0.6 pt 7
set linestyle 6 lc rgb 'orange'	 lw 3 dt 1 ps 0.6 pt 7

set title "$t=0.190$ eV, $\\hbar/\\tau=t/25$, $t_z=-0.07t$"


x="(column('p0'))"

set xlabel '$p$'
set xrange [-1.1:1.1]
set xtics .2


y1="(column('sigmaxx')/10)"
y2="(column('sigmaxy')/200)"
y3="(10*column('alphaxx'))"
y4="(column('alphaxy')/10)"
y5="(column('betaxx')/40)"
y6="(column('betaxy')/400)"

set ylabel '$\sigma$, $\alpha$, $\beta$' offset -1,0
set yrange [-2:4]
set ytics 1


# set key at graph 0.7,0.95 top left samplen 1 title '$T\approx$'
set key top left samplen 2 spacing 0.8 maxrows 2
# set label 1 '\scriptsize $\tau=25$' at first 0.21,4.3
# set label 2 '\scriptsize $\tau=250$' at first 0.24,6

set palette defined (0 'blue', 1 'orange')
unset colorbox

N = 1

plot 0 ls 1 lc 'black' lw 1 notitle,\
'transport_vs_mu.dat' i 0 u @x:@y1 w lp ls 1 t "$\\sigma_{xx}/10$"	,\
'transport_vs_mu.dat' i 0 u @x:@y2 w lp ls 2 t "$\\sigma_{xy}/200$"	,\
'transport_vs_mu.dat' i 0 u @x:@y3 w lp ls 3 t "$10\\alpha_{xx}$"	,\
'transport_vs_mu.dat' i 0 u @x:@y4 w lp ls 4 t "$\\alpha_{xy}/10$"	,\
'transport_vs_mu.dat' i 0 u @x:@y5 w lp ls 5 t "$\\beta_{xx}/40$"	,\
'transport_vs_mu.dat' i 0 u @x:@y6 w lp ls 6 t "$\\beta_{xy}/400$" ,\
'transport_vs_mu.dat' i 0 u @x:(column('dos')) w l lw 3 lc 'black' t ""

