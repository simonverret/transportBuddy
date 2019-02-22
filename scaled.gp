w = 800
h = 600
set term wxt font "Arial,16" persist size w,h
bind "w" "unset output; exit gnuplot"

set title '' #"$t=0.190$ eV, $\\hbar/\\tau=t/25$, $t_z=-0.07t$"

x="(column('p'))"
set xlabel 'p0'
# set xrange [-1:1]
# set xtics .2

sfac=1.845e-3
# sfac=1

sxx="(column('sigmaxx')*column('eta')*2)"
rxx="(1/(column('sigmaxx')*column('eta')*2))"
sxy="(column('sigmaxy')*column('eta')*column('eta')*4)"
axx="(column('alphaxx')*column('eta')*2/column('T'))"
axy="(column('alphaxy')*column('eta')*column('eta')*4/column('T'))"
bxx="(column('betaxx')*column('eta')*2/pi)"
bxy="(column('betaxy')*column('eta')*column('eta')*4/pi)"
cv="(column('Cv')/column('T')/pi/pi*3)"
dos="(column('dos'))"
set ylabel 's, a, b, cv, dos, doc' offset -1,0
# set yrange [-0.3:0.6]
# set ytics 1
set y2tics scale 0.5

set key top right Left samplen 2 spacing 0.8 reverse
set palette defined (0 'blue', 1 'orange')
unset colorbox


set linestyle 1  lc rgb 'blue'       lw 4 dt 1 ps 0.6 pt 7
set linestyle 2  lc rgb 'purple'     lw 4 dt 1 ps 0.6 pt 7
set linestyle 3  lc rgb '#009900'    lw 3 dt 1 ps 0.6 pt 7
set linestyle 4  lc rgb 'green'      lw 3 dt 1 ps 0.6 pt 7
set linestyle 5  lc rgb 'red'        lw 3 dt 1 ps 0.6 pt 7
set linestyle 6  lc rgb 'orange'     lw 3 dt 1 ps 0.6 pt 7
set linestyle 7  lc rgb 'black'      lw 4 dt 1 ps 0.6 pt 7
set linestyle 9  lc rgb '#555555'       lw 3 dt 1 ps 0.6 pt 7
set linestyle 10 lc rgb '#888888'       lw 3 dt 1 ps 0.6 pt 7
set linestyle 11 lc rgb '#aaaaaa'       lw 3 dt 1 ps 0.6 pt 7

plot 0 ls 1 lc 'black' lw 1 notitle,\
'transport_vs_mu.dat' i 0 u @x:@sxx w l ls 1 t "sxx/tau" ,\
'transport_vs_mu.dat' i 0 u @x:@sxy w l ls 2 t "sxy/tau^2"    ,\
'transport_vs_mu.dat' i 0 u @x:@axx w l ls 3 t "axx/T*tau"    ,\
'transport_vs_mu.dat' i 0 u @x:@axy w l ls 4 t "axy/T*tau^2"    ,\
'transport_vs_mu.dat' i 0 u @x:@bxx w l ls 5 t "bxx/pi*tau"    ,\
'transport_vs_mu.dat' i 0 u @x:@bxy w l ls 6 t "bxy/pi*tau^2" ,\
'transport_vs_mu.dat' i 0 u @x:@cv  w l ls 7 t "3cv/pi^2",\
'transport_vs_mu.dat' i 0 u @x:@dos w l ls 8 t "dos"

pause -1