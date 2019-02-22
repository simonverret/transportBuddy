w = 800
h = 600
set term wxt font "Arial,16" persist size w,h
bind "w" "unset output; exit gnuplot"

set title ''

echarge = 1.602e-19
Tfac = 2204.8
sfac = 1.845e5   # (Ohm * m)^-1
afac = 8.617e-5*sfac  # A * (K * m)^-1
nHfac = 1/sfac
CvFac = 2.866 # mJ (K^2 mol)-1 
RHfac = 0.0002137  # m^3 / C

T  = "(column('T')*Tfac)"
mu = "(column('mu'))"
p  = "(column('p0'))"

sxx     = "(sfac*column('sigmaxx'))"            # (Ohm * m)^-1
rxx     = '(10e7/'.sxx.')'                      #  muOhm * cmm
sxy     = "(sfac*column('sigmaxy'))"            # (Ohm * m)^-1   / B []
rxy     =  '('.sxy.'/('.sxx.'*'.sxx.'))'       # Ohm * m
nH      =  '(nHfac/'.rxy.')'                    # 1+p
RH      =  '(1e9*RHfac*'.rxy.')'                # mm^3 / C
axx     = "(afac*column('alphaxx'))"            # A * (Ohm * m)^-1
Seebeck = '(-10e6*'.axx.'/'.sxx.'/'.T.')'       # muV / K^2

axy     = "(column('alphaxy'))"
bxx     = "(column('betaxx'))"
bxy     = "(column('betaxy'))"
Cv      = "(CvFac*column('Cv')/column('T'))"    # mJ (K^2 mol)-1 
dos     = "(column('dos'))"                     # at E_F

y = Seebeck

set ylabel 'see file' offset -1,0
# set yrange [-0.3:0.6]
# set ytics 1

x = p
set xlabel 'p'
# set xrange [-1:1]
# set xtics .2
# set xrange [0.0004:0.12]
# set logscale x

set key top right Left samplen 2 spacing 0.8 reverse maxrows 10
set palette defined (0 'blue', 1 'orange')
unset colorbox

set linestyle 1  lc rgb 'blue'       lw 4 dt 1 ps 0.6 pt 7
set linestyle 2  lc rgb 'purple'     lw 4 dt 1 ps 0.6 pt 7
set linestyle 3  lc rgb '#009900'    lw 3 dt 1 ps 0.6 pt 7
set linestyle 4  lc rgb 'green'      lw 3 dt 1 ps 0.6 pt 7
set linestyle 5  lc rgb 'red'        lw 3 dt 1 ps 0.6 pt 7
set linestyle 6  lc rgb 'orange'     lw 3 dt 1 ps 0.6 pt 7
set linestyle 7  lc rgb 'black'      lw 4 dt 1 ps 0.6 pt 7
set linestyle 8  lc rgb '#888888'    lw 3 dt 1 ps 0.6 pt 7

set linestyle 9                      lw 2 dt 1 ps 0.3 pt 7
# set palette defined (0 'red',0.5 'black' , 1 'green')
set palette defined (0 'blue', 1 'orange')
unset colorbox

stats 'transport_vs_mu.dat' u @x:@y nooutput
N = STATS_blocks-1
Ncurves = 20
if (Ncurves>N) step=1; else step = N/(Ncurves-1);
plot 0 w l lc 'black' lw 1 t '',\
for[in=0:(N-1)/step] 'transport_vs_mu.dat'  i step*in u @x:@y w lp ls 9 lc palette frac (in/(1.0*N/step)) t ''

pause -1