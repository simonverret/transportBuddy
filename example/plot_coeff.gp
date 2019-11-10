w = 800
h = 600
set term qt font "Arial,16" persist size w,h

# w = 10.5
# h = 8
# set terminal epslatex size w cm,h cm color dashlength 0.3
# set output 'plot.tex'

borders_width = 0.5
# set termopt enhanced 


set title \
'$\Gamma_{\bm k}= 0.12t + 0.2t\omega  + 0.52\cos(\phi_{\bm k})^{12}$'
# "$t=0.190$ eV\n$t'=-0.14t$\n$t''=0.07t$\n$t_z=-0.07t$\n".\

echarge = 1.602e-19
hopping = 0.190           # eV  -- DO NOT CHANGE, because implicit in other factors below
boltzmann = 8.62e-5       # eV/K
Tfac = hopping/boltzmann
sfac = 1.845e5            # (Ohm * m)^-1
afac = 8.617e-5*sfac      # A * (K * m)^-1
nHfac = 1/sfac
CvFac = 2.866             # mJ (K^2 mol)-1 
RHfac = 0.0002137         # m^3 / C

T  = "(column('T')*Tfac)"
mu = "(column('mu'))"
p  = "(column('p0'))"
n  = "(column('n0'))"

sxx     = "(sfac*column('sigmaxx'))"            # (Ohm * m)^-1
szz     = "(sfac*column('sigmazz'))"            # (Ohm * m)^-1
sxy     = "(sfac*column('sigmaxy'))"            # (Ohm * m)^-1   / B [T ?]
axx     = "(afac*column('alphaxx'))"            # A * (Ohm * m)^-1
azz     = "(afac*column('alphazz'))"            # A * (Ohm * m)^-1
axy     = "(column('alphaxy'))"         #TODO
bxx     = "(column('betaxx'))"          #TODO
bxy     = "(column('betaxy'))"          #TODO

dos     = "(column('dos'))"                     # at E_F
Cv      = "(CvFac*column('Cv')/column('T'))"    # mJ (K^2 mol)-1 
rxx     = '(10e7/'.sxx.')'                      # muOhm * cm
rxy     =  '('.sxy.'/('.sxx.'*'.sxx.'))'        # Ohm * m
nH      =  '(nHfac/'.rxy.')'                    # 1+p
RH      =  '(1e9*RHfac*'.rxy.')'                # mm^3 / C
Sx      = '(-10e6*'.axx.'/'.sxx.'/'.T.')'       # muV / K^2
Sz      = '(-10e6*'.azz.'/'.szz.'/'.T.')'       # muV / K^2





x1 = p
set xrange [-1:1]
set xtics .2
set xlabel 'doping'


nsxx     = "(10e-7 * sfac*column('sigmaxx'))"           # (Ohm * m)^-1
nszz     = "(10e-4 * sfac*column('sigmazz'))"           # (Ohm * m)^-1
naxx     = "(10e0 * afac*column('alphaxx'))"            # A * (Ohm * m)^-1
nazz     = "(10e3 * afac*column('alphazz'))"            # A * (Ohm * m)^-1


y1 = nsxx ; title1 = '$10^{-7}\sigma_{xx}$' 
y2 = nszz ; title2 = '$10^{-4}\sigma_{zz}$' 
y3 = naxx ; title3 = '$10^{0}\alpha_{xx}$' 
y4 = nazz ; title4 = '$10^{3}\alpha_{zz}$' 
set yrange [-0.3:0.3] 
set ytics .1
set ylabel '$\sigma$ [($\Omega$ m)$^{-1}$], $\alpha$ [A($\Omega$ m)$^{-1}$]'

set linestyle 1  lc rgb 'blue'       lw 3 dt 1 ps 0.7 pt 7
set linestyle 2  lc rgb 'red'        lw 3 dt 1 ps 0.7 pt 7
set linestyle 3  lc rgb '#009900'    lw 3 dt 1 ps 0.7 pt 7
set linestyle 4  lc rgb 'orange'    lw 3 dt 1 ps 0.7 pt 7
set linestyle 5  lc rgb 'black'      lw 3 dt 1 ps 0.7 pt 7
set linestyle 6  lc rgb '#888888'    lw 3 dt 1 ps 0.7 pt 7

set key bottom right samplen 1 

plot \
0 w l ls 6 t ''\
, 'transport_vs_Mu.dat' u @x1:@y1 w l ls 1 t title1 \
, 'transport_vs_Mu.dat' u @x1:@y2 w l ls 2 t title2 \
, 'transport_vs_Mu.dat' u @x1:@y3 w l ls 3 t title3 \
, 'transport_vs_Mu.dat' u @x1:@y4 w l ls 4 t title4 \



pause -1