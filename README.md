# TRANSPORT BUDDY (tbd)
This C program computes electronic transport coefficients for the 3D non-interacting one-band models of cuprates high temperature superconductors.

## INSTALL
Installation only requires a C compiler (gcc). In the main directory, run:
```
$ make
```
This should produce the `tbd` executable in the main directory. To copy this executable to your `$HOME/bin/` directory (given you have one) run:
```
$ make install
```
You can increase the speed of the computation by adding flags `-march=native` and `-ffast-math` in the makfile. Be aware that using `-ffast-math` [break strict IEEE compliance](https://stackoverflow.com/questions/7420665/what-does-gccs-ffast-math-actually-do).

## EXAMPLES
You can go in the `example` directory to test your installation:
```
$cd example
$ ../tbd
```
the program will run and output two files:
- `transport_vs_Mu.dat`
- `transport_vs_T.dat`

Both files contain the same data. The first contains blocks for each temperature values, in which the doping changes. The other contains blocks for each doping value, in which the temperature changes. This redundance is intended to make plotting with [gnuplot](http://www.gnuplot.info) easier. Since the loop on Mu is the outer loop, you must wait until the program fnishes to get `transport_vs_Mu.dat`, whereas `transport_vs_T.dat` is continuously updated during the computation.

## Input Parameters
The code reads all adjustable parameters of the model from the 'model.dat' file. These parameters are:

#### band parameters:
- `t`      : first neighbour hopping, defines the energy scale (t=1)
- `tp`     : second neighbour hopping
- `tpp`    : third neighbour hopping
- `tppp`   : fourth neighbour hpping
- `tz`     : interplane hopping

For example, the parameters for Nd-LSCO are t=1, tp=-0.14, tpp=0.07, tppp=0, and tz=0.07.

#### Scattering rates contributions, with [typical ranges]:
- `ETA`    : baseline, should always be non-zero, can be very small : [0.001,0.1]
- `ETAell` : constant mean-free path, proportional to the velocity : [0,0.1]
- `ETAdos` : density of states (DOS), proportional to 1/|v_k| : [0,0.1]
- `ETAaFL` : Fermi-Liquid, proportional to (omega^2+pi^2*T^2) : [0,1]
- `ETAbFL` : Fermi-Liquid and DOS, proportional to (omega^2+pi^2*T^2)/|v_k| : [0,1]
- `ETAaPL` : Planckian Limit, proportional to T : [0,0.3]
- `ETAk`   : Large scattering at the antinode, proportional to cos(2*phi_k)^12 : [0,1]  
- `ETAw`   : Linear in frequency, proportional to omega: [0,1]
- `ETAkw`  : Inverse in frequency, proportional to omega: [0,10] 

WARNING: at the moment there are no safeguards that prevent the scattering to be negative at low energies (because of the linear dependence on omega). This might cause problems at high temperature.

#### loop on chemical potential (doping) and temperature, with [typical ranges]
- `muMin  `: starting chemical potential [-6,6] 
- `muMax  `: ending chemical potential [-6,6] larger than muMin
- `nMu    `: number of chemical potential values. `nMu = 1` will make `muMin` the only value considered.
- `Tmin   `: lowest temperature in units of t, [0.005,0.1], For t=1 referencing 250 meV, T=0.005 is roughly equivalent to 6 K
- `Tmax   `: highest temperature, [0.005,0.1], For t=1 referencing 250 meV, T=0.1 is roughly equivalent to 300 K
- `nT     `: number of temerature considered. `nT = 1` will make `Tmin` the only value considered.
- `logT   `: boolean. 0 will yield equally spaced values of T, 1 will produce equally spaced values on a log scale.

### Riemann Integration parameters (to tradeoff speed and accuracy)
- `nK     `: nb of kx and ky considered.  401
- `nKz    `: nb of kz considered in kz, if tz=0, it will be automatically set to 1
- `nOmega `: nb of omega values evaluated. Ex: 31
- `amplitudeCutoff `: The range of the omega integral is different for each temperature. It ends when the value of fermi function derivative reaches this cutoff, suggested value : 0.005.

## Plotting and managing Units
The Gnuplot scripts provided in the `plotting` and `example` directory contain a few exemples of calculated units to compare results with experiments.
