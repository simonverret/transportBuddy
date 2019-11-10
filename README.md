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

#### Scattering rates contributions (with typical ranges):
- `ETA`    : baseline, should always be non-zero, can be very small : [0.001,0.1]
- `ETAell` : constant mean-free path, proportional to the velocity : [0,0.1]
- `ETAdos` : density of states (DOS), proportional to 1/|v_k| : [0,0.1]
- `ETAaFL` : Fermi-Liquid, proportional to (omega^2+pi^2*T^2) : [0,1]
- `ETAbFL` : Fermi-Liquid and DOS, proportional to (omega^2+pi^2*T^2)/|v_k| : [0,1]
- `ETAaPL` : Planckian Limit, proportional to T : [0,0.3]
- `ETAk`   : Large scattering at the antinode, proportional to cos(2*phi_k)^12 : [0,1]  
- `ETAw`   : Linear in frequency, proportional to omega: [0,1]
- `ETAkw`  : Inverse in frequency, proportional to omega: [0,10] 

#### loop on chemical potential (doping) and temperature
- `muMin  `: starting chemical potential. Between -6 and 6 
- `muMax  `: ending chemical potential. Between -6 and 6, larger than muMin
- `nMu    `: number of chemical potential values. 1 will make muMin the only value considered.
- `Tmin   `: lowest temperature in units of t. For t=1 referencing 250 meV, T=0.005 is roughly equivalent to 6 K
- `Tmax   `: highest temperature. For t=1 referencing 250 meV, T=0.1 is roughly equivalent to 300 K
- `nT     `: number of temerature considered
- `logT   `: boolean. 0 will yield equally spaced values of T, 1 will produce equally spaced values on a log scale.

### Riemann Integration parameters (tradeoff between speed and accuracy)
- `nK     `: nb of kx and ky considered.  401
- `nKz    `: nb of kz considered in kz, if tz=0, it will be automatically set to 1
- `nOmega `: nb of omega values evaluated. Ex: 31
- `amplitudeCutoff `: The range of the omega integral is different for each temperature. It ends when the value of fermi function derivative reaches this cutoff, suggested value : 0.005.

## Plotting and units
The Gnuplot plotting scripts provided in the `expected_results` directory contain exemples of how to translate the results obtained with transportBuddy in units that compare to experiments. For example, the `Seebeck_vs_p.gp` file allows to plot the Seebeck coefficient as a function of doping (tested with Gnuplot 5.2 patchlevel 6). However, even without gnuplot, the `.gp` files required all numerical factors to go from the theoretical units used in the code to more usual units used in experiments.
