# TRANSPORT BUDDY (tbd)
This C program computes electronic transport coefficients for a non-interacting model of single-band high temperature superconductor cuprates like Nd-LSCO.

## REQUIREMENTS
C compiler (gcc)

## INSTALL
In the main directory, run:
 `$ make`
This should produce an executable `tbd`in the main directory. To copy this executable to your `$HOME/bin/` directory (given you have one) run:
 `$ make install`

## EXAMPLES
You can go the the `example` directory to test your installation:
 `$ cd example`
 `$ ../tbd`
the program will run quickly and output two files:
 `transport_vs_Mu.dat`
 `transport_vs_T.dat`
Both files contain the same data. The first contains blocks for every temperature values, in which the doping changes. The other contains blocks for each doping value, in which the temperature changes. This is only to make plotting with gnuplot easier.

## PLOT EXAMPLES and HANDLING UNITS
Among the most difficult aspects of electronic transport is the handling of units. The Gnuplot plotting scripts provided in the `expected_results` directory contain many exemples of how to translate the results obtained with transportBuddy in units that compare to experiments. For example, the `Seebeck_vs_p.gp` file allows to plot the Seebeck coefficient as a function of doping (tested with Gnuplot 5.2 patchlevel 6). However, even without gnuplot, the `.gp` files required all numerical factors to go from the theoretical units used in the code to more usual units used in experiments.

## MORE INFO: Input Parameters in 'model.dat'
All adjustable parameters of the model are read from the 'model.dat' file.
Anything else in that file will cause errors
These parameters are:

### band parameters:
- `t      `: first neighbour hopping. Defines the energy scale. Ex: 1.0
- `tp     `: second neighbour hopping. Ex: Nd-LSCO: -0.14
- `tpp    `: third neighbour hopping. Ex: Nd-LSCO: 0.07
- `tppp   `: fourth neighbour hpping. Ex: Nd-LSCO: 0.07
- `tz     `: interplane hoppin. Ex: Nd-LSCO: 0.07

### lifetime parameters (scattering rates):
- `ETA    `: baseline, should always be non-zero, can be very small. Ex: 0.00001 to 0.1
- `ETAell `: constant mean-free path - proportional to the velocity. Ex: 0 to 0.1
- `ETAdos `: proportional to the density of states - inversely proportional to the velocity Ex: 0 to 0.1
- `ETAaFL `: Fermi-Liquid a - proportional to omega square. Ex: 0 to 1
- `ETAbFL `: Fermi-Liquid b - proportional to temperature square. Ex: 0 to 1
- `ETAaPL `: Planckian Limit - proportional to temperature Ex: 0 to 0.3

### loop on chemical potential (doping)
- `muMin  `: starting chemical potential. Between -6 and 6 
- `muMax  `: ending chemical potential. Between -6 and 6, larger than muMin
- `nMu    `: number of chemical potential values. 1 will make muMin the only value considered.

### loop on temperatures
- `Tmin   `: lowest temperature in units of t. For t=1 referencing 250 meV, T=0.005 is roughly equivalent to 6 K
- `Tmax   `: highest temperature. For t=1 referencing 250 meV, T=0.1 is roughly equivalent to 300 K
- `nT     `: number of temerature considered
- `logT   `: boolean. 0 will yield equally spaced values of T, 1 will produce equally spaced values on a log scale.

### k integral (large influence on the speed of the simulation)
- `nK     `: nb of k considered in kx and ky401
- `nKz    `: nb of k considered in z, if tz=0, it will be ignored 1

### energy integral (large influence on the speed of the simulation)
- `nOmega `: nb of omega values evaluated. Ex: 31
- `amplitudeCutoff `: spacing between omega is modulated with temperature up the the value where the derivative of the fermi function reaches this cutoff, suggested value : 0.005