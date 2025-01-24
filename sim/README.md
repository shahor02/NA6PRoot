<!-- doxy
\page refSIM sim
/doxy -->

# Simulation

Use executable 
```
na6psim <options> --configKeyValues "<semicolon separated configurable params>"
```

### Layout definition
The detector layout is defined by the ConfigurableParam-based [NA6PLayoutParam](../base/include/NA6PLayoutParam.h) class. Any parameter defined in this class can be modified from the command line. For instance, 
to set the `Z` of the 1-st VerTel station to `8.0` cm instead of its default value and then to move the whole VerTel in X and Z globally by `0.1` and `-0.2` cm respectively and change the 2nd target material to Iron (make sure that
corresponding material and medium are defined in the `createMaterials()` method) one can pass an option:
```
na6psim <options> --configKeyValues "layout.shiftVerTel[0]=0.1;layout.shiftVerTel[2]=-0.2;layout.posVerTelPlaneZ[0]=8.0;layout.medTarget[1]=Iron;<other settings>"
```
After the execution a corresponding `geometry.root` is generated in the working directory as well as `na6pLayout.ini` text file with the settings used for layout (the directory where this file is generated can be modified by `keyval.output_dir=<directory>` configurable).
Instead of passing long configKeyValues options string on the command line one can generate this `na6pLayout.ini` once (even with default options) and then edit the values of the `[layout]` block.
To run simulations with this modified layout (as well as with other configurable parames) one can use
```
na6psim <options> --load-ini <user-file.ini> --configKeyValues "<semicolon separated configurable params>"
```
Note that even in this case, one can still modify some parameters via `--configKeyValues ...` settings: the priority is given to the value provided from the command line, then loaded from the ini file (if any), and then to the class default setting.

Writing of the ini file can be disable by the option `--disable-write-ini`.


## Available generators

### GenBox

### GenParam

Allows to generate single type of particle with its transverse (pt or Mt) and longitudinal (Y or eta) distributions provided via string which are converted to TFormula.
The number of particles to generate per event is either provided explicitly (`setNTracks`) or estimated from provided `setdNdY` (or `setdNdEta`). This number is treated as
as Poisson mean in `setPoisson` was called.
See an example [test/genParamPi.C](../test/genParamPi.C).
