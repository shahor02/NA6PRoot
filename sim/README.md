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

Generators should be configured in root compilable macros with the generator producer function name being the same as the macro name.
The function signature should be `NA6PGenerator* <name>()`. The macro file name must be passed as `na6psim -g <macro_file_name>+ ...`

### GenBox

### GenParam

Allows to generate single type of particle with its transverse (pt or Mt) and longitudinal (Y or eta) distributions provided via string which are converted to TFormula.
The number of particles to generate per event is either provided explicitly (`setNTracks`) or estimated from provided `setdNdY` (or `setdNdEta`). This number is treated as
as Poisson mean in `setPoisson` was called.
See an example [test/genParamPi.C](../test/genParamPi.C).
The [test/genBgEvent.C](../test/genBgEvent.C) shows an example of CockTail generator production PbPb hadronic (pi, K, p) event for specific SPS energies.

## User hooks

User can provide a macro with the function which is executed at each entry and exit of certain methods, e.g. `NA60PMC::AddParticles` (which allows to modify know particle table and their decays)
or `NA60PMC::selectTracksToSave` (which defines which particles will be saved). In the latter case the function e.g. may get the stack and enforce storing some particles.
The signature of the function must be `int <name>(int entryPoint, bool entering_or_exiting)`. The negative returned value signals a fatal error. E.g. macro `testHooks.C` below can be passed as
`na6psim -u testHooks.C+`:

```
#include <fairlogger/Logger.h>
#include <TVirtualMC.h>
#include "NA6PMCStack.h"
#include "NA6PMC.h"
#include "NA6PSimMisc.h"

int testHooks(int arg, bool inout)
{
  if (arg == UserHook::ADDParticles) {
    LOGP(info, "Calling user hook from {} of the AddParticles method", inout ? "entry":"exit");
  } else if (arg == UserHook::SelectParticles) {
    if (!inout) {
      return 0; // do nothing
    }
    LOGP(info, "Calling user hook from selectTracksToSave method entrance");
    auto mc = (NA6PMC*)TVirtualMC::GetMC();
    auto stack = mc->getMCStack();
    int ntr = stack->GetNtrack(), nPrim = stack->GetNprimary();
    for (int i = nPrim; i < ntr; i++) {
      auto* part = stack->GetParticle(i);
      int idMoth = part->GetFirstMother();
      int mothPdg = -1;
      while (idMoth >= 0) {
	auto* currMoth = stack->GetParticle(idMoth);
	int absPdg = std::abs(currMoth->GetPdgCode());
	if (absPdg == 11 || absPdg == 22 || absPdg == 211 || absPdg == 130 || absPdg == 321 || absPdg == 2212 || absPdg == 2112) {  // stop if a "stable" particle is found in the ancestors
	  break;
	}
	if ((absPdg > 400 && absPdg < 600) || (absPdg > 4000 && absPdg < 6000) || absPdg == 100443 || absPdg == 20443) {
	  part->SetBit(UserHook::KeepParticleBit); // force saving particle i (and its ancestors)
	  break;
	}
	idMoth = currMoth->GetFirstMother();
      }
    }
  } else {
    LOGP(error, "Unknown hook ID {}", arg);
    return -1;
  }
  return 0;
}
```


