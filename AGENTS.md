# Shell guidance

Safe commands:
- rg
- cmake
- ninja
- ctest
- root -l -b -q

Require approval:
- apt
- pip install
- sudo
- git reset --hard
- git clean -fdx



# Tracking conventions

Units:
- cm
- GeV/c
- kGauss

Coordinate system:
- Z: beam direction
- X: bending plane
- Y: non-bending plane

State vector:
[x, y, tx, ty, q/pXZ]

Definitions:
tx = px / pXZ
ty = py / pXZ
pXZ = sqrt(px^2 + pz^2)

Do NOT replace with:
- q/p
- q/pt
- px/pz conventions
unless explicitly requested.

Track propagation assumes:
- pz > 0
- forward geometry
- detector planes orthogonal to Z


# Covariance layout

Covariance stored as packed lower triangle:
mC[15]

Index helpers:
see NA6PTrackParCov.h / CovarMap and getCovMatElem, setCovMatElem, incCovMatElem helpers

Ordering:
C00
C10 C11
C20 C21 C22
...

Magnetic field:
Jacobian / linear transport matrix assumes only By != 0
Full propagation uses nonzero Bx, By and Bz

Rules:
- preserve packed lower-triangle covariance layout
- no heap allocations in propagation
- prefer float state, double Jacobian evaluation
- avoid Fresnel-integral formulations
- preserve CPU efficiency

Preferences:
- avoid divisions by tiny curvature terms
- when possible, avoid divisions by difference of two large numbers of similar scale

# Preferred coding patterns

- use explicit loops in hot code
- precompute repeated factors
- prefer free functions over virtual dispatch
- preserve ROOT dictionary compatibility
- use constexpr for constants

# Build
Setup prebuilt O2 environment via alias o2d:
o2d

Configure:
cd NA6PRoot;
NA6PRootSrc=`pwd`;
BUILDDIR="$NA6PRootSrc/build"                                     # assuming that we want the build directory to be here ...
INSTALLDIR="$NA6PRootSrc/install"                                 # assuming that we want the installation directory to be here
mkdir -p $BUILDDIR;
cd $BUILDDIR;
cmake -DCMAKE_INSTALL_PREFIX=$INSTALLDIR $NA6PRootSrc

Build
cd $BUILDDIR;
make -j5 install

Setup include and library paths
source $INSTALLDIR/init.sh
