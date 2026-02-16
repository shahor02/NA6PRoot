<!-- doxy
\page refREC rec
/doxy -->

# Reconstruction
---

Use the executable:

```
na6prec <options> --configKeyValues "<semicolon-separated configurable params>"
```

### Layout definition

The detector layout is defined by the `ConfigurableParam`-based [NA6PRecoParam](../rec/include/NA6PRecoParam.h) class. Any parameter defined in this class can be modified from the command line:

```
na6prec <options> --configKeyValues "rec.msNIterationsTrackerCA=3;rec.zMatching=170;<other settings>"
```

After execution,  a `na6pRecoParam.ini` file will be created containing the reconstruction settings.

Instead of passing a long `--configKeyValues` string on the command line, you can generate the `na6pRecoParam.ini` file once (even with default options) and then edit the values in the `[reco]` block.

To run the reconstruction with modified parameters (and optionally override some of them from the command line), use:

```
na6prec <options> --load-ini <user-file.ini> --configKeyValues "<semicolon-separated configurable params>"
```

The priority order of parameters is:

1. values provided via `--configKeyValues`
2. values loaded from the `.ini` file
3. class default values

Writing the `.ini` file can be disabled with:

```
--disable-write-ini
```
### Execution control options

In addition to the parameters in `NA6PRecoParam`, the executable provides options to control the event range.

#### Event range

```
--firstevent, -f <int>
--lastevent,  -l <int>
```

* `--firstevent` (`-f`): first event to process. Default: `0`.
* `--lastevent` (`-l`): last event to process (inclusive). Default: `-1` (process until the end of the file).

Example:

```
na6prec <options> -f 100 -l 999
```

Processes events from 100 to 999.

#### Reconstruction step switches

Each main reconstruction stage can be enabled or disabled from the command line:

```
--doHitsToRecPoints, -cl <bool>
--doTrackletVertex,  -vert <bool>
--doVTTracking,      -vt <bool>
--doMSTracking,      -ms <bool>
--doMatching,        -mt <bool>
```

* `--doHitsToRecPoints` (`-cl`): run the hits → clusters (rec-points) step.
  Default: `true`.

* `--doTrackletVertex` (`-vert`): run the tracklet-based vertex reconstruction.
  Default: `true`.

* `--doVTTracking` (`-vt`): run the Vertex Telescope (VT) tracking.
  Default: `true`.

* `--doMSTracking` (`-ms`): run the Muon Spectrometer (MS) tracking.
  Default: `true`.

* `--doMatching` (`-mt`): run the matching between VT and MS tracks.
  Default: `true`.

Example:

```
na6prec <options> --doMSTracking false --doMatching false
```

This runs the cluster creation, trackelt vertexing and VT reconstruction without MS tracking and VT–MS matching.

---
## Clusters

NA6PBaseCluster class defines the basic properties of cluster objects.
It contains the data members and the functions that are needed for the vertexing and tracking classes.
Derived classes with specific additional information for VT and MS clusters can be implemented if needed.

