#!/bin/bash

g++ GenPYTHIA.cxx -o GenPYTHIA \
    $(pythia8-config --cflags --libs) \
    $(root-config --cflags --libs) \
    $(HepMC3-config --cflags --libs) \
    -lHepMC3 -lHepMC3rootIO
