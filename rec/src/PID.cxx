// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   PID.cxx
/// @author Ruben Shahoyan
/// @brief  particle ids, masses, names class implementation

#include "PID.h"
#include <cassert>
#include <fairlogger/Logger.h>

//_______________________________
PID::PID(const char* name) : mID(nameToID(name, First))
{
  // construct from the name
  assert(mID < NIDs);
}

//_______________________________
PID PID::PDG2PID(int pdg)
{
  PID pid;
  switch (std::abs(pdg)) {
    case 211:
      pid = PID(Pion);
      break;
    case 321:
      pid = PID(Kaon);
      break;
    case 2212:
      pid = PID(Proton);
      break;
    case 11:
      pid = PID(Electron);
      break;
    case 13:
      pid = PID(Muon);
      break;
    case 1000010020:
      pid = PID(Deuteron);
      break;
    case 1000010030:
      pid = PID(Triton);
      break;
    case 1000020030:
      pid = PID(Helium3);
      break;
    case 1000020040:
      pid = PID(Alpha);
      break;
    case 111:
      pid = PID(PI0);
      break;
    case 22:
      pid = PID(Photon);
      break;
    case 311:
      pid = PID(K0);
      break;
    case 3122:
      pid = PID(Lambda);
      break;
    case 1010010030:
      pid = PID(HyperTriton);
      break;
    case 1010010040:
      pid = PID(Hyperhydrog4);
      break;
    case 3312:
      pid = PID(XiMinus);
      break;
    case 3334:
      pid = PID(OmegaMinus);
      break;
    case 1010020040:
      pid = PID(HyperHelium4);
      break;
    case 1010020050:
      pid = PID(HyperHelium5);
      break;
    default:
      LOGP(error, "pdg code {} not valid, falling back to {}", pdg, pid.getName());
  }
  return pid;
}
