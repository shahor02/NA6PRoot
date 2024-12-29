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

/// \author ruben.shahoyan@cern.ch
/// \brief params for ConfigurableParam

#ifndef NA6P_CONFIGURABLE_KEYVAL_PARAM_H_
#define NA6P_CONFIGURABLE_KEYVAL_PARAM_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

namespace na6p
{
namespace conf
{
struct KeyValParam : public na6p::conf::ConfigurableParamHelper<KeyValParam> {
  std::string input_dir = "none";
  std::string output_dir = "none";

  NA6PParamDef(KeyValParam, "keyval");
};
} // namespace conf
} // namespace na6p

#endif
