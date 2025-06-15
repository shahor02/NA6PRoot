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

/// @brief Class for creating debug root trees with std::cout like intervace
/// @author Marian Ivanov, marian.ivanov@cern.ch (original code in AliRoot)
///         Ruben Shahoyan, ruben.shahoyan@cern.ch (porting to O2)

#ifndef NA6P_TREESTREAMREDIRECTOR_H
#define NA6P_TREESTREAMREDIRECTOR_H

#include <Rtypes.h>
#include <TDirectory.h>
#include "NA6PTreeStream.h"

/// The NA6PTreeStreamRedirector class manages one or few NA6PTreeStream objects to be written to
/// the same output file.
/// NA6PTreeStreamRedirector myNA6PTreeStreamRedirector("myOutFile.root","recreate");
/// myNA6PTreeStreamRedirector<<"myStream0"<<"brName00="<<obj00<<"brName01="<<obj01<<"\n";
/// ...
/// myNA6PTreeStreamRedirector<<"myStream2"<<"brName10="<<obj10<<"brName11="<<obj11<<"\n";
/// ...
/// will create ouput file with 2 trees stored.
///
/// The flushing of trees to the file happens on NA6PTreeStreamRedirector::Close() call
/// or at its desctruction.
///

class NA6PTreeStreamRedirector
{
 public:
  NA6PTreeStreamRedirector(const char* fname = "", const char* option = "recreate");
  virtual ~NA6PTreeStreamRedirector();
  void Close();
  TFile* GetFile() { return mDirectory->GetFile(); }
  TDirectory* GetDirectory() { return mDirectory; }
  virtual NA6PTreeStream& operator<<(Int_t id);
  virtual NA6PTreeStream& operator<<(const char* name);
  void SetDirectory(TDirectory* sfile);
  void SetFile(TFile* sfile);
  static void FixLeafNameBug(TTree* tree);

 private:
  NA6PTreeStreamRedirector(const NA6PTreeStreamRedirector& tsr);
  NA6PTreeStreamRedirector& operator=(const NA6PTreeStreamRedirector& tsr);

  std::unique_ptr<TDirectory> mOwnDirectory;             // own directory of the redirector
  TDirectory* mDirectory = nullptr;                      // output directory
  std::vector<std::unique_ptr<NA6PTreeStream>> mDataLayouts; // array of data layouts

  ClassDefNV(NA6PTreeStreamRedirector, 0);
};

#endif
