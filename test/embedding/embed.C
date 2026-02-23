#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "NA6PMCEventHeader.h"
#include <cstdlib>
#include <iostream>

void embed(float& x, float& y, float& z)
{
  const char* evStr = gSystem->Getenv("BKG_EVENT");
  int iEvent = evStr ? atoi(evStr) : 0;

  static TFile* fKBkg = TFile::Open("MCKine_bck.root");
  static TTree* tBkgKin = (TTree*)fKBkg->Get("mckine");
  static NA6PMCEventHeader* bkgmcHead = nullptr;

  if (!bkgmcHead) {
    tBkgKin->SetBranchAddress("header", &bkgmcHead);
  }
  tBkgKin->GetEntry(iEvent);

  x = bkgmcHead->getVX();
  y = bkgmcHead->getVY();
  z = bkgmcHead->getVZ();

  std::cout << "Using background event " << iEvent << std::endl;
  std::cout << "with vertex at x " << x << " y " << y << " z " << z << endl;
}
