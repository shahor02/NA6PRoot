#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TRandom3.h>
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecHit.h"
#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include "NA6PVerTelDigitizer.h"
#endif

void hitsToDigits()
{

  // control histos
  TH2F** hLocX = new TH2F*[5];
  TH2F** hLocY = new TH2F*[5];
  TH2F** hModID = new TH2F*[5];
  for (int jl = 0; jl < 5; jl++) {
    hLocX[jl] = new TH2F(Form("hLocX%d", jl), Form("Local x vs. global coord. Layer %d;x_{glo} (cm);y_{glo} (cm); x_{loc} (cm)", jl), 81, -15., 15., 81, -15., 15.);
    hLocY[jl] = new TH2F(Form("hLocY%d", jl), Form("Local y vs. global coord. Layer %d;x_{glo} (cm);y_{glo} (cm); y_{loc} (cm)", jl), 81, -15., 15., 81, -15., 15.);
    hModID[jl] = new TH2F(Form("hModId%d", jl), Form("Module ID vs. global coord. Layer %d;x_{glo} (cm);y_{glo} (cm)", jl), 8, -15., 15., 8, -15., 15.);
    hLocX[jl]->SetStats(0);
    hLocY[jl]->SetStats(0);
    hModID[jl]->SetStats(0);
  }

  // Process VerTel hits
  TFile* fhVT = TFile::Open("HitsVerTel.root");
  if (!fhVT || fhVT->IsZombie()) {
    LOGP(error, "Cannot open HitsVerTel.root");
    if (fhVT)
      delete fhVT;
  } else {
    TTree* thVT = (TTree*)fhVT->Get("hitsVerTel");
    if (!thVT) {
      LOGP(error, "Cannot find tree 'hitsVerTel' in HitsVerTel.root");
    } else {
      std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
      thVT->SetBranchAddress("VerTel", &vtHitsPtr);
      int nEvVT = thVT->GetEntriesFast();
      NA6PVerTelDigitizer dig;
      dig.init("geometry.root");
      for (int jEv = 0; jEv < nEvVT; jEv++) {
        thVT->GetEvent(jEv);
        for (const auto& hit : vtHits) {
          double xyzLocS[3], xyzLocE[3];
          dig.getHitLocalCoord(hit, xyzLocS, xyzLocE);
          int lay = dig.detID2Layer(hit.getDetectorID());
          double xg = hit.getXIn();
          double yg = hit.getYIn();
          double zg = hit.getZIn();
          hLocX[lay]->SetBinContent(hLocX[lay]->FindBin(xg, yg), xyzLocS[0]);
          hLocY[lay]->SetBinContent(hLocY[lay]->FindBin(xg, yg), xyzLocS[1]);
          hModID[lay]->SetBinContent(hModID[lay]->FindBin(xg, yg), hit.getDetectorID());
        }
        dig.process(vtHits);
      }
      dig.closeDigitsOutput();
    }
  }

  TCanvas* c = new TCanvas("c", "", 1400, 800);
  c->Divide(5, 2);
  for (int jl = 0; jl < 5; jl++) {
    c->cd(jl + 1);
    gPad->SetRightMargin(0.14);
    hLocX[jl]->Draw("colz");
    c->cd(jl + 6);
    gPad->SetRightMargin(0.14);
    hLocY[jl]->Draw("colz");
  }
}
