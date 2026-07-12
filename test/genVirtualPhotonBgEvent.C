#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenParam.h"
#include "NA6PGenHisto.h"
#include "NA6PBeamParam.h"
#include "NA6PGenCocktail.h"
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3D.h>
#include <cstdlib>
#include <fairlogger/Logger.h>
#endif

bool getParams(float& T_piM, float& y0_piM, float& ysig_piM, float& dNdY_piM,
               float& T_piP, float& y0_piP, float& ysig_piP, float& dNdY_piP,
               float& T_KM, float& y0_KM, float& ysig_KM, float& dNdY_KM,
               float& T_KP, float& y0_KP, float& ysig_KP, float& dNdY_KP,
               float& T_Prot, float& y0_Prot, float& ysig_Prot, float& dNdY_Prot)
{
  float en = NA6PBeamParam::Instance().energyPerNucleon;

  // Background parametrization is available only for 40 and 150 Gev and it can be switched on/off through the argument Bg

  if (std::abs(en - 40.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 40. GeV/c");
    // pions and Kaons from  NA49 nucl-ex/0205002
    T_piM = T_piP = 0.169;
    T_KM = 0.232;
    T_KP = 0.226;
    T_Prot = 0.25;

    y0_piM = y0_piP = 0.666; // pions 2-gaussian poles half distance
    y0_KM = 0.694;           // K- 2-gaussian poles half distance
    y0_KP = 0.694;           // K+ 2-gaussian poles half distance
    y0_Prot = 0.907;         // ???

    ysig_piM = ysig_piP = 0.872;
    ysig_KM = 0.635;
    ysig_KP = 0.725;
    ysig_Prot = 0.798;

    dNdY_piM = 106.1;
    dNdY_piP = 96.6;
    dNdY_KM = 7.58;
    dNdY_KP = 20.1;
    dNdY_Prot = 39.9;
  } else if (std::abs(en - 50.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 50. GeV/c");
    LOGP(info, "No hadronic event parameters for PbPb at 50. GeV/c");
  } else if (std::abs(en - 70.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 70. GeV/c");
    LOGP(info, "No hadronic event parameters for PbPb at 70. GeV/c");
  } else if (std::abs(en - 90.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 90. GeV/c");
    LOGP(info, "No hadronic event parameters for PbPb at 90. GeV/c");
  } else if (std::abs(en - 110.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 110. GeV/c");
    LOGP(info, "No hadronic event parameters for PbPb at 110. GeV/c");
  } else if (std::abs(en - 130.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 130. GeV/c");
    LOGP(info, "No hadronic event parameters for PbPb at 130. GeV/c");
  } else if (std::abs(en - 158.) < 3 || std::abs(en - 150.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 158. GeV/c");

    T_piM = T_piP = 0.18;
    T_KM = 0.226;
    T_KP = 0.232;
    T_Prot = 0.31;

    y0_piM = y0_piP = 0.72; // pions 2-gaussian poles half distance
    y0_KM = 0.727;          // K- 2-gaussian poles half distance
    y0_KP = 0.839;          // K+ 2-gaussian poles half distance
    y0_Prot = 39.8;         // ???

    ysig_piM = ysig_piP = 1.18;
    ysig_KM = 0.81;
    ysig_KP = 0.88;
    ysig_Prot = 8.07;

    dNdY_piM = 175.4;
    dNdY_piP = 170.1;
    dNdY_KM = 16.8;
    dNdY_KP = 29.6;
    dNdY_Prot = 27.;
  } else {
    LOGP(fatal, "No parameters defined for beam energy {}", en);
  }
  return true;
}

NA6PGenerator* addBgEventGenerator(NA6PGenCocktail* genCocktail, TH3* h3, int pdg, int NSignalinAcc, float ptMin, float ptMax, float yMin, float yMax, bool Bg)
{
  float T_piM = 0., y0_piM = 0., ysig_piM = 0., dNdY_piM = 0.;
  float T_piP = 0., y0_piP = 0., ysig_piP = 0., dNdY_piP = 0.;
  float T_KM = 0., y0_KM = 0., ysig_KM = 0., dNdY_KM = 0.;
  float T_KP = 0., y0_KP = 0., ysig_KP = 0., dNdY_KP = 0.;
  float T_Prot = 0., y0_Prot = 0., ysig_Prot = 0., dNdY_Prot = 0.;

  const float en = NA6PBeamParam::Instance().energyPerNucleon;
  const float ycm = NA6PBeamParam::Instance().getYCM();
  LOGF(info, "en = %3.1f GeV, ycm = %3.2f", en, ycm);

  std::string dndptBgFun = "x*exp(-sqrt(x*x+[0]*[0])/[1])";
  std::string dndyBgFun = "exp(-0.5*pow((x-[0]-[1])/[2],2))+exp(-0.5*pow((x-[0]+[1])/[2],2))";

  getParams(T_piM, y0_piM, ysig_piM, dNdY_piM,
            T_piP, y0_piP, ysig_piP, dNdY_piP,
            T_KM, y0_KM, ysig_KM, dNdY_KM,
            T_KP, y0_KP, ysig_KP, dNdY_KP,
            T_Prot, y0_Prot, ysig_Prot, dNdY_Prot);

  bool IsPoisson = true; 	

  if (Bg) {
    auto genpiM = new NA6PGenParam("pi-", -211, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genpiM->setParametersTrans({0.1396f, T_piM});
    genpiM->setParametersLong({ycm, y0_piM, ysig_piM});
    genpiM->setdNdY(dNdY_piM);
    genCocktail->addGenerator(genpiM);

    auto genpiP = new NA6PGenParam("pi+", 211, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genpiP->setParametersTrans({0.1396f, T_piP});
    genpiP->setParametersLong({ycm, y0_piP, ysig_piP});
    genpiP->setdNdY(dNdY_piP);
    genCocktail->addGenerator(genpiP);

    auto genKM = new NA6PGenParam("K-", -321, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genKM->setParametersTrans({0.4937f, T_KM});
    genKM->setParametersLong({ycm, y0_KM, ysig_KM});
    genKM->setdNdY(dNdY_KM);
    genCocktail->addGenerator(genKM);

    auto genKP = new NA6PGenParam("K+", 321, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genKP->setParametersTrans({0.4937f, T_KP});
    genKP->setParametersLong({ycm, y0_KP, ysig_KP});
    genKP->setdNdY(dNdY_KP);
    genCocktail->addGenerator(genKP);

    auto genProt = new NA6PGenParam("Proton", 2212, 0, dndptBgFun, dndyBgFun, ptMin, ptMax, yMin, yMax, true, true, IsPoisson);
    genProt->setParametersTrans({0.9383f, T_Prot});
    genProt->setParametersLong({ycm, y0_Prot, ysig_Prot});
    genProt->setdNdY(dNdY_Prot);
    genCocktail->addGenerator(genProt);
  }

  auto gen = new NA6PGenHisto("VirtualPhoton", 23, NSignalinAcc, false, h3, ycm);
  genCocktail->addGenerator(gen);
  return genCocktail;
}

NA6PGenerator *genVirtualPhotonBgEvent(const char* fileName = "dy.root", const char* histName = "h3MPtY", int nSignalinAcc = 1, float ptMin = 0, float ptMax = 5, float yMin = 0, float yMax = 6, bool Bg = false)  
{
  // macro to generate a fixed number of signals for virtual photon (= prompt dimon) based on theoretical model in d3sigma/dmdpTdy or d3N/dmdpTdy.

  TFile* rootfile = TFile::Open(fileName, "READ");
  TH3F *h3 = (TH3F*)rootfile->Get(histName); // in CM frame.
  h3->SetDirectory(0);
  rootfile->Close();

  std::cout << "Before revmoving negative bins: Integral = " << h3->Integral() << '\n';

  for (int ix = 1; ix <= h3->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h3->GetNbinsY(); ++iy) {
      for (int iz = 1; iz <= h3->GetNbinsZ(); ++iz) {
        double content = h3->GetBinContent(ix, iy, iz);
        if (content < 0.0) {
          double xCenter = h3->GetXaxis()->GetBinCenter(ix);
          double yCenter = h3->GetYaxis()->GetBinCenter(iy);
          double zCenter = h3->GetZaxis()->GetBinCenter(iz);
          std::cout << "Negative bin:" << " ix=" << ix << " iy=" << iy << " iz=" << iz << " | x=" << xCenter << " y=" << yCenter << " z=" << zCenter << " | content=" << content << '\n';
          h3->SetBinContent(ix, iy, iz, 0);
        }
      }
    }
  }

  std::cout << "After revmoving negative bins: Integral = " << h3->Integral() << '\n';

  if (h3->Integral() <= 0.0) {
    throw std::runtime_error("TH3 integral is not positive");
  }

  NA6PGenCocktail* genCockt = new NA6PGenCocktail("testCocktail");
  return addBgEventGenerator(genCockt, h3, 23, nSignalinAcc, ptMin, ptMax, yMin, yMax, Bg);
}
