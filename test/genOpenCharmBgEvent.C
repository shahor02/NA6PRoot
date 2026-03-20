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

  // macro to generate a fixed number of signals (D0, D+, Ds, Lambda_c) in a kinematic range
  //
  // Open charm y and pT parametrization are from fragmentation spectra stored in TH3D histograms
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

  } else if (std::abs(en - 158.) < 3) {
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

NA6PGenerator* addBgEventGenerator(NA6PGenCocktail* genCocktail, int NSignalinAcc, float ptMin, float ptMax, float yMin, float yMax, bool Bg, const char* Part)
{
  float T_piM = 0., y0_piM = 0., ysig_piM = 0., dNdY_piM = 0.;
  float T_piP = 0., y0_piP = 0., ysig_piP = 0., dNdY_piP = 0.;
  float T_KM = 0., y0_KM = 0., ysig_KM = 0., dNdY_KM = 0.;
  float T_KP = 0., y0_KP = 0., ysig_KP = 0., dNdY_KP = 0.;
  float T_Prot = 0., y0_Prot = 0., ysig_Prot = 0., dNdY_Prot = 0.;

  const float ycm = NA6PBeamParam::Instance().getYCM();

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
  // Open charm particles using TH3D parametrization
  const char* na6pRoot = getenv("NA6PROOT_ROOT");
  std::string sourceFile = std::string(na6pRoot) + "/../data/pp0_frag-PtSpectra-Boost.root";

  TFile f(sourceFile.c_str());
  if (!f.IsOpen()) {
    LOGP(fatal, "Failed to open file {}, did you download it? Check out the data README", sourceFile);
  }

  auto addCharmFromTH3 = [&](const char* partName, int pdg, const char* hName, const char* h2Name) {
    TH3D* hPtYEtaHisto = dynamic_cast<TH3D*>(f.Get(hName));
    if (!hPtYEtaHisto) {
      LOGP(fatal, "Histogram {} not found in {}", hName, sourceFile);
    }
    TH2* hPtYProj = dynamic_cast<TH2*>(hPtYEtaHisto->Project3D("yx"));
    if (!hPtYProj) {
      LOGP(fatal, "Failed to project histogram {} to pt-y", hName);
    }
    TH1* hPt = dynamic_cast<TH1*>(hPtYEtaHisto->Project3D("x"));
    TH1* hY = dynamic_cast<TH1*>(hPtYEtaHisto->Project3D("y"));
    TH2* hPtY = dynamic_cast<TH2*>(hPtYProj->Clone(h2Name));
    if (!hPtY) {
      LOGP(fatal, "Failed to clone projected histogram {}", hName);
    }
    hPtY->SetDirectory(nullptr);
    hY->SetDirectory(nullptr);
    hPt->SetDirectory(nullptr);
    
    delete hPtYProj;

    auto gen = new NA6PGenHisto(partName, pdg, NSignalinAcc, false, hPt, hY, 0); //in this example the yCM shift is set to 0 because the y-distribution is already shifted
    genCocktail->addGenerator(gen);
  };

  if (strcmp(Part, "D0") == 0) {
    addCharmFromTH3("D0", 421, "hptyeta421", "h2_421");
  } else if (strcmp(Part, "D+") == 0) {
    addCharmFromTH3("D+", 411, "hptyeta411", "h2_411");
  } else if (strcmp(Part, "Ds") == 0) {
    addCharmFromTH3("Ds", 431, "hptyeta431", "h2_431");
  } else if (strcmp(Part, "LambdaC") == 0) {
    addCharmFromTH3("LambdaC", 4122, "hptyeta4122", "h2_4122");
  } else if (strcmp(Part, "antiD0") == 0) {
    addCharmFromTH3("antiD0", -421, "hptyetam421", "h2_m421");
  } else if (strcmp(Part, "D-") == 0) {
    addCharmFromTH3("D-", -411, "hptyetam411", "h2_m411");
  } else if (strcmp(Part, "antiDs") == 0) {
    addCharmFromTH3("antiDs", -431, "hptyetam431", "h2_m431");
  } else if (strcmp(Part, "antiLambdaC") == 0) {
    addCharmFromTH3("antiLambdaC", -4122, "hptyetam4122", "h2_m4122");
  } else {
    LOGP(fatal, "Unknown particle option '{}'. Supported: D0, D+, Ds, LambdaC, antiD0, D-, antiDs, antiLambdaC", Part);
  }

  return genCocktail;
}

NA6PGenerator* genOpenCharmBgEvent(int nSignalinAcc = 1, const char* Part = "D0", float ptMin = 0., float ptMax = 5., float yMin = 0, float yMax = 6, bool Bg = false)
{
  NA6PGenCocktail* genCockt = new NA6PGenCocktail("openCharmBgCocktail");
  return addBgEventGenerator(genCockt, nSignalinAcc, ptMin, ptMax, yMin, yMax, Bg, Part);
}