#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenParam.h"
#include "NA6PBeamParam.h"
#include "NA6PGenCocktail.h"
#include <TMath.h>
#include <TString.h>
#include <fairlogger/Logger.h>
#include <cmath>
#include <string>
#endif

int findEnergyIndex(float en, const float Elab[], int NEnergy)
{
  for (int i = 0; i < NEnergy; ++i) {
    if (std::abs(en - Elab[i]) < 3) {
      return i;
    }
  }
  LOGP(fatal, "No parameterization defined for E_lab={} GeV/A", en);
  return -1;
}

bool getLightHadronParams(float& T_piM, float& y0_piM, float& ysig_piM, float& dNdY_piM,
                          float& T_piP, float& y0_piP, float& ysig_piP, float& dNdY_piP,
                          float& T_KM, float& y0_KM, float& ysig_KM, float& dNdY_KM,
                          float& T_KP, float& y0_KP, float& ysig_KP, float& dNdY_KP,
                          float& T_Prot, float& y0_Prot, float& ysig_Prot, float& dNdY_Prot)
{
  float en = NA6PBeamParam::Instance().energyPerNucleon;
  // pions and Kaons from  NA49 nucl-ex/0205002
  if (std::abs(en - 40.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 40 GeV/A");
    T_piM = T_piP = 0.169;
    T_KM = 0.232;
    T_KP = 0.226;
    T_Prot = 0.25;

    y0_piM = y0_piP = 0.666;
    y0_KM = 0.694;
    y0_KP = 0.694;
    y0_Prot = 0.907;

    ysig_piM = ysig_piP = 0.872;
    ysig_KM = 0.635;
    ysig_KP = 0.725;
    ysig_Prot = 0.798;

    dNdY_piM = 106.1;
    dNdY_piP = 96.6;
    dNdY_KM = 7.58;
    dNdY_KP = 20.1;
    dNdY_Prot = 39.9;
    return true;
  }

  if (std::abs(en - 158.) < 3) {
    LOGP(info, "Setting hadronic event parameters for PbPb at 158 GeV/A");
    T_piM = T_piP = 0.18;
    T_KM = 0.226;
    T_KP = 0.232;
    T_Prot = 0.31;

    y0_piM = y0_piP = 0.72;
    y0_KM = 0.727;
    y0_KP = 0.839;
    y0_Prot = 0.0;

    ysig_piM = ysig_piP = 1.18;
    ysig_KM = 0.81;
    ysig_KP = 0.88;
    ysig_Prot = 0.0;

    dNdY_piM = 175.4;
    dNdY_piP = 170.1;
    dNdY_KM = 16.8;
    dNdY_KP = 29.6;
    dNdY_Prot = 27.;
    return true;
  }

  LOGP(fatal, "No pi/K/p parameterization defined for E_lab={} GeV/A", en);
  return false;
}

NA6PGenerator* addStrangenessBgEventGenerator(NA6PGenCocktail* genCocktail, float ptMin, float ptMax, float yMin, float yMax, bool Bg)
{
  float T_piM = 0., y0_piM = 0., ysig_piM = 0., dNdY_piM = 0.;
  float T_piP = 0., y0_piP = 0., ysig_piP = 0., dNdY_piP = 0.;
  float T_KM = 0., y0_KM = 0., ysig_KM = 0., dNdY_KM = 0.;
  float T_KP = 0., y0_KP = 0., ysig_KP = 0., dNdY_KP = 0.;
  float T_Prot = 0., y0_Prot = 0., ysig_Prot = 0., dNdY_Prot = 0.;

  const int NEnergy = 5;
  const int NParticles = 5;
  const float Elab[NEnergy] = {20, 30, 40, 80, 158};
  const char* particleName[NParticles] = {"phi", "K0s", "Lambda0", "Omega-", "Xi-"};
  const int pdgStrangePart[NParticles] = {333, 310, 3122, 3334, 3312};
  const float massStrangePart[NParticles] = {1.01946, 0.497611, 1.115683, 1.67245, 1.32171};
  // Parameters for phi, K0S, Lambda, Xi, Omega from NA49 papers:
  // [7] NA49 collaboration, Energy and centrality dependence of anti-p and p production and the
  // anti-Lambda/anti-p ratio in Pb+Pb collisions between 20 AGeV and 158 AGeV, Phys. Rev.
  // C73 (2006) 044910.
  // [8] C. Alt, T. Anticic, B. Baatar, D. Barna, J. Bartke, L. Betev et al., Energy dependence of 𝜙
  // meson production in central pb+ pb collisions at s nn= 6 to 17 gev, Physical Review C 78
  // (2008) 044907.
  // [9] C. Alt, T. Anticic, B. Baatar, D. Barna, J. Bartke, L. Betev et al., Energy dependence of 𝜆 and
  // 𝜉 production in central pb+ pb collisions at 20 a, 30 a, 40 a, 80 a, and 158 a gev measured
  // at the cern super proton synchrotron, Physical Review C 78 (2008) 034918.
  // [10] M. Mitrovski, Omega and anti-omega production in pb+ pb and p+ p collisions at 30-a-gev,
  // 40-a-gev and 158-a-gev, J. Phys. G 30 (2003) S357.
  // Omega- and Omega+ at 40-a-gev were measured together, so we generate only Omega-
  const float Tslope[2][NParticles][NEnergy] = {
    {{196.8, 237.4, 244.6, 239.8, 298.7}, {0, 0, 229, 223.1, 229}, {244, 249, 258, 265, 301}, {0, 0, 218, 0, 267}, {221, 233, 222, 227, 277}},
    {{0, 0, 0, 0, 0}, {0, 0, 226, 217, 226}, {339, 284, 301, 292, 303}, {0, 0, 218, 0, 259}, {311, 277, 255, 321, 0}}};
  const float sigmaRapidity[2][NParticles][NEnergy] = {
    {{0.425, 0.538, 0.696, 0.658, 1.451}, {0, 0, 0.725, 0.792, 0.88}, {0.51, 0.66, 0.91, 0.87, 0}, {0, 0, 0.6, 0, 1.2}, {0.45, 0.56, 0.76, 0.71, 1.18}},
    {{0, 0, 0, 0, 0}, {0, 0, 0.635, 0.705, 0.81}, {0.62, 0.69, 0.77, 0.83, 1.00}, {0, 0, 0.6, 0, 1.0}, {0, 0.76, 0.65, 0.87, 0.73}}};
  const float y0Rapidity[2][NParticles][NEnergy] = {
    {{0.425, 0.538, 0.487, 0.682, 0.0}, {0, 0, 0.694, 0.742, 0.839}, {0.49, 0.59, 0.65, 0.94, 0}, {0, 0, 0, 0, 0}, {0.45, 0.47, 0.54, 0.68, 0}},
    {{0, 0, 0, 0, 0}, {0, 0, 0.569, 0.668, 0.727}, {0, 0, 0, 0, 0}, {0, 0, 0.0, 0, 0}, {0, 0, 0, 0, 0}}};
  const float multiplicity[2][NParticles][NEnergy] = {
    {{1.89, 1.84, 2.55, 4.04, 8.46}, {0, 0, 59.1, 76.9, 103.0}, {27.1, 36.9, 43.1, 50.1, 44.9}, {0, 0, 0.14, 0, 0.43}, {1.50, 2.42, 2.96, 3.80, 4.04}},
    {{0, 0, 0, 0, 0}, {0, 0, 19.2, 32.4, 51.9}, {0.16, 0.39, 0.68, 1.82, 3.07}, {0, 0, 0, 0, 0.19}, {0, 0.12, 0.13, 0.58, 0.66}}};

  const float ycm = NA6PBeamParam::Instance().getYCM();
  const float en = NA6PBeamParam::Instance().energyPerNucleon;
  int iE = findEnergyIndex(en, Elab, NEnergy);
  LOGP(info, "Setting strangeness cocktail at E_lab={} GeV/A", en);

  std::string dndptFun = "x*exp(-sqrt(x*x+[0]*[0])/[1])";
  std::string dndyFun = "exp(-0.5*pow((x-[0]-[1])/[2],2))+exp(-0.5*pow((x-[0]+[1])/[2],2))";

  if (getLightHadronParams(T_piM, y0_piM, ysig_piM, dNdY_piM,
                           T_piP, y0_piP, ysig_piP, dNdY_piP,
                           T_KM, y0_KM, ysig_KM, dNdY_KM,
                           T_KP, y0_KP, ysig_KP, dNdY_KP,
                           T_Prot, y0_Prot, ysig_Prot, dNdY_Prot) &&
      Bg) {
    auto genpiM = new NA6PGenParam("pi-", -211, 0, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, true);
    genpiM->setParametersTrans({0.1396f, T_piM});
    genpiM->setParametersLong({ycm, y0_piM, ysig_piM});
    genpiM->setdNdY(dNdY_piM);
    genCocktail->addGenerator(genpiM);

    auto genpiP = new NA6PGenParam("pi+", 211, 0, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, true);
    genpiP->setParametersTrans({0.1396f, T_piP});
    genpiP->setParametersLong({ycm, y0_piP, ysig_piP});
    genpiP->setdNdY(dNdY_piP);
    genCocktail->addGenerator(genpiP);

    auto genKM = new NA6PGenParam("K-", -321, 0, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, true);
    genKM->setParametersTrans({0.4937f, T_KM});
    genKM->setParametersLong({ycm, y0_KM, ysig_KM});
    genKM->setdNdY(dNdY_KM);
    genCocktail->addGenerator(genKM);

    auto genKP = new NA6PGenParam("K+", 321, 0, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, true);
    genKP->setParametersTrans({0.4937f, T_KP});
    genKP->setParametersLong({ycm, y0_KP, ysig_KP});
    genKP->setdNdY(dNdY_KP);
    genCocktail->addGenerator(genKP);

    auto genProt = new NA6PGenParam("Proton", 2212, 0, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, true);
    genProt->setParametersTrans({0.9383f, T_Prot});
    genProt->setParametersLong({ycm, y0_Prot, ysig_Prot});
    genProt->setdNdY(dNdY_Prot);
    genCocktail->addGenerator(genProt);
  }

  for (int iMatter = 0; iMatter < 2; ++iMatter) {
    for (int iPart = 0; iPart < NParticles; ++iPart) {
      if (iMatter == 1 && (pdgStrangePart[iPart] == 333 || pdgStrangePart[iPart] == 310)) {
        continue;
      }

      const bool isK0s = (pdgStrangePart[iPart] == 310);
      const float tMeV = isK0s ? 0.5 * (Tslope[0][iPart][iE] + Tslope[1][iPart][iE]) : Tslope[iMatter][iPart][iE];
      const float y0 = isK0s ? 0.5 * (y0Rapidity[0][iPart][iE] + y0Rapidity[1][iPart][iE]) : y0Rapidity[iMatter][iPart][iE];
      const float sigmaY = isK0s ? 0.5 * (sigmaRapidity[0][iPart][iE] + sigmaRapidity[1][iPart][iE]) : sigmaRapidity[iMatter][iPart][iE];
      const float mult = isK0s ? 0.5 * (multiplicity[0][iPart][iE] + multiplicity[1][iPart][iE]) : multiplicity[iMatter][iPart][iE];

      if (mult <= 0.) {
        LOGP(warning, "Skipping {} (pdg={}): no parameterization available for E_lab={} GeV/A", particleName[iPart], pdgStrangePart[iPart], en);
        continue;
      }

      const int pdg = (iMatter == 0) ? pdgStrangePart[iPart] : -pdgStrangePart[iPart];
      const float tGeV = 1.e-3 * tMeV;
      TString name = particleName[iPart];
      if (iMatter == 1 && pdgStrangePart[iPart] != 333 && pdgStrangePart[iPart] != 310) {
        name += "bar";
      }

      auto gen = new NA6PGenParam(name.Data(), pdg, 0, dndptFun, dndyFun, ptMin, ptMax, yMin, yMax, true, true, true);
      gen->setParametersTrans({massStrangePart[iPart], tGeV});
      gen->setParametersLong({ycm, y0, sigmaY});
      gen->setMultiplicity(mult);
      genCocktail->addGenerator(gen);

      LOGP(info, "Added {} (pdg={}) T={} GeV sigmaY={} y0={} multiplicity={}", name.Data(), pdg, tGeV, sigmaY, y0, mult);
    }
  }

  return genCocktail;
}

NA6PGenerator* genStrangenessBgEvent(float ptMin = 0., float ptMax = 10., float yMin = 0, float yMax = 6, bool Bg = true)
{
  NA6PGenCocktail* genCockt = new NA6PGenCocktail("strangessBgCocktail");
  return addStrangenessBgEventGenerator(genCockt, ptMin, ptMax, yMin, yMax, Bg);
}