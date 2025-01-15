// NA6PCCopyright
#include "NA6PTGeoHelper.h"
#include <TMath.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <fairlogger/Logger.h>

using EProc = NA6PTGeoHelper::EProc;
using ECut = NA6PTGeoHelper::ECut;

const std::unordered_map<EProc, const char*> NA6PTGeoHelper::ProcessIDToName = {
  {EProc::kPAIR, "PAIR"},
  {EProc::kCOMP, "COMP"},
  {EProc::kPHOT, "PHOT"},
  {EProc::kPFIS, "PFIS"},
  {EProc::kDRAY, "DRAY"},
  {EProc::kANNI, "ANNI"},
  {EProc::kBREM, "BREM"},
  {EProc::kHADR, "HADR"},
  {EProc::kMUNU, "MUNU"},
  {EProc::kDCAY, "DCAY"},
  {EProc::kLOSS, "LOSS"},
  {EProc::kMULS, "MULS"},
  {EProc::kCKOV, "CKOV"},
  {EProc::kRAYL, "RAYL"},
  {EProc::kLABS, "LABS"}};

const std::unordered_map<ECut, const char*> NA6PTGeoHelper::CutIDToName = {
  {ECut::kCUTGAM, "CUTGAM"},
  {ECut::kCUTELE, "CUTELE"},
  {ECut::kCUTNEU, "CUTNEU"},
  {ECut::kCUTHAD, "CUTHAD"},
  {ECut::kCUTMUO, "CUTMUO"},
  {ECut::kBCUTE, "BCUTE"},
  {ECut::kBCUTM, "BCUTM"},
  {ECut::kDCUTE, "DCUTE"},
  {ECut::kDCUTM, "DCUTM"},
  {ECut::kPPCUTM, "PPCUTM"},
  {ECut::kTOFMAX, "TOFMAX"}};

void NA6PTGeoHelper::addMedium(const std::string& medName, const std::string& matName, Color_t col)
{
  const auto& matN = matName.empty() ? medName : matName;
  if (mMedPool.find(medName) != mMedPool.end()) {
    throw std::runtime_error(Form("Medium %s was already created", medName.c_str()));
  }
  if (mMatPool.find(matN) == mMatPool.end()) {
    throw std::runtime_error(Form("Material %s does not exist", matN.c_str()));
  }
  mMedPool[medName] = new TGeoMedium(matN.c_str(), mMedPool.size(), mMatPool[matN]);
  mColorPool[medName] = col;
}

TGeoMedium* NA6PTGeoHelper::getMedium(const std::string& medName, bool fatalIfMissing) const
{
  auto m = mMedPool.find(medName);
  if (m == mMedPool.end()) {
    if (fatalIfMissing) {
      LOGP(fatal, "Medium {} was not created", medName);
    } else {
      LOGP(error, "Medium {} was not created", medName);
      return nullptr;
    }
  }
  return m->second;
}

Color_t NA6PTGeoHelper::getMediumColor(const std::string& medName, bool fatalIfMissing) const
{
  auto m = mColorPool.find(medName);
  if (m == mColorPool.end()) {
    if (fatalIfMissing) {
      LOGP(fatal, "Medium {} was not created", medName);
    } else {
      LOGP(error, "Medium {} was not created", medName);
      return kWhite;
    }
  }
  return m->second;
}

TGeoRotation* NA6PTGeoHelper::rotAroundVector(float uX, float uY, float uZ, float ddelta)
{
  ddelta *= TMath::DegToRad();
  double sinDelta = std::sin(ddelta), cosDelta = std::cos(ddelta);
  double oneMinusCosDelta = 1.0 - cosDelta;
  double rmat[9] = {
    oneMinusCosDelta * uX * uX + cosDelta,
    oneMinusCosDelta * uX * uY - sinDelta * uZ,
    oneMinusCosDelta * uX * uZ + sinDelta * uY,
    oneMinusCosDelta * uY * uX + sinDelta * uZ,
    oneMinusCosDelta * uY * uY + cosDelta,
    oneMinusCosDelta * uY * uZ - sinDelta * uX,
    oneMinusCosDelta * uZ * uX - sinDelta * uY,
    oneMinusCosDelta * uZ * uY + sinDelta * uX,
    oneMinusCosDelta * uZ * uZ + cosDelta};
  auto rot = new TGeoRotation();
  rot->SetMatrix(rmat);
  return rot;
}

void NA6PTGeoHelper::loadCutsAndProcessesFromFile(const std::string& fname)
{
  // Based on O2 MaterialManager::loadCutsAndProcessesFromFile
  // Implementation of a method to set cuts and processes as done in AliRoot.
  // The file is expected to contain columns of the form
  // MODNAME LOCALMEDIUMID CUT0 ... CUT9 FLAG0 ... FLAG11
  // where cuts and flags correspond to keys denoted in cutnames and procnames below.

  const int NCUTS = 10;  // number of cut columns expected in file
  const int NFLAGS = 12; // number of process flag columns expected in file

  /// list of cut enumerated in ascending column mode as written in file
  ECut cutnames[NCUTS] = {ECut::kCUTGAM,
                          ECut::kCUTELE,
                          ECut::kCUTNEU,
                          ECut::kCUTHAD,
                          ECut::kCUTMUO,
                          ECut::kBCUTE,
                          ECut::kBCUTM,
                          ECut::kDCUTE,
                          ECut::kDCUTM,
                          ECut::kPPCUTM};

  /// list of process flags enumerated in ascending column mode as written in file
  // missing STRA for the moment
  EProc procnames[NFLAGS - 1] = {EProc::kANNI,
                                 EProc::kBREM,
                                 EProc::kCOMP,
                                 EProc::kDCAY,
                                 EProc::kDRAY,
                                 EProc::kHADR,
                                 EProc::kLOSS,
                                 EProc::kMULS,
                                 EProc::kPAIR,
                                 EProc::kPHOT,
                                 EProc::kRAYL};

  std::ifstream cutfile(fname);
  if (!cutfile.is_open()) {
    LOGP(fatal, "File {} does not exist", fname);
  }
  // reading from file
  float cut[NCUTS]; // to store cut values
  int flag[NFLAGS]; // to store flags
  int itmed, iret;
  char line[500];
  char medName[100];
  std::unordered_map<int, bool> treated;

  while (cutfile.getline(line, 500)) { // Initialise cuts and flags for this line
    for (int i = 0; i < NCUTS; i++) {
      cut[i] = -99;
    }
    for (int i = 0; i < NFLAGS; i++) {
      flag[i] = -99;
    }
    if (strlen(line) == 0) {
      continue;
    }
    // ignore comments marked by *
    if (line[0] == '*') {
      continue;
    }
    // Read the numbers
    iret = sscanf(line, "%s %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d %d",
                  medName, &cut[0], &cut[1], &cut[2], &cut[3], &cut[4], &cut[5], &cut[6], &cut[7], &cut[8],
                  &cut[9], &flag[0], &flag[1], &flag[2], &flag[3], &flag[4], &flag[5], &flag[6], &flag[7],
                  &flag[8], &flag[9], &flag[10], &flag[11]);
    if (!iret) { // nothing read
      continue;
    }

    auto med = getMedium(medName, false);
    if (!med) {
      LOGP(error, "Cannot set cuts for medium {}", medName);
      continue;
    }
    int id = med->GetId();
    if (treated.find(id) != treated.end()) {
      LOGP(fatal, "Cuts for medium {} were already provided", medName);
    }
    treated[id] = true;

    for (int i = 0; i < NCUTS; ++i) {
      if (cut[i] >= 0.) {
        auto it = CutIDToName.find(ECut(i));
        if (it == CutIDToName.end()) {
          continue;
        }
        TVirtualMC::GetMC()->Gstpar(id, it->second, cut[i]);
      }
    }
    // apply process flags
    for (int i = 0; i < NFLAGS - 1; ++i) {
      if (flag[i] >= 0) {
        auto it = ProcessIDToName.find(EProc(i));
        if (it == ProcessIDToName.end()) {
          continue;
        }
        TVirtualMC::GetMC()->Gstpar(id, it->second, flag[i]);
      }
    }
    LOGP(info, "Assigned special cuts to medium {}", medName);
  } // end loop over lines

  auto lst = gGeoManager->GetListOfMedia();
  for (const auto m : *lst) {
    if (treated.find(((TGeoMedium*)m)->GetId()) == treated.end()) {
      LOGP(warn, "Medium {} was not assigned special cuts", m->GetName());
    }
  }
}

const char* NA6PTGeoHelper::getMediumCutName(ECut cut) const
{
  auto it = CutIDToName.find(cut);
  if (it != CutIDToName.end()) {
    return it->second;
  }
  return "UNKNOWN";
}

const char* NA6PTGeoHelper::getPhysicsProcessName(EProc process) const
{
  auto it = ProcessIDToName.find(process);
  if (it != ProcessIDToName.end()) {
    return it->second;
  }
  return "UNKNOWN";
}
