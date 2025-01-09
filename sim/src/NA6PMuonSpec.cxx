// NA6PCCopyright

#include "NA6PMuonSpec.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>

void NA6PMuonSpec::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Silicon");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kCyan + 1);
  }
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
}

void NA6PMuonSpec::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();
  const float EnvelopDXH = 0.5;
  const float EnvelopDYH = 0.5;
  const float EnvelopDZH = 0.5;
  createMaterials();

  for (int ist = 0; ist < param.nMSPlanes; ist++) {
    TGeoShape *station = nullptr, *stationEnv = nullptr;
    auto stnm = fmt::format("MS{}", ist);
    if (param.dimYMSPlane[ist] > 0) {
      stationEnv = new TGeoBBox((stnm + "SHEnv").c_str(), param.dimXMSPlane[ist] / 2 + EnvelopDXH, param.dimYMSPlane[ist] / 2 + EnvelopDYH, param.thicknessMSPlane[ist] / 2 + EnvelopDZH);
      station = new TGeoBBox((stnm + "SH").c_str(), param.dimXMSPlane[ist] / 2, param.dimYMSPlane[ist] / 2, param.thicknessMSPlane[ist] / 2);
    } else {
      stationEnv = new TGeoTube((stnm + "SHEnv").c_str(), 0.0, param.dimXMSPlane[ist] / 2 + EnvelopDXH, param.thicknessMSPlane[ist] / 2 + EnvelopDZH);
      station = new TGeoTube((stnm + "SH").c_str(), 0.0, param.dimXMSPlane[ist] / 2, param.thicknessMSPlane[ist] / 2);
    }
    if (param.dimXMSPlaneHole[ist] > 0) {
      TGeoShape* hole = nullptr;
      if (param.dimYMSPlaneHole[ist] > 0) {
        hole = new TGeoBBox((stnm + "HL").c_str(), param.dimXMSPlaneHole[ist] / 2, param.dimYMSPlaneHole[ist] / 2, param.thicknessMSPlane[ist]);
      } else {
        hole = new TGeoTube((stnm + "HL").c_str(), 0.0, param.dimXMSPlaneHole[ist] / 2, param.thicknessMSPlane[ist]);
      }
      station = new TGeoCompositeShape((stnm + "_HS").c_str(), new TGeoSubtraction(station, hole)); // station with plug hole
    }
    auto stationSensVol = new TGeoVolume(stnm.c_str(), station, NA6PTGeoHelper::instance().getMedium(addName(param.medMSPlane[ist])));
    stationSensVol->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName(param.medMSPlane[ist])));
    auto stationSensVolEnv = new TGeoVolume((stnm + "Env").c_str(), station, NA6PTGeoHelper::instance().getMedium(addName("Air")));
    // put physical station into dummy envelope to avoid precision problems with stepping (very thin objects are not well recognised with large steps)
    stationSensVolEnv->AddNode(stationSensVol, composeSensorVolID(ist));
    world->AddNode(stationSensVolEnv, composeNonSensorVolID(ist), new TGeoTranslation(param.shiftMS[0] + param.posMSPlaneX[ist], param.shiftMS[1] + param.posMSPlaneY[ist], param.shiftMS[2] + param.posMSPlaneZ[ist]));
  }
}

bool NA6PMuonSpec::stepManager(int volID)
{
  int sensID = NA6PModule::volID2SensID(volID);
  if (sensID < 0) {
    LOGP(fatal, "Non-sensor volID={} was provided to stepManager of {}", volID, getName());
  }
  auto mc = TVirtualMC::GetMC();
  bool startHit = false, stopHit = false;
  unsigned char status = 0;
  if (mc->IsTrackEntering()) {
    status |= NA6PBaseHit::kTrackEntering;
  }
  if (mc->IsTrackInside()) {
    status |= NA6PBaseHit::kTrackInside;
  }
  if (mc->IsTrackExiting()) {
    status |= NA6PBaseHit::kTrackExiting;
  }
  if (mc->IsTrackOut()) {
    status |= NA6PBaseHit::kTrackOut;
  }
  if (mc->IsTrackStop()) {
    status |= NA6PBaseHit::kTrackStopped;
  }
  if (mc->IsTrackAlive()) {
    status |= NA6PBaseHit::kTrackAlive;
  }

  // track is entering or created in the volume
  if ((status & NA6PBaseHit::kTrackEntering) || (status & NA6PBaseHit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (NA6PBaseHit::kTrackExiting | NA6PBaseHit::kTrackOut | NA6PBaseHit::kTrackStopped))) {
    stopHit = true;
  }
  // increment energy loss at all steps except entrance
  if (!startHit) {
    mTrackData.mEnergyLoss += mc->Edep();
  }
  if (!(startHit | stopHit)) {
    return false; // do noting
  }
  auto stack = (NA6PMCStack*)mc->GetStack();
  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    mc->TrackMomentum(mTrackData.mMomentumStart);
    mc->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }
  if (stopHit) {
    TLorentzVector positionStop;
    mc->TrackPosition(positionStop);
    // Retrieve the indices with the volume path
    auto* p = addHit(stack->GetCurrentTrackNumber(), sensID, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                     mTrackData.mMomentumStart.Vect(), positionStop.T(),
                     mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    if (mVerbosity > 0) {
      LOGP(info, "{} Tr{} {}", getName(), stack->GetCurrentTrackNumber(), p->asString());
    }
    // register det points in TParticle
    stack->addHit(getActiveIDBit());
    return true;
  }
  return false;
}

NA6PMuonSpecHit* NA6PMuonSpec::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom,
                                      float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

void NA6PMuonSpec::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PMuonSpec::closeHitsOutput()
{
  if (mHitsTree && mHitsFile) {
    mHitsFile->cd();
    mHitsTree->Write();
    delete mHitsTree;
    mHitsTree = 0;
    mHitsFile->Close();
    delete mHitsFile;
    mHitsFile = 0;
  }
}

void NA6PMuonSpec::writeHits(const std::vector<int>& remapping)
{
  int nh = mHits.size();
  for (int i = 0; i < nh; i++) {
    auto& h = mHits[i];
    if (remapping[h.getTrackID()] < 0) {
      LOGP(error, "Track {} hit {} in {} was not remapped!", h.getTrackID(), i, getName());
    }
    h.setTrackID(remapping[h.getTrackID()]);
  }
  if (mHitsTree) {
    mHitsTree->Fill();
  }
}
