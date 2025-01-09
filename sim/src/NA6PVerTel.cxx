// NA6PCCopyright

#include "NA6PVerTel.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoParaboloid.h>
#include <TGeoTrd1.h>
#include <TGeoTrd2.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>

void NA6PVerTel::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Silicon");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kCyan + 1);
  }
  nameM = addName("CarbonFoam");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 0.5);
    mixt->AddElement(12.01, 6, 1.0); // Carbon-only mixture
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 6);
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

void NA6PVerTel::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  float frameDX = 30.0f;
  float frameDY = 30.0f;
  float frameDZ = 0.50f;
  float frameHoleR = 0.424f;

  float pixChipHoleDX = 12.5f;
  float pixChipHoleDY = 13.0f;
  float dxyCut = 1.0f;

  float pixChipContainerDX = frameDX;
  float pixChipContainerDY = frameDY;
  float pixChipContainerDz = 0.05f;

  float pixChipDX = 14.69f;
  float pixChipDY = 14.69f;
  float pixChipDz = 50e-4f;
  float pixChipOffsX = 0.29f;
  float pixChipOffsY = 0.31f;

  float boxDZMargin = pixChipContainerDz + 0.5f;
  float boxDZ = param.posVerTelPlaneZ[param.nVerTelPlanes - 1] - param.posVerTelPlaneZ[0] + 2 * boxDZMargin;

  std::vector<float> chipHoleX = {pixChipDX / 2 + dxyCut, -pixChipDX / 2, -pixChipDX / 2 - dxyCut, pixChipDX / 2};
  std::vector<float> chipHoleY = {pixChipDY / 2, pixChipDY / 2 + dxyCut, -pixChipDY / 2, -pixChipDY / 2 - dxyCut};

  // Container
  auto* vtShape = new TGeoBBox("VTContainer", (frameDX + 2.f) / 2, (frameDY + 2.f) / 2, boxDZ / 2);
  TGeoVolume* vtContainer = new TGeoVolume("VTContainer", vtShape, NA6PTGeoHelper::instance().getMedium(addName("Air")));
  auto* vtTransform = new TGeoCombiTrans(param.shiftVerTel[0],
                                         param.shiftVerTel[1],
                                         param.shiftVerTel[2] + (param.posVerTelPlaneZ[0] + boxDZ / 2 - boxDZMargin),
                                         NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(vtContainer, composeNonSensorVolID(0), vtTransform);

  // pixel station Frame with holes (box with subtracted holes)
  auto* pixStFrameBox = new TGeoBBox("PixStFrameBox", frameDX / 2, frameDY / 2, frameDZ / 2);
  auto* beamPipeHole = new TGeoTube("PixStFrameBoxBPHole", 0, frameHoleR, frameDZ);
  auto* frameSubtraction = new TGeoSubtraction(pixStFrameBox, beamPipeHole);
  auto* pixStFrameShape = new TGeoCompositeShape("PixStFrameBoxHole0", frameSubtraction);
  for (size_t ii = 0; ii < chipHoleX.size(); ++ii) {
    auto* pixChipHole = new TGeoBBox("PixChipHole", pixChipHoleDX / 2, pixChipHoleDY / 2, frameDZ);
    auto* holeTransform = new TGeoTranslation(chipHoleX[ii], chipHoleY[ii], 0);
    frameSubtraction = new TGeoSubtraction(pixStFrameShape, pixChipHole, nullptr, holeTransform);
    pixStFrameShape = new TGeoCompositeShape(Form("PixStFrameBoxHole0%zu", ii), frameSubtraction);
  }
  TGeoVolume* pixStFrame = new TGeoVolume("PixStFrame", pixStFrameShape, NA6PTGeoHelper::instance().getMedium(addName("CarbonFoam")));
  pixStFrame->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("CarbonFoam")));

  // Silicon Tracker Station
  auto* pixelStationShape = new TGeoBBox("PixelStationShape", pixChipContainerDX / 2, pixChipContainerDY / 2, pixChipContainerDz / 2);
  auto* sensorShape = new TGeoBBox("SensorShape", pixChipDX / 2, pixChipDY / 2, pixChipDz / 2);
  auto* pixelStationVol = new TGeoVolume("PixelStationVol", pixelStationShape, NA6PTGeoHelper::instance().getMedium(addName("Air")));
  TGeoVolume* pixelSensor = new TGeoVolume("PixelSensor", sensorShape, NA6PTGeoHelper::instance().getMedium(addName("Silicon")));
  pixelSensor->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Silicon")));

  // place sensors to station
  std::vector<float> alpdx{pixChipDX / 2 + pixChipOffsX, -pixChipDX / 2 + pixChipOffsY, -pixChipDX / 2 - pixChipOffsX, pixChipDX / 2 - pixChipOffsY};
  std::vector<float> alpdy{pixChipDY / 2 - pixChipOffsY, pixChipDY / 2 + pixChipOffsX, -pixChipDY / 2 + pixChipOffsY, -pixChipDY / 2 - pixChipOffsX};
  for (size_t ii = 0; ii < alpdx.size(); ++ii) {
    auto* sensorTransform = new TGeoTranslation(alpdx[ii], alpdy[ii], 0);
    pixelStationVol->AddNode(pixelSensor, composeSensorVolID(ii), sensorTransform);
  }
  // place frames + sensor stations
  float zoffs = param.posVerTelPlaneZ[0] + boxDZ / 2 - boxDZMargin; // offset to be added due to the placement of stations to the VT box
  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    auto* stationTransform = new TGeoCombiTrans(param.posVerTelPlaneX[ll], param.posVerTelPlaneY[ll], param.posVerTelPlaneZ[ll] - zoffs,
                                                NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
    vtContainer->AddNode(pixelStationVol, composeNonSensorVolID(ll), stationTransform);
    auto* frameTransform = new TGeoCombiTrans(param.posVerTelPlaneX[ll], param.posVerTelPlaneY[ll], param.posVerTelPlaneZ[ll] + 0.5 * (frameDZ + pixChipDz) - zoffs,
                                              NA6PTGeoHelper::rotAroundVector(0, 0.0, 0.0, 0.0));
    vtContainer->AddNode(pixStFrame, composeNonSensorVolID(ll + 20), frameTransform);
  }
}

bool NA6PVerTel::stepManager(int volID)
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
    int stationID(-1);
    mc->CurrentVolOffID(1, stationID);
    stationID = volID2NonSensID(stationID);

    int chipindex = NChipsPerStation * stationID + sensID;
    auto* p = addHit(stack->GetCurrentTrackNumber(), chipindex, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
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

NA6PVerTelHit* NA6PVerTel::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom,
                                  float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

void NA6PVerTel::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PVerTel::closeHitsOutput()
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

void NA6PVerTel::writeHits(const std::vector<int>& remapping)
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
