// NA6PCCopyright

#include "NA6PMuonSpecModular.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoManager.h>
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

// place sensors to station - generalized for n modules per side
void NA6PMuonSpecModular::placeSensors(int modulesPerSide, float pixChipDX, float pixChipDY, float pixChipOffsX, float pixChipOffsY, TGeoVolume* pixelStationVol, TGeoVolume* pixelSensor) {
    
    std::vector<float> alpdx, alpdy;
    
    // Generate positions for n x n grid
    float baseStartX = -(modulesPerSide - 1) * pixChipDX / 2.0f;
    float baseStartY = -(modulesPerSide - 1) * pixChipDY / 2.0f;
    int midPoint = modulesPerSide / 2;

    for (int row = 0; row < modulesPerSide; ++row) {
        for (int col = 0; col < modulesPerSide; ++col) {
            // Determine offset signs based on quadrant
            float offsetX = (row + 1 > midPoint) ? pixChipOffsX / 2.0f : -pixChipOffsX / 2.0f;
            float offsetY = (col + 1 > midPoint) ? -pixChipOffsY / 2.0f : pixChipOffsY / 2.0f;
            
            // Calculate position
            float x = baseStartX + offsetX + col * pixChipDX;
            float y = baseStartY + offsetY + row * pixChipDY;
            
            alpdx.push_back(x);
            alpdy.push_back(y);
            LOGP(info, "Row {} Col {}: sensor at ({}, {})", row, col, x, y);
        }
    }
    
    // Place the sensors
    for (size_t ii = 0; ii < alpdx.size(); ++ii) {
        auto* sensorTransform = new TGeoTranslation(alpdx[ii], alpdy[ii], 0);
        pixelStationVol->AddNode(pixelSensor, composeSensorVolID(ii), sensorTransform);
    }
}


void NA6PMuonSpecModular::createMaterials()
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

void NA6PMuonSpecModular::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  float pixChipDX = 20.f;
  float pixChipDY = 20.f;

  const float EnvelopDXH = 1;
  const float EnvelopDYH = 1;
  const float EnvelopDZH = 1;
  float pixChipDz = param.thicknessMSPlane[0];

  createMaterials();

  auto* sensorShape = new TGeoBBox("SensorShape", pixChipDX / 2, pixChipDY / 2, pixChipDz / 2);
  TGeoVolume* pixelSensor = new TGeoVolume("PixelSensor", sensorShape, NA6PTGeoHelper::instance().getMedium(addName("Silicon")));
  pixelSensor->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Silicon")));

  for (int ist = 0; ist < param.nMSPlanes; ist++) {
    auto stnm = fmt::format("MS{}", ist);

    auto* station = new TGeoBBox((stnm + "SH").c_str(), param.dimXMSPlane[ist] / 2 + EnvelopDXH, param.dimYMSPlane[ist] / 2 + EnvelopDYH, param.thicknessMSPlane[ist] / 2 + EnvelopDZH);
    auto stationSensVol = new TGeoVolume(stnm.c_str(), station, NA6PTGeoHelper::instance().getMedium(addName(param.medMSPlane[ist])));

    LOGP(info, "Creating MS station {} with dimensions: X={} Y={} Z={}", stnm, param.dimXMSPlane[ist], param.dimYMSPlane[ist], param.thicknessMSPlane[ist]);
    int nModulesPerSide = int(param.dimXMSPlane[ist] / (pixChipDX));
    LOGP(info, "N={} DX={} DY={}", nModulesPerSide, pixChipDX, pixChipDY);
    placeSensors(nModulesPerSide, pixChipDX, pixChipDY, param.dimXMSPlaneHole[ist], param.dimYMSPlaneHole[ist], stationSensVol, pixelSensor);

    auto stationSensVolEnv = new TGeoVolume((stnm + "Env").c_str(), station, NA6PTGeoHelper::instance().getMedium(addName("Air")));
    stationSensVolEnv->AddNode(stationSensVol, composeSensorVolID(ist));
    world->AddNode(stationSensVolEnv, composeNonSensorVolID(ist), new TGeoTranslation(param.shiftMS[0] + param.posMSPlaneX[ist], param.shiftMS[1] + param.posMSPlaneY[ist], param.shiftMS[2] + param.posMSPlaneZ[ist]));
  }
}

void NA6PMuonSpecModular::setAlignableEntries()
{
  const auto& param = NA6PLayoutParam::Instance();
  int svolCnt = 0;
  for (int ll = 0; ll < param.nMSPlanes; ++ll) {
    int id = getActiveID() * 100 + svolCnt;
    std::string nm = fmt::format("MS_Lr{}_Sens{}", ll, 0);
    std::string path = fmt::format("/World_1/MS{}Env_{}/MS{}_{}", ll, composeNonSensorVolID(ll), ll, composeSensorVolID(ll));
    gGeoManager->SetAlignableEntry(nm.c_str(), path.c_str(), id);
    LOGP(info, "Adding {} {} as alignable sensor {}", nm, path, id);
    svolCnt++;
  }
}

bool NA6PMuonSpecModular::stepManager(int volID)
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

NA6PMuonSpecModularHit* NA6PMuonSpecModular::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom,
                                      float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

void NA6PMuonSpecModular::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PMuonSpecModular::closeHitsOutput()
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

void NA6PMuonSpecModular::writeHits(const std::vector<int>& remapping)
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
