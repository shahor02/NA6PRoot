// NA6PCCopyright

#include "NA6PVerTel.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoXtru.h>
#include <TGeoParaboloid.h>
#include <TGeoTrd1.h>
#include <TGeoTrd2.h>
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

static TGeoXtru* CreateRoundedRect(const char* name, double dx, double dy, double r, double dz)
{
  const int nSegments = 20;
  int nPoints = 4 * nSegments;
  double* x = new double[nPoints];
  double* y = new double[nPoints];

  for (int i = 0; i < nSegments; i++) {
    double phi = i * (TMath::Pi() / 2) / (nSegments - 1);
    x[i] = (dx - r) + r * std::cos(phi);
    y[i] = (dy - r) + r * std::sin(phi);
    x[i + nSegments] = -(dx - r) + r * std::cos(phi + TMath::Pi() / 2);
    y[i + nSegments] = (dy - r) + r * std::sin(phi + TMath::Pi() / 2);
    x[i + 2 * nSegments] = -(dx - r) + r * std::cos(phi + TMath::Pi());
    y[i + 2 * nSegments] = -(dy - r) + r * std::sin(phi + TMath::Pi());
    x[i + 3 * nSegments] = (dx - r) + r * std::cos(phi + 3 * TMath::Pi() / 2);
    y[i + 3 * nSegments] = -(dy - r) + r * std::sin(phi + 3 * TMath::Pi() / 2);
  }

  TGeoXtru* xtru = new TGeoXtru(2);
  xtru->SetName(name);
  xtru->DefinePolygon(nPoints, x, y);
  xtru->DefineSection(0, -dz);
  xtru->DefineSection(1, dz);
  delete[] x;
  delete[] y;
  return xtru;
}

static TGeoXtru* CreateStadium(const char* name, double halfLen, double halfH, double dz)
{
  const int nSemi = 30;
  int nPoints = 2 * nSemi;
  double* x = new double[nPoints];
  double* y = new double[nPoints];

  // Right semicircle
  for (int i = 0; i < nSemi; i++) {
    double phi = -TMath::Pi() / 2 + i * TMath::Pi() / (nSemi - 1);
    x[i] = halfLen + halfH * std::cos(phi);
    y[i] = halfH * std::sin(phi);
  }
  // Left semicircle
  for (int i = 0; i < nSemi; i++) {
    double phi = TMath::Pi() / 2 + i * TMath::Pi() / (nSemi - 1);
    x[nSemi + i] = -halfLen + halfH * std::cos(phi);
    y[nSemi + i] = halfH * std::sin(phi);
  }

  TGeoXtru* xtru = new TGeoXtru(2);
  xtru->SetName(name);
  xtru->DefinePolygon(nPoints, x, y);
  xtru->DefineSection(0, -dz);
  xtru->DefineSection(1, dz);
  delete[] x;
  delete[] y;
  return xtru;
}

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
  nameM = addName("CarbonFiber");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 1.91);
    mixt->AddElement(12.01, 6, 1.0); // Carbon-only mixture
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 1);
  }
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
  nameM = addName("Aluminium");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 26.98, 13, 2.7);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kAzure - 9);
  }
}

void NA6PVerTel::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();
  // ── pixel station dimensions
  float pixChipContainerDX = 30.0f;
  float pixChipContainerDY = 30.0f;
  float pixChipContainerDz = 0.1f;
  float pixChipDX = 14.69f;
  float pixChipDY = 14.69f;
  float pixChipDz = 50e-4f;
  float pixChipOffsX = 0.29f;
  float pixChipOffsY = 0.31f;
  float carbonPlateDX = 2 * pixChipDX + 2 * pixChipOffsX;
  float carbonPlateDY = 2 * pixChipDY + 2 * pixChipOffsY;
  float carbonPlateDz = 400e-4f;

  // Al frame dimensions
  const double alFrameHX = 19.0; // base plate half-x
  const double alFrameHY = 16.0; // base plate half-y
  const double alFrameHZ = 0.4;  // base plate half-z
  const double alCutHZ = 0.5;    // cut solids slightly thicker than plate
  const double alCornerR = 1;    // corner rounding radius
  const double alSlotR = 2.5;    // stadium slot half-height
  const double alSlotL = 10.5;   // stadium slot straight half-length

  float boxDZMargin = pixChipContainerDz + 0.5f;
  float boxDZ = param.posVerTelPlaneZ[param.nVerTelPlanes - 1] - param.posVerTelPlaneZ[0] + 2 * boxDZMargin + 2 * static_cast<float>(alFrameHZ);

  auto* vtShape = new TGeoBBox("VTContainer",
                               alFrameHX + 1.f,
                               alFrameHY + 1.f,
                               boxDZ / 2);
  TGeoVolume* vtContainer = new TGeoVolume("VTContainer", vtShape,
                                           NA6PTGeoHelper::instance().getMedium(addName("Air")));
  auto* vtTransform = new TGeoCombiTrans(
    param.shiftVerTel[0],
    param.shiftVerTel[1],
    param.shiftVerTel[2] + (param.posVerTelPlaneZ[0] + boxDZ / 2 - boxDZMargin),
    NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(vtContainer, composeNonSensorVolID(0), vtTransform);

  // Pixel sensor station
  auto* pixelStationShape = new TGeoBBox("PixelStationShape",
                                         pixChipContainerDX / 2,
                                         pixChipContainerDY / 2,
                                         pixChipContainerDz / 2);
  auto* sensorShape = new TGeoBBox("SensorShape",
                                   pixChipDX / 2, pixChipDY / 2, pixChipDz / 2);

  auto* pixelStationVol = new TGeoVolume("PixelStationVol", pixelStationShape,
                                         NA6PTGeoHelper::instance().getMedium(addName("Air")));
  TGeoVolume* pixelSensor = new TGeoVolume("PixelSensor", sensorShape,
                                           NA6PTGeoHelper::instance().getMedium(addName("Silicon")));
  pixelSensor->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Silicon")));

  // Carbon-fiber plate with central beam hole 
  auto* carbonplateFullShape = new TGeoBBox("CarbonPlateFullShape",
                                            carbonPlateDX / 2, carbonPlateDY / 2, carbonPlateDz / 2);
  auto* beamHole = new TGeoBBox("CarbonPlateBeamHole", pixChipOffsX, pixChipOffsY, carbonPlateDz);
  auto* holeRemoval = new TGeoSubtraction(carbonplateFullShape, beamHole);
  auto* cbPlateWithHoleShape = new TGeoCompositeShape("CarbonPlateWithHoleShape", holeRemoval);
  TGeoVolume* cbPlate = new TGeoVolume("CarbonPlateWithHole", cbPlateWithHoleShape,
                                       NA6PTGeoHelper::instance().getMedium(addName("CarbonFiber")));
  cbPlate->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("CarbonFiber")));

  // Base plate
  auto* alBase = new TGeoBBox("AlFrameBase", alFrameHX, alFrameHY, alFrameHZ);

  // Central hole
  TGeoXtru* hC = CreateRoundedRect("AlHoleC", 6.5, 6.5, alCornerR, alCutHZ);

  // Horizontal stadium slots (top & bottom)
  TGeoXtru* hH = CreateStadium("AlHoleH", alSlotL, alSlotR, alCutHZ);
  auto* tTop = new TGeoTranslation("AlTTop", 0, 10.5, 0);
  auto* tBot = new TGeoTranslation("AlTBot", 0, -10.5, 0);
  tTop->RegisterYourself();
  tBot->RegisterYourself();

  // Vertical rounded-rectangle holes (left & right)
  TGeoXtru* hV = CreateRoundedRect("AlHoleV", 2.5, 5.5, alCornerR, alCutHZ);
  auto* tL = new TGeoTranslation("AlTL", -10.5, -0.5, 0);
  auto* tR = new TGeoTranslation("AlTR", 10.5, -0.5, 0);
  tL->RegisterYourself();
  tR->RegisterYourself();

  auto* alFrameShape = new TGeoCompositeShape("AlFrameShape",
                                              "AlFrameBase - (AlHoleC + AlHoleH:AlTTop + AlHoleH:AlTBot + AlHoleV:AlTL + AlHoleV:AlTR)");

  TGeoVolume* alFrame = new TGeoVolume("AlFrame", alFrameShape,
                                       NA6PTGeoHelper::instance().getMedium(addName("Aluminium")));
  alFrame->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Aluminium")));

  // place sensors into the station volume
  std::vector<float> alpdx{pixChipDX / 2 + pixChipOffsX, -pixChipDX / 2 + pixChipOffsX,
                           -pixChipDX / 2 - pixChipOffsX, pixChipDX / 2 - pixChipOffsX};
  std::vector<float> alpdy{pixChipDY / 2 - pixChipOffsY, pixChipDY / 2 + pixChipOffsY,
                           -pixChipDY / 2 + pixChipOffsY, -pixChipDY / 2 - pixChipOffsY};
  for (size_t ii = 0; ii < alpdx.size(); ++ii) {
    auto* sensorTransform = new TGeoTranslation(alpdx[ii], alpdy[ii], 0);
    pixelStationVol->AddNode(pixelSensor, composeSensorVolID(ii), sensorTransform);
  }

  // Carbon plate
  auto* cbTransform = new TGeoCombiTrans(0., 0.,
                                         pixChipDz / 2 + carbonPlateDz / 2,
                                         NA6PTGeoHelper::rotAroundVector(0, 0.0, 0.0, 0.0));
  pixelStationVol->AddNode(cbPlate, composeNonSensorVolID(20), cbTransform);

  // place station planes + Al frames into the VT container
  float zoffs = param.posVerTelPlaneZ[0] + boxDZ / 2 - boxDZMargin;
  float alFrameRelZ = static_cast<float>(pixChipDz / 2 + carbonPlateDz + alFrameHZ);

  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    // Pixel station
    auto* stationTransform = new TGeoCombiTrans(
      param.posVerTelPlaneX[ll],
      param.posVerTelPlaneY[ll],
      param.posVerTelPlaneZ[ll] - zoffs,
      NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
    vtContainer->AddNode(pixelStationVol, composeNonSensorVolID(ll), stationTransform);

    if (param.useAluminumPlate) {
      // Al frame behind the carbon plate
      auto* alFrameTransform = new TGeoCombiTrans(
        param.posVerTelPlaneX[ll],
        param.posVerTelPlaneY[ll],
        param.posVerTelPlaneZ[ll] - zoffs + alFrameRelZ,
        NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
      vtContainer->AddNode(alFrame, composeNonSensorVolID(ll + 40), alFrameTransform);
    }
  }
}

void NA6PVerTel::setAlignableEntries()
{
  const auto& param = NA6PLayoutParam::Instance();
  int svolCnt = 0;
  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    for (int ii = 0; ii < 4; ++ii) {
      int id = getActiveID() * 100 + svolCnt;
      std::string nm = fmt::format("VT_Lr{}_Sens{}", ll, ii);
      std::string path = fmt::format("/World_1/VTContainer_{}/PixelStationVol_{}/PixelSensor_{}", composeNonSensorVolID(0), composeNonSensorVolID(ll), composeSensorVolID(ii));
      gGeoManager->SetAlignableEntry(nm.c_str(), path.c_str(), id);
      LOGP(info, "Adding {} {} as alignable sensor {}", nm, path, id);
      svolCnt++;
    }
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
    TLorentzVector positionStop, momentumStop;
    mc->TrackMomentum(momentumStop);
    mc->TrackPosition(positionStop);
    int stationID(-1);
    mc->CurrentVolOffID(1, stationID);
    stationID = volID2NonSensID(stationID);

    int chipindex = NChipsPerStation * stationID + sensID;
    auto* p = addHit(stack->GetCurrentTrackNumber(), chipindex, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                     mTrackData.mMomentumStart.Vect(), momentumStop.Vect(), positionStop.T(),
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

NA6PVerTelHit* NA6PVerTel::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                                  float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endMom, endTime, eLoss, startStatus, endStatus);
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
