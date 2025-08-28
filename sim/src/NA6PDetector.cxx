// NA6PCCopyright
#include "NA6PDetector.h"
#include "NA6PModule.h"
#include "NA6PTGeoHelper.h"
#include "NA6PDipoleVT.h"
#include "NA6PDipoleMS.h"
#include "NA6PTarget.h"
#include "NA6PVerTel.h"
#include "NA6PAbsorber.h"
#include "NA6PMuonSpec.h"
#include "NA6PMuonSpecModular.h"
#include "NA6PLayoutParam.h"
#include "StringUtils.h"
#include <TGeoManager.h>
#include <fairlogger/Logger.h>

NA6PDetector::NA6PDetector()
{
  // just declare modules
  addModule(new NA6PDipoleVT());
  addModule(new NA6PTarget());
  addModule(new NA6PVerTel());
  addModule(new NA6PAbsorber());
  addModule(new NA6PDipoleMS());
  //addModule(new NA6PMuonSpec());
  addModule(new NA6PMuonSpecModular());
}

void NA6PDetector::createCommonMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Air") == matPool.end()) {
    auto mixt = new TGeoMixture("Air", 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool["Air"] = mixt;
    NA6PTGeoHelper::instance().addMedium("Air");
  }
}

void NA6PDetector::addModule(NA6PModule* m)
{
  if (!m) {
    LOGP(fatal, "null pointer is provided for a module");
  }
  if (getModule(m->getName())) {
    LOGP(fatal, "module {} was already defined", m->getName());
  }
  if (m->isActive()) {
    for (const auto mo : mModulesVec) {
      if (mo->getActiveID() == m->getActiveID()) {
        LOGP(fatal, "Redefinition of ActiveModuleID:{} {}, previosly defined module: {}", m->getActiveID(), m->getName(), mo->getName());
      }
    }
  }
  mModulesMap[m->getName()] = m;
  m->setID((int)mModulesVec.size() + 1); // ID should be counted from 1!!!
  mModulesVec.push_back(m);
  if (m->isActive()) {
    mActiveModulesVec.push_back(m);
  }
  LOGP(info, "Adding module:{} {} ActiveID:{}", m->getID(), m->getName(), m->getActiveID());
}

void NA6PDetector::createGeometry(const std::string& name)
{
  TGeoManager* geom = new TGeoManager(name.c_str(), (name + " TGeometry").c_str());
  createCommonMaterials();
  TGeoVolume* world = geom->MakeBox("World", NA6PTGeoHelper::instance().getMedium("Air"), 1000, 1000, 1000);
  geom->SetTopVolume(world);
  //
  // build modules
  for (auto m : mModulesVec) {
    m->createMaterials();
    m->createGeometry(world);
  }

  // Close the geometry
  geom->CloseGeometry();
  for (auto m : mModulesVec) {
    m->setAlignableEntries();
  }
  // Export to file
  geom->Export("geometry.root");
}

void NA6PDetector::setVerbosity(int v)
{
  mVerbosity = v;
  for (auto m : mModulesVec) {
    m->setVerbosity(v);
  }
}

void NA6PDetector::createHitsOutput(const std::string& dir)
{
  for (auto* det : mActiveModulesVec) {
    det->createHitsOutput(dir);
  }
}

void NA6PDetector::writeHits(const std::vector<int>& remapping)
{
  for (auto* det : mActiveModulesVec) {
    det->writeHits(remapping);
  }
}

void NA6PDetector::closeHitsOutput()
{
  for (auto* det : mActiveModulesVec) {
    det->closeHitsOutput();
  }
}
