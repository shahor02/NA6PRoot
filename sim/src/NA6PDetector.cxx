// NA6PCCopyright
#include "NA6PDetector.h"
#include "NA6PModule.h"
#include "NA6PTGeoHelper.h"
#include "NA6PDipoleIP.h"
#include "NA6PTarget.h"
#include "NA6PVerTel.h"
#include "NA6PLayoutParam.h"

#include <TGeoManager.h>
#include <fairlogger/Logger.h>


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
  mModules[m->getName()] = m;
}

void NA6PDetector::createGeometry(const std::string& name)
{
  TGeoManager *geom = new TGeoManager(name.c_str(), (name + " TGeometry").c_str());
  createCommonMaterials();
  TGeoVolume *world = geom->MakeBox("World", NA6PTGeoHelper::instance().getMedium("Air"), 1000, 1000, 1000);
  geom->SetTopVolume(world);
  //
  // build modules
  addModule(new NA6PDipoleIP());
  addModule(new NA6PTarget());
  addModule(new NA6PVerTel());
  
  for (auto m : mModules) {
    m.second->createMaterials();
    m.second->createGeometry(world);
  }
  
  // Close the geometry
  geom->CloseGeometry();

  // Export to file
  geom->Export("geometry.root");
}
