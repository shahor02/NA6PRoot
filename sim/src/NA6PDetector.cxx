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
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <cstdlib>

namespace
{
// Global switch: false = analytic C++ magnets (NA6PDipoleVT/MS),
// true = import magnets from GDML (BothMagnets.gdml).
bool gUseGDMLMagnets = false;
} // namespace

void NA6PDetector::setUseGDMLMagnets(bool v)
{
  gUseGDMLMagnets = v;
}

NA6PDetector::NA6PDetector()
{
  // just declare modules
  if (!gUseGDMLMagnets) {
    // Analytic magnets from C++ (original behaviour)
    addModule(new NA6PDipoleVT());
  }

  addModule(new NA6PTarget());
  addModule(new NA6PVerTel());
  addModule(new NA6PAbsorber());

  if (!gUseGDMLMagnets) {
    addModule(new NA6PDipoleMS());
  }
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
  if (!gUseGDMLMagnets) {
    // === Analytic geometry: build magnets in C++ (NA6PDipoleVT/MS) ===
    TGeoManager* geom =
      new TGeoManager(name.c_str(), (name + " TGeometry").c_str());

    // Make sure common materials (e.g. Air) exist
    createCommonMaterials();

    // World box in cm (half-sizes 1000x1000x1000)
    TGeoVolume* world =
      geom->MakeBox("World",
                    NA6PTGeoHelper::instance().getMedium("Air"),
                    1000., 1000., 1000.);
    geom->SetTopVolume(world);

    // build modules (including NA6PDipoleVT/MS if they were added)
    for (auto m : mModulesVec) {
      m->createMaterials();
      m->createGeometry(world);
    }

    // Recolour volumes: all Iron* media blue, all Copper* media yellow.
    {
      auto* vols = geom->GetListOfVolumes();
      if (vols) {
        for (int i = 0; i < vols->GetEntriesFast(); ++i) {
          auto* v = static_cast<TGeoVolume*>(vols->At(i));
          if (!v) {
            continue;
          }
          auto* med = v->GetMedium();
          if (!med) {
            continue;
          }
          const TGeoMaterial* mat = med->GetMaterial();
          const char* mname = nullptr;
          if (mat && mat->GetName()) {
            mname = mat->GetName();
          } else if (med->GetName()) {
            mname = med->GetName();
          }
          if (!mname) {
            continue;
          }
          std::string nm = mname;
          if (nm.find("Iron") != std::string::npos || nm.find("Fe") != std::string::npos) {
            v->SetLineColor(kGray);
          } else if (nm.find("Copper") != std::string::npos || nm.find("Cu") != std::string::npos) {
            v->SetLineColor(kRed);
          }
        }
      }
    }

    // Close the geometry
    geom->CloseGeometry();
    for (auto m : mModulesVec) {
      m->setAlignableEntries();
    }
    // Export to file
    geom->Export("geometry.root");
    return;
  }

  // === GDML-based magnets: import BothMagnets.gdml and attach remaining modules ===

  // Path to the GDML file can be overridden via NA6P_BOTHMAGNETS_GDML.
  const char* gdmlPathEnv = std::getenv("NA6P_BOTHMAGNETS_GDML");
  std::string gdmlPath = gdmlPathEnv ? gdmlPathEnv : "/home/access/Alice/NA6PRoot/data/BothMagnets.gdml";

  LOGP(info, "Importing geometry from GDML file: {}", gdmlPath);
  TGeoManager::Import(gdmlPath.c_str());

  TGeoManager* geom = gGeoManager;
  if (!geom) {
    LOGP(fatal, "Failed to create TGeoManager from GDML file {}", gdmlPath);
  }

  // Ensure common materials (e.g. Air) exist for modules that use NA6PTGeoHelper.
  createCommonMaterials();

  TGeoVolume* world = geom->GetTopVolume();
  if (!world) {
    LOGP(fatal, "No top volume found after importing GDML file {}", gdmlPath);
  }

  {
    const auto& param = NA6PLayoutParam::Instance();
    const double shiftX = param.posDipIP[0];
    const double shiftY = param.posDipIP[1];
    const double shiftZ = param.posDipIP[2];

    TGeoIterator it(world);
    TGeoNode* node = nullptr;
    TGeoNode* mep48Node = nullptr;

    while ((node = it())) {
      const char* nname = node->GetName();
      const char* vname = node->GetVolume()->GetName();
      if ((strcmp(nname, "PV_MEP48") == 0) || (strcmp(vname, "MEP48") == 0)) {
        mep48Node = node;
        break;
      }
    }

    if (mep48Node) {
      const Double_t* tr = it.GetCurrentMatrix()->GetTranslation();
      // Apply layout shift relative to the baseline GDML placement.
      Double_t newTr[3] = {tr[0] + shiftX, tr[1] + shiftY, tr[2] + shiftZ};

      auto* origMat = mep48Node->GetMatrix();
      auto* newMat = new TGeoHMatrix(*origMat);
      newMat->SetTranslation(newTr);

      if (auto* nodeMat = dynamic_cast<TGeoNodeMatrix*>(mep48Node)) {
        nodeMat->SetMatrix(newMat);
        LOGP(info, "Shifted first magnet (PV_MEP48) by d(posDipIP) = (%g,%g,%g) cm",
             shiftX, shiftY, shiftZ);
      } else {
        LOGP(warn, "Found PV_MEP48 / MEP48, but it is not a TGeoNodeMatrix; cannot reset matrix");
        delete newMat;
      }
    } else {
      LOGP(warn, "Could not find PV_MEP48 / MEP48 node to align first magnet to posDipIP");
    }
  }

  // --- Second magnet (MNP33) : PV_MNP33 / MNP33 ---
  {
    const auto& param = NA6PLayoutParam::Instance();
    const double targetX = param.posDipMS[0];
    const double targetY = param.posDipMS[1];
    const double targetZ = param.posDipMS[2];

    TGeoIterator it(world);
    TGeoNode* node = nullptr;
    TGeoNode* mnp33Node = nullptr;

    while ((node = it())) {
      // Heuristics: either node name "PV_MNP33" (from GDML physvol),
      // or volume name "MNP33" (local world of the second magnet).
      const char* nname = node->GetName();
      const char* vname = node->GetVolume()->GetName();
      if ((strcmp(nname, "PV_MNP33") == 0) || (strcmp(vname, "MNP33") == 0)) {
        mnp33Node = node;
        break;
      }
    }

    if (mnp33Node) {
      const Double_t* tr = it.GetCurrentMatrix()->GetTranslation();
      double dx = targetX - tr[0];
      double dy = targetY - tr[1];
      double dz = targetZ - tr[2];

      // Create a new matrix with updated translation and attach it to the node,
      // if it is a TGeoNodeMatrix (only that subclass owns a mutable matrix).
      auto* origMat = mnp33Node->GetMatrix();
      auto* newMat = new TGeoHMatrix(*origMat);
      Double_t newTr[3] = {tr[0] + dx, tr[1] + dy, tr[2] + dz};
      newMat->SetTranslation(newTr);

      if (auto* nodeMat = dynamic_cast<TGeoNodeMatrix*>(mnp33Node)) {
        nodeMat->SetMatrix(newMat);
        LOGP(info, "Shifted second magnet (PV_MNP33) by dx=%g dy=%g dz=%g cm to align with posDipMS",
             dx, dy, dz);
      } else {
        LOGP(warn, "Found PV_MNP33 / MNP33, but it is not a TGeoNodeMatrix; cannot reset matrix");
        delete newMat;
      }
    } else {
      LOGP(warn, "Could not find PV_MNP33 / MNP33 node to align second magnet to posDipMS");
    }
  }

  {
    auto& helper = NA6PTGeoHelper::instance();
    TGeoMedium* airMed = helper.getMedium("Air", /*fatalIfMissing=*/true);
    if (airMed) {
      if (world->GetMedium() != airMed) {
        LOGP(info, "Assigning Air medium to top volume {}", world->GetName());
        world->SetMedium(airMed);
      }
      if (auto v = geom->FindVolumeFast("World")) {
        if (v->GetMedium() != airMed) {
          LOGP(info, "Assigning Air medium to volume {}", v->GetName());
          v->SetMedium(airMed);
        }
      }

      if (auto v = geom->FindVolumeFast("MNP33")) {
        if (v->GetMedium() != airMed) {
          LOGP(info, "Assigning Air medium to volume {}", v->GetName());
          v->SetMedium(airMed);
        }
      }
      // First magnet local world (assembly MEP48) in GDML
      if (auto v = geom->FindVolumeFast("MEP48")) {
        if (v->GetMedium() != airMed) {
          LOGP(info, "Assigning Air medium to volume {}", v->GetName());
          v->SetMedium(airMed);
        }
      }
    }

    // Try to reuse media created from GDML materials for iron and copper.
    auto getMedIfExists = [geom](const char* name) -> TGeoMedium* {
      auto* m = geom->GetMedium(name);
      return m;
    };

    TGeoMedium* feMed = getMedIfExists("G4_Fe");
    if (!feMed) {
      feMed = getMedIfExists("G4_Fe0x39eca690");
    }
    TGeoMedium* cuMed = getMedIfExists("G4_Cu");
    if (!cuMed) {
      cuMed = getMedIfExists("G4_Cu0x39e8a110");
    }

    // assign names Iron_DipoleVT / Copper_DipoleVT / Iron_DipoleMS /
    // Copper_DipoleMS to GDML magnet parts.
    {
      NA6PDipoleVT vtDummy;
      vtDummy.createMaterials();
      NA6PDipoleMS msDummy;
      msDummy.createMaterials();
    }

    TGeoMedium* feVT = helper.getMedium("Iron_DipoleVT", /*fatalIfMissing=*/false);
    TGeoMedium* cuVT = helper.getMedium("Copper_DipoleVT", /*fatalIfMissing=*/false);
    TGeoMedium* feMS = helper.getMedium("Iron_DipoleMS", /*fatalIfMissing=*/false);
    TGeoMedium* cuMS = helper.getMedium("Copper_DipoleMS", /*fatalIfMissing=*/false);

    auto assignMediumIfFound = [&](const char* volName, TGeoMedium* med) {
      if (!med) {
        return;
      }
      if (auto* v = geom->FindVolumeFast(volName)) {
        LOGP(info, "Assigning medium {} to volume {}", med->GetName(), v->GetName());
        v->SetMedium(med);
      }
    };

    // Magnet volumes in BothMagnets.gdml

    assignMediumIfFound("gap",      feVT ? feVT : feMed);
    assignMediumIfFound("yokeUp",   feVT ? feVT : feMed);
    assignMediumIfFound("yokeDown", feVT ? feVT : feMed);
    assignMediumIfFound("poleUp",   feVT ? feVT : feMed);
    assignMediumIfFound("poleDown", feVT ? feVT : feMed);
    assignMediumIfFound("coilUp",   cuVT ? cuVT : cuMed);
    assignMediumIfFound("coilDown", cuVT ? cuVT : cuMed);
    assignMediumIfFound("yoke",       feMS ? feMS : feMed);
    assignMediumIfFound("coil_big",   cuMS ? cuMS : cuMed);
    assignMediumIfFound("coil_small", cuMS ? cuMS : cuMed);
  }

  {
    if (auto* v = geom->FindVolumeFast("MEP48")) {
      LOGP(info, "Renaming GDML volume MEP48 to DipoleVT for consistency with analytic geometry");
      v->SetName("DipoleVT");
    }
    if (auto* v = geom->FindVolumeFast("MNP33")) {
      LOGP(info, "Renaming GDML volume MNP33 to DipoleMS for consistency with analytic geometry");
      v->SetName("DipoleMS");
    }
  }

  // build remaining modules (magnets are already in the GDML)
  for (auto m : mModulesVec) {
    m->createMaterials();
    m->createGeometry(world);
  }

  // Recolour volumes: all Iron* media blue, all Copper* media yellow.
  {
    auto* vols = geom->GetListOfVolumes();
    if (vols) {
      for (int i = 0; i < vols->GetEntriesFast(); ++i) {
        auto* v = static_cast<TGeoVolume*>(vols->At(i));
        if (!v) {
          continue;
        }
        auto* med = v->GetMedium();
        if (!med) {
          continue;
        }
        const TGeoMaterial* mat = med->GetMaterial();
        const char* mname = nullptr;
        if (mat && mat->GetName()) {
          mname = mat->GetName();
        } else if (med->GetName()) {
          mname = med->GetName();
        }
        if (!mname) {
          continue;
        }
        std::string nm = mname;
        if (nm.find("Iron") != std::string::npos || nm.find("Fe") != std::string::npos) {
          v->SetLineColor(kGray);
        } else if (nm.find("Copper") != std::string::npos || nm.find("Cu") != std::string::npos) {
          v->SetLineColor(kRed);
        }
      }
    }
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
