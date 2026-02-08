// NA6PCCopyright
#include "NA6PMagnetsGDML.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PDipoleVT.h"
#include "NA6PDipoleMS.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <cstdlib>

std::string NA6PMagnetsGDML::getGDMLPath() const
{
  const char* gdmlPathEnv = std::getenv("NA6P_BOTHMAGNETS_GDML");
  // Default to BothMagnets.gdml in user's Alice directory
  // TODO: should be moved to NA6PRoot/data/ directory for production
  return gdmlPathEnv ? gdmlPathEnv : "/home/access/Alice/NA6PRoot/data/BothMagnets.gdml";
}

void NA6PMagnetsGDML::createMaterials()
{
  // For GDML-based magnets we need to ensure that the analytic magnet
  // materials (Iron_DipoleVT, Copper_DipoleVT, Iron_DipoleMS, Copper_DipoleMS)
  // exist so we can map GDML volumes to them for consistency.
  {
    NA6PDipoleVT vtDummy;
    vtDummy.createMaterials();
    NA6PDipoleMS msDummy;
    msDummy.createMaterials();
  }
}

void NA6PMagnetsGDML::createGeometry(TGeoVolume* world)
{
  // The GDML import and magnet setup is handled by NA6PDetector::createGeometry
  // which calls the helper methods (assignMediaToGDMLVolumes, alignMagnetsToLayout,
  // harmoniseVolumeNames) directly. This method is left empty because the GDML
  // file already contains the magnet geometry.
}

void NA6PMagnetsGDML::assignMediaToGDMLVolumes()
{
  // Assign proper TGeoMedium to volumes imported from GDML, which initially
  // have no medium or dummy medium. This is needed for Geant4 VMC to work.
  auto& helper = NA6PTGeoHelper::instance();
  auto* geom = gGeoManager;

  TGeoMedium* airMed = helper.getMedium("Air", /*fatalIfMissing=*/false);
  if (airMed) {
    auto assignAirIfNeeded = [&](const char* vname) {
      if (auto* v = geom->FindVolumeFast(vname)) {
        if (v->GetMedium() != airMed) {
          LOGP(info, "Assigning Air medium to volume {}", v->GetName());
          v->SetMedium(airMed);
        }
      }
    };
    assignAirIfNeeded("worldVOL");
    assignAirIfNeeded("World");
    assignAirIfNeeded("MNP33");
    assignAirIfNeeded("MEP48");
  }

  // Get GDML materials (fallback if analytic ones are not available)
  auto getMedIfExists = [geom](const char* name) -> TGeoMedium* {
    return geom->GetMedium(name);
  };

  TGeoMedium* feMed = getMedIfExists("G4_Fe");
  if (!feMed) {
    feMed = getMedIfExists("G4_Fe0x39eca690");
  }
  TGeoMedium* cuMed = getMedIfExists("G4_Cu");
  if (!cuMed) {
    cuMed = getMedIfExists("G4_Cu0x39e8a110");
  }

  // Prefer analytic magnet media for consistency
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

  // First magnet (DipoleVT / MEP48): gap, yokes, poles, coils
  assignMediumIfFound("gap",      feVT ? feVT : feMed);
  assignMediumIfFound("yokeUp",   feVT ? feVT : feMed);
  assignMediumIfFound("yokeDown", feVT ? feVT : feMed);
  assignMediumIfFound("poleUp",   feVT ? feVT : feMed);
  assignMediumIfFound("poleDown", feVT ? feVT : feMed);
  assignMediumIfFound("coilUp",   cuVT ? cuVT : cuMed);
  assignMediumIfFound("coilDown", cuVT ? cuVT : cuMed);

  // Second magnet (DipoleMS / MNP33): yoke, coils
  assignMediumIfFound("yoke",       feMS ? feMS : feMed);
  assignMediumIfFound("coil_big",   cuMS ? cuMS : cuMed);
  assignMediumIfFound("coil_small", cuMS ? cuMS : cuMed);
}

void NA6PMagnetsGDML::alignMagnetsToLayout(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();
  auto* geom = gGeoManager;

  // --- First magnet (MEP48) : align to posDipIP[...] ---
  {
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

  // --- Second magnet (MNP33) : align to posDipMS[...] ---
  {
    const double targetX = param.posDipMS[0];
    const double targetY = param.posDipMS[1];
    const double targetZ = param.posDipMS[2];

    TGeoIterator it(world);
    TGeoNode* node = nullptr;
    TGeoNode* mnp33Node = nullptr;

    while ((node = it())) {
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
}

void NA6PMagnetsGDML::harmoniseVolumeNames()
{
  // Rename GDML volumes to match analytic magnet names:
  //   MEP48 (first magnet) -> DipoleVT
  //   MNP33 (second magnet) -> DipoleMS
  auto* geom = gGeoManager;
  if (auto* v = geom->FindVolumeFast("MEP48")) {
    LOGP(info, "Renaming GDML volume MEP48 to DipoleVT for consistency with analytic geometry");
    v->SetName("DipoleVT");
  }
  if (auto* v = geom->FindVolumeFast("MNP33")) {
    LOGP(info, "Renaming GDML volume MNP33 to DipoleMS for consistency with analytic geometry");
    v->SetName("DipoleMS");
  }
}
