// NA6PCCopyright

#include "NA6PMuonSpec.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"

#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>

void NA6PMuonSpec::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Silicon") == matPool.end()) {
    matPool["Silicon"] = new TGeoMaterial("Silicon", 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium("Silicon");
  }
}

void NA6PMuonSpec::createGeometry(TGeoVolume *world)
{
  const auto& param = NA6PLayoutParam::Instance();
  
  createMaterials();

  for (int ist=0;ist<param.nMSPlanes;ist++) {
    TGeoShape* station = nullptr;
    auto stnm = fmt::format("MS{}", ist);
    if (param.dimYMSPlane[ist] > 0) {
      station = new TGeoBBox((stnm+"SH").c_str(), param.dimXMSPlane[ist]/2, param.dimYMSPlane[ist]/2, param.thicknessMSPlane[ist]/2);
    } else {
      station = new TGeoTube((stnm+"SH").c_str(), 0.0, param.dimXMSPlane[ist]/2, param.thicknessMSPlane[ist]/2);
    }
    if (param.dimXMSPlaneHole[ist]>0) {
      TGeoShape* hole = nullptr;
      if (param.dimYMSPlaneHole[ist]>0) {
	hole = new TGeoBBox((stnm + "HL").c_str(), param.dimXMSPlaneHole[ist]/2, param.dimYMSPlaneHole[ist]/2, param.thicknessMSPlane[ist]);
      } else {
	hole = new TGeoTube((stnm + "HL").c_str(), 0.0, param.dimXMSPlaneHole[ist]/2, param.thicknessMSPlane[ist]);
      }
      station = new TGeoCompositeShape((stnm + "_HS").c_str(), new TGeoSubtraction(station, hole));   // station with plug hole
    }
    auto stationVol = new TGeoVolume(stnm.c_str(), station, NA6PTGeoHelper::instance().getMedium(param.medMSPlane[ist]));
    world->AddNode(stationVol, ist + 1, new TGeoTranslation(param.shiftMS[0] + param.posMSPlaneX[ist],
							    param.shiftMS[1] + param.posMSPlaneY[ist],
							    param.shiftMS[2] + param.posMSPlaneZ[ist]));
  }

}
