// NA6PCCopyright

#include "MagneticField.h"
#include "NA6PLayoutParam.h"
#include <fairlogger/Logger.h>
#include <TVirtualMagField.h>
#include <TSystem.h>
#include "NA6PLayoutParam.h"

MagneticField::~MagneticField()
{ 
  if (TGeoGlobalMagField::Instance()->GetField() == this) {
    LOGP(error, "This field {} was set as global field, cannot be deleted explicitly", GetName());
  }
}

void MagneticField::loadFlukaField()
{
  const auto& param = NA6PLayoutParam::Instance();
  mDipoleVT.loadFlukaField(gSystem->ExpandPathName(param.flukaInpDipVT.c_str()));
  mDipoleMS.loadFlukaField(gSystem->ExpandPathName(param.flukaInpDipMS.c_str()));
  mDipoleVT.setRefPosition(param.posDipIP[0], param.posDipIP[1], param.posDipIP[2]);
  mDipoleMS.setRefPosition(param.posDipMS[0], param.posDipMS[1], param.posDipMS[2]);
  mInitDone = true;
}

void MagneticField::setAsGlobalField()
{
  if (!isInitDone()) {
    LOGP(fatal, "The magnetic field was not initialized yet");
  }
  TGeoGlobalMagField::Instance()->SetField(this);
}

void MagneticField::dumpMagFieldMap(double xmin, double xmax,
                                    double ymin, double ymax,
                                    double zmin, double zmax,
                                    double step, std::string fieldMapFile,
                                    std::string rotatedFieldMapFile)
{
  std::ofstream fieldMap(fieldMapFile.c_str());
  std::ofstream rotatedFieldMap(rotatedFieldMapFile.c_str());

  Double_t x[3];
  Double_t b[3];

  for (double ix = xmin; ix <= xmax; ix += step) {
    for (double iy = ymin; iy <= ymax; iy += step) {
      for (double iz = zmin; iz <= zmax; iz += step) {
        x[0] = ix;
        x[1] = iy;
        x[2] = iz;

        Field(x, b);
        
        // ACTS uses mm instead of cm
        fieldMap << ix*10 << " " << iy*10 << " " << iz*10 << " "
            << b[0]/mScale2Unit << " " << b[1]/mScale2Unit << " " << b[2]/mScale2Unit << "\n";
      }
    }
  }


  for (double ix = xmin; ix <= xmax; ix += step) {
    for (double iz = ymin; iz <= ymax; iz += step) {
      for (double iy = zmin; iy <= zmax; iy += step) {
        x[0] = ix;
        x[1] = iy;
        x[2] = iz;

        Field(x, b);
        
        // ACTS uses mm instead of cm
        rotatedFieldMap << ix*10 << " " << iz*10 << " " << -iy*10 << " "
            << b[0]/mScale2Unit << " " << b[1]/mScale2Unit << " " << b[2]/mScale2Unit << "\n";
      }
    }
  }

  fieldMap.close();
  rotatedFieldMap.close();
}

void MagneticField::dumpMagFieldFromConfig()
{
  const auto& param = NA6PLayoutParam::Instance();

  dumpMagFieldMap(param.xMinFieldDump, param.xMaxFieldDump,
                  param.yMinFieldDump, param.yMaxFieldDump,
                  param.zMinFieldDump, param.zMaxFieldDump,
                  param.stepFieldDump, param.fieldMap,
                  param.rotatedFieldMap);
}