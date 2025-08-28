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

void MagneticField::loadField()
{
    const auto& param = NA6PLayoutParam::Instance();
    auto loadOne = [](auto& dipole, const std::string& path, const float pos[3], bool flipSign) {
        std::string fname = gSystem->ExpandPathName(path.c_str());
        if (fname.find(".inp") != std::string::npos) {
            // FLUKA format
            dipole.loadFlukaField(fname, flipSign);
        } else {
            // Alternative format: OPERA3D or Ansys 
            dipole.loadOPERA3DField(fname, flipSign);
        }
        dipole.setRefPosition(pos[0], pos[1], pos[2]);
    };
    loadOne(mDipoleVT, param.flukaInpDipVT, param.posDipIP, param.flipSignVT);
    loadOne(mDipoleMS, param.flukaInpDipMS, param.posDipMS, param.flipSignMS);
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