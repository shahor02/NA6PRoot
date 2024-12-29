// NA6PCCopyright

#include "MagneticField.h"
#include "NA6PLayoutParam.h"
#include <fairlogger/Logger.h>
#include <TVirtualMagField.h>
#include <TSystem.h>

MagneticField::~MagneticField()
{
  if (TGeoGlobalMagField::Instance()->GetField() == this) {
    LOGP(error, "This field {} was set as global field, cannot be deleted explicitly", GetName());
  }
}

void MagneticField::loadFlukaField()
{
  const auto& param = NA6PLayoutParam::Instance();
  mDipoleVT.loadFlukaField(gSystem->ExpandPathName(param.flukaInpDipIP.c_str()));  
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
