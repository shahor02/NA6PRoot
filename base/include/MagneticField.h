// NA6PCCopyright
#ifndef _MAGFIELD_H
#define _MAGFIELD_H

#include "MagneticFieldRegion.h"
#include <TGeoGlobalMagField.h>
#include <TVirtualMagField.h>

class MagneticField : public TVirtualMagField
{
 public:

  MagneticField() : TVirtualMagField("NA6PMagField") {}
  ~MagneticField() override;
  void Field(const Double_t *x, Double_t *B) override
  {
    getField(x, B);
  }

  bool isInitDone() const { return mInitDone; }
  void loadFlukaField();

  template<typename T = float>
  bool getField(const T* xyz, T* bxbybz) const
  {
    bxbybz[0] = bxbybz[1] = bxbybz[2] = 0;
    return mDipoleIP.getField(xyz, bxbybz) || mDipoleMS.getField(xyz, bxbybz); // we know that the fields are not overlapping
  }
  void setAsGlobalField();
  
 private:
  
  MagneticFieldRegion mDipoleIP{};
  MagneticFieldRegion mDipoleMS{};
  bool mInitDone = false;
  
  ClassDefOverride(MagneticField, 1);
};

#endif
