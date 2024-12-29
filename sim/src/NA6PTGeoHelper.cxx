// NA6PCCopyright
#include "NA6PTGeoHelper.h"
#include <TMath.h>

void NA6PTGeoHelper::addMedium(const std::string& medName, const std::string& matName, Color_t col )
{
  const auto& matN = matName.empty() ? medName : matName;
  if (mMedPool.find(medName) != mMedPool.end()) {
    throw std::runtime_error(Form("Medium %s was already created", medName.c_str()));
  }
  if (mMatPool.find(matN) == mMatPool.end()) {
    throw std::runtime_error(Form("Material %s does not exist", matN.c_str()));
  }
  mMedPool[medName] = new TGeoMedium(matN.c_str(), mMedPool.size(), mMatPool[matN]);
  mColorPool[medName] = col;
}

TGeoMedium* NA6PTGeoHelper::getMedium(const std::string& medName) const
{
  auto m = mMedPool.find(medName);
  if (m == mMedPool.end()) {
    throw std::runtime_error(Form("Medium %s was not created", medName.c_str()));
  }
  return m->second;
}

Color_t NA6PTGeoHelper::getMediumColor(const std::string& medName) const
{
  auto m = mColorPool.find(medName);
  if (m == mColorPool.end()) {
    throw std::runtime_error(Form("Medium %s was not created", medName.c_str()));
  }
  return m->second;
}

TGeoRotation* NA6PTGeoHelper::rotAroundVector(float uX, float uY, float uZ, float ddelta)
{
  ddelta *= TMath::DegToRad();
  double sinDelta = std::sin(ddelta), cosDelta = std::cos(ddelta);
  double oneMinusCosDelta = 1.0 - cosDelta;
  double rmat[9] = {
    oneMinusCosDelta * uX * uX + cosDelta,
    oneMinusCosDelta * uX * uY - sinDelta * uZ,
    oneMinusCosDelta * uX * uZ + sinDelta * uY,
    oneMinusCosDelta * uY * uX + sinDelta * uZ,
    oneMinusCosDelta * uY * uY + cosDelta,
    oneMinusCosDelta * uY * uZ - sinDelta * uX,
    oneMinusCosDelta * uZ * uX - sinDelta * uY,
    oneMinusCosDelta * uZ * uY + sinDelta * uX,
    oneMinusCosDelta * uZ * uZ + cosDelta
  };
  auto rot = new TGeoRotation();
  rot->SetMatrix(rmat);
  return rot;
}

