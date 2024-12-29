// NA6PCCopyright


#ifndef NA6P_LAYOUT_PARAMS_H_
#define NA6P_LAYOUT_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PLayoutParam : public na6p::conf::ConfigurableParamHelper<NA6PLayoutParam>
{
  static constexpr int MaxVTPlanes = 10;
  static constexpr int MaxTargers = 10;

  // VerTel Dipole 
  float posDipIP[3] = {0.f, 0.f, 0.f};  // Vertex dipole position
  std::string flukaInpDipIP = "$(NA6PROOT_ROOT)/share/data/MEP48_field_map.inp";
  
  // MS Dipole 
  float posDipMS[3] = {0.f, 0.f, 430.f};  // Vertex dipole position
  std::string flukaInpDipMS = "$(NA6PROOT_ROOT)/share/data/MNP33_field_map.inp";
  
  // Target
  int nTargets = 5;
  float shiftTargets[3] = {0.f, 0.f, 0.f}; // Targer box global shift, added to posTargetX,Y,Z
  float posTargetX[MaxVTPlanes] = {};
  float posTargetY[MaxVTPlanes] = {};
  float posTargetZ[MaxVTPlanes] = {-2.f, -1.f, 0.f, 1.f, 2.f};
  float thicknessTarget[MaxVTPlanes] = {0.15f, 0.15f, 0.15f, 0.15f, 0.15f};
  float radTarget[MaxVTPlanes] = {0.3f, 0.1f, 0.1f, 0.1f, 0.1f};
  std::string medTarget[MaxVTPlanes] = {"Lead", "Lead", "Lead", "Lead", "Lead"};

  // VerTel
  int nVerTelPlanes = 5;                // number of stations
  float shiftVerTel[3] = {0.f, 0.f, 0.f}; // VerTel box global shift, added to posVerTelPlaneX,Y,Z
  float posVerTelPlaneX[MaxVTPlanes] = {};
  float posVerTelPlaneY[MaxVTPlanes] = {};
  float posVerTelPlaneZ[MaxVTPlanes] = {7.1175f, 15.1175f, 20.1175f, 25.1175f, 38.1175f};
  
  
  O2ParamDef(NA6PLayoutParam, "layout");
};

#endif
