// NA6PCCopyright


#ifndef NA6P_LAYOUT_PARAMS_H_
#define NA6P_LAYOUT_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PLayoutParam : public na6p::conf::ConfigurableParamHelper<NA6PLayoutParam>
{
  static constexpr int MaxVTPlanes = 7;
  static constexpr int MaxMSPlanes = 7;
  static constexpr int MaxTargers = 5;
  static constexpr int MaxAbsorberSlices = 7;
 
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

  // Muon Stations
  int nMSPlanes = 6;                   // number of stations
  float shiftMS[3] = {0.f, 0.f, 0.f};  // MS global shift, added to posMSPlaneX,Y,Z
  float posMSPlaneX[MaxMSPlanes] = {};
  float posMSPlaneY[MaxMSPlanes] = {};
  float posMSPlaneZ[MaxMSPlanes] = {300.f, 360.f, 530.f, 590.f, 810.f, 850.f};
  float thicknessMSPlane[MaxMSPlanes] = {0.06f, 0.06f, 0.06f, 0.06f, 0.06f, 0.06f};
  float dimXMSPlane[MaxMSPlanes] = {162.f, 240.f, 240.f, 318.f, 417.f, 459.f};   // size in X of rectangular plane. or diameter if dimYMSPlane is 0
  float dimYMSPlane[MaxMSPlanes] = {162.f, 240.f, 240.f, 318.f, 417.f, 458.f};   // size in Y if >0
  float dimXMSPlaneHole[MaxMSPlanes] = {5.f, 5.f, 5.f, 5.f, 5.f, 5.f};   // size in X of rectangular hole. or diameter if dimYMSPlaneHole is 0
  float dimYMSPlaneHole[MaxMSPlanes] = {5.f, 5.f, 5.f, 5.f, 5.f, 5.f};   // size in Y if >0  
  std::string medMSPlane[MaxMSPlanes] = {"Silicon", "Silicon", "Silicon", "Silicon", "Silicon", "Silicon"};
  
  // Absorbers
  int nAbsorberSlices = 6;
  float posZStartAbsorber = 82.25f;   // start of the absorber
  float posMuonWall = 610.f;          // start of the muon wall (if!=0 : last absorber slice)
  float thicknessAbsorber[MaxAbsorberSlices] = {39.5f, 35.0f, 30.5, 65.f, 65.f, 180.f};
  float dimXAbsorber[MaxAbsorberSlices] = {52.f, 80.f, 120.f, 170.f, 220.f, 550.f};   // size in X of rectangular abs. or diameter if dimYAbsorber is 0
  float dimYAbsorber[MaxAbsorberSlices] = {52.f, 80.f, 120.f, 170.f, 220.f, 550.f};   // size in Y if >0
  float radPlug[MaxAbsorberSlices] = {2.47f, 4.0f, 5.48f, 11.f,  11.f, 0.f};       // radius of the plug
  std::string medAbsorber[MaxAbsorberSlices] = {"BeO", "BeO", "BeO", "Graphite", "Graphite", "Graphite"};  
  std::string medAbsorberPlug[MaxAbsorberSlices] = {"Tungsten", "Tungsten", "Tungsten", "Tungsten", "Tungsten", ""};
  
  O2ParamDef(NA6PLayoutParam, "layout");
};

#endif
