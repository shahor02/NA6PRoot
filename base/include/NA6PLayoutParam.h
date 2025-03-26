// NA6PCCopyright

#ifndef NA6P_LAYOUT_PARAMS_H_
#define NA6P_LAYOUT_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PLayoutParam : public na6p::conf::ConfigurableParamHelper<NA6PLayoutParam> {
  static constexpr int MaxVTPlanes = 7;
  static constexpr int MaxMSPlanes = 7;
  static constexpr int MaxTargets = 5;
  static constexpr int MaxAbsorberSlices = 10;

  // VerTel Dipole
  float posDipIP[3] = {0.f, 0.f, 0.f}; // Vertex dipole position
  std::string flukaInpDipVT = "$(NA6PROOT_ROOT)/share/data/MEP48_field_map.inp";

  // MS Dipole
  float posDipMS[3] = {0.f, 0.f, 430.f}; // Vertex dipole position
  std::string flukaInpDipMS = "$(NA6PROOT_ROOT)/share/data/MNP33_field_map.inp";

  // Target
  int nTargets = 5;
  float shiftTargets[3] = {0.f, 0.f, 0.f}; // Target box global shift, added to posTargetX,Y,Z
  float posTargetX[MaxTargets] = {};
  float posTargetY[MaxTargets] = {};
  float posTargetZ[MaxTargets] = {-2.f, -1.f, 0.f, 1.f, 2.f};
  float thicknessTarget[MaxTargets] = {0.15f, 0.15f, 0.15f, 0.15f, 0.15f};
  float radTarget[MaxTargets] = {0.3f, 0.1f, 0.1f, 0.1f, 0.1f};
  float xsectionTarget[MaxTargets] = {6.6f, 6.6f, 6.6f, 6.6f, 6.6f}; // X-sections in barns
  std::string medTarget[MaxTargets] = {"Lead", "Lead", "Lead", "Lead", "Lead"};

  // VerTel
  int nVerTelPlanes = 5;                  // number of stations
  float shiftVerTel[3] = {0.f, 0.f, 0.f}; // VerTel box global shift, added to posVerTelPlaneX,Y,Z
  float posVerTelPlaneX[MaxVTPlanes] = {};
  float posVerTelPlaneY[MaxVTPlanes] = {};
  float posVerTelPlaneZ[MaxVTPlanes] = {7.1175f, 15.1175f, 20.1175f, 25.1175f, 38.1175f};

  // Muon Stations
  int nMSPlanes = 6;                  // number of stations
  float shiftMS[3] = {0.f, 0.f, 0.f}; // MS global shift, added to posMSPlaneX,Y,Z
  float posMSPlaneX[MaxMSPlanes] = {};
  float posMSPlaneY[MaxMSPlanes] = {};
  float posMSPlaneZ[MaxMSPlanes] = {300.f, 360.f, 530.f, 590.f, 810.f, 850.f};
  float thicknessMSPlane[MaxMSPlanes] = {0.06f, 0.06f, 0.06f, 0.06f, 0.06f, 0.06f};
  float dimXMSPlane[MaxMSPlanes] = {162.f, 240.f, 240.f, 318.f, 417.f, 458.f}; // size in X of rectangular plane. or diameter if dimYMSPlane is 0
  float dimYMSPlane[MaxMSPlanes] = {162.f, 240.f, 240.f, 318.f, 417.f, 458.f}; // size in Y if >0
  float dimXMSPlaneHole[MaxMSPlanes] = {5.f, 5.f, 5.f, 5.f, 5.f, 5.f};         // size in X of rectangular hole. or diameter if dimYMSPlaneHole is 0
  float dimYMSPlaneHole[MaxMSPlanes] = {5.f, 5.f, 5.f, 5.f, 5.f, 5.f};         // size in Y if >0
  std::string medMSPlane[MaxMSPlanes] = {"Silicon", "Silicon", "Silicon", "Silicon", "Silicon", "Silicon"};

  // Absorbers
  int nAbsorberSlices = 10;       // total number of slices
  int nAbsorberSlicesWall = 3;    // of which nAbsorberSlicesWall are the Muon wall
  float posZStartAbsorber = 45.f; // start of the absorber
  float posMuonWall = 610.f;      // start of the muon wall (if!=0 : slice nAbsorberSlices - nAbsorberSlicesWall)
  float thicknessAbsorber[MaxAbsorberSlices] = {10.f, 31.5f, 33.0f, 30.5, 50.f, 50.f, 30.f, /*wall*/ 60.f, 60.f, 60.f};
  float dimXAbsorber[MaxAbsorberSlices] = {39.5f, 52.f, 80.f, 120.f, 170.f, 200.f, 220.f, /*wall*/ 550.f, 550.f, 550.f};      // size in X of rectangular abs. or diameter if dimYAbsorber is 0
  float dimYAbsorber[MaxAbsorberSlices] = {39.5f, 52.f, 80.f, 120.f, 170.f, 200.f, 220.f, /*wall*/ 550.f, 550.f, 550.f};      // size in Y if >0
  float radPlug[MaxAbsorberSlices] = {2.47f, 2.47f, 4.0f, 5.48f, 11.f, 11.f, 11.f, 0.f, 0.f, 0.f};                            // radius of the plug
  float killParticlesBelowEAbsorber[MaxAbsorberSlices] = {0.1, 0.1, 0.1, 0.1, 0.05, 0.01, 0.005, /*wall*/ 0.01, 0.01, 0.001}; // kill entering particles with momentum below this (for speed up)
  std::string medAbsorber[MaxAbsorberSlices] = {"BeO", "BeO", "BeO", "BeO", "Graphite", "Graphite", "Graphite", /*wall*/ "Graphite", "Graphite", "Graphite"};
  std::string medAbsorberPlug[MaxAbsorberSlices] = {"Tungsten", "Tungsten", "Tungsten", "Tungsten", "Tungsten", "Tungsten", ""};

  std::string materialsCutsFile = "$(NA6PROOT_ROOT)/share/data/simcuts.dat";

  NA6PParamDef(NA6PLayoutParam, "layout");
};

#endif
