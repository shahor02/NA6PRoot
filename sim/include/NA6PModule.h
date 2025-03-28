// NA6PCCopyright
#ifndef NA6P_MODULE_H_
#define NA6P_MODULE_H_

#include <string>
#include "TLorentzVector.h" // for TLorentzVector
#include "NA6PSimMisc.h"    // for TLorentzVector

class TGeoVolume;

// Base (pasive or active) module of NA6P

class NA6PModule
{
 public:
  static constexpr int MaxActiveID = 9;    // max number of active modules
  static constexpr int MaxNonSensID = 100; // non-sensors can be assigned mVolIDOffset <= volID < mVolIDOffset + MaxNonSensID
  static constexpr int MaxVolID = 1000;    // sensor can be assigne mVolIDOffset + MaxNonSensID <= volID < MaxVolID

  NA6PModule(std::string name) : mName(name) {}
  virtual ~NA6PModule() = default;
  virtual void createMaterials() = 0;
  virtual void createGeometry(TGeoVolume* base) = 0;
  virtual bool stepManager(int); // called by NA6PMC Stepping, passing the volume copyID (from which sensID can be extracted using volID2SensID
  virtual size_t getNHits() const { return 0; }
  virtual void setAlignableEntries() {}

  auto isActive() const { return mActiveID >= 0; }
  auto getActiveID() const { return mActiveID; }

  void setID(int v) { mVolIDOffset = (mID = v) * MaxVolID; }
  auto getID() const { return mID; }
  auto getVolIDOffset() const { return mVolIDOffset; }

  const std::string& getName() const { return mName; }

  static bool isSensor(int v) { return v % MaxVolID >= MaxNonSensID; }
  static int volID2SensID(int v) { return v % MaxVolID - MaxNonSensID; } // if < 0 : non-sensor
  static int volID2NonSensID(int v) { return v % MaxVolID; }
  static int volID2ModuleID(int v) { return v / MaxVolID - 1; }
  int composeNonSensorVolID(int id) const;
  int composeSensorVolID(int id) const;

  void setVerbosity(int v) { mVerbosity = v; }
  auto getVerbosity() const { return mVerbosity; }

  virtual void clearHits() {}

  static int getActiveIDBit(int id) { return 0x1 << (id + 14); } // convert ActiveID to TObject user bit
  static int testActiveIDBits(const TObject& obj) { return obj.TestBits(((0x1 << MaxActiveID) - 1) << 14) >> 14; }
  static int testActiveIDOrKeepBits(const TObject& obj) { return obj.TestBits((((0x1 << MaxActiveID) - 1) << 14) | UserHook::KeepParticleBit) >> 14; }

  bool testActiveIDBit(const TObject& obj) const { return mActiveID >= 0 ? (testActiveIDBits(obj) & (0x1 << mActiveID)) != 0 : false; }
  int getActiveIDBit() const { return mActiveID >= 0 ? getActiveIDBit(mActiveID) : 0; }

  virtual void createHitsOutput(const std::string&);
  virtual void closeHitsOutput();
  virtual void writeHits(const std::vector<int>&);

 protected:
  void setActiveID(int i);
  std::string addName(const std::string& n);

  int mActiveID = -1; // uinuqeID of the active detector
  int mID = -1;       // module ID
  int mVolIDOffset = -1;
  ; // sensitive volumes IDs of this module will have this offsets
  int mVerbosity = 0;
  std::string mName{};

  // transient track data for hit creation
  struct TrackData {                   // this is transient
    bool mHitStarted = false;          //! hit creation started
    unsigned char mTrkStatusStart = 0; //! track status flag
    TLorentzVector mPositionStart{};   //! position at entrance
    TLorentzVector mMomentumStart{};   //! momentum
    double mEnergyLoss;                //! energy loss
  } mTrackData;                        //!
};

#endif
