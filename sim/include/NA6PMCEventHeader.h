// NA6PCCopyright

#ifndef NA6P_MCEVENT_HEADER_H
#define NA6P_MCEVENT_HEADER_H

#include <Rtypes.h>
#include <vector>
#include "NA6PMCGenHeader.h"

class NA6PMCEventHeader
{
 public:

  uint32_t getEventID() const { return mEventID; }
  uint32_t getRunNumber() const { return mRunNumber; }
  float getVX() const { return mVX; }
  float getVY() const { return mVY; }
  float getVZ() const { return mVZ; }
  uint32_t getNTracks() const { return mNTracks; }
  uint32_t getNPrimaries() const { return mNPrimaries; }
  const std::vector<NA6PMCGenHeader>& getGenHeaders() const { return mGenHeaders; }
  std::vector<NA6PMCGenHeader>& getGenHeaders() { return mGenHeaders; }
  size_t getNGenHeaders() const { return mGenHeaders.size(); }
  const NA6PMCGenHeader& getGenHeader(int i) const { return mGenHeaders[i]; }
  NA6PMCGenHeader& getGenHeader(int i) { return mGenHeaders[i]; }
  
  void setEventID(uint32_t eventID) { mEventID = eventID; }
  void setRunNumber(uint32_t runNumber) { mRunNumber = runNumber; }
  void setVX(float vx) { mVX = vx; }
  void setVY(float vy) { mVY = vy; }
  void setVZ(float vz) { mVZ = vz; }
  void setNTracks(uint32_t nTrack) { mNTracks = nTrack; }
  void setNPrimaries(uint32_t nPrimaries) { mNPrimaries = nPrimaries; }
  void setGenHeaders(const std::vector<NA6PMCGenHeader>& genHeaders) { mGenHeaders = genHeaders; }
  void addGenHeader(const NA6PMCGenHeader& h) { mGenHeaders.push_back(h); }

  void clear();
  void print(bool genh = true) const;
  std::string asString(bool genh = true) const;
  
 protected:
  uint32_t mEventID = 0;
  uint32_t mRunNumber = 0;
  float mVX = 0.f;   // Primary vertex x
  float mVY = 0.f;   // Primary vertex y
  float mVZ = 0.f;   // Primary vertex z

  uint32_t mNTracks = 0;      // total number of tracks
  uint32_t mNPrimaries = 0;   // total number of primaries from generators

  std::vector<NA6PMCGenHeader> mGenHeaders; // generator headers
  
  ClassDefNV(NA6PMCEventHeader,1)
};

#endif
