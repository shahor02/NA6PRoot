// NA6PCCopyright

#ifndef NA6P_GENBOX_H
#define NA6P_GENBOX_H

#include "NA6PGenerator.h"

class NA6PGenBox : public NA6PGenerator
{
 public:
  NA6PGenBox(int pdg=211, int npart=1, const std::string& name = "genbox") : NA6PGenerator(name), mPDGCode(pdg), mNTracks(npart) {}
  ~NA6PGenBox() override = default;
  void generate() override;

  void init() override;
  void setPDGCode(int c) { mPDGCode = c; }
  auto getPDGCode() const { return mPDGCode; }
  void setNTracks(int n) { mNTracks = n; }
  auto getNTracks() const { return mNTracks; }

protected:

  int mPDGCode = 211;     // particle type
  int mNTracks = 1;       // number of particles to generate
  
  ClassDefOverride(NA6PGenBox,1);  // Square box (in variables of NA6PGenCutParam) random generator  
};

#endif
