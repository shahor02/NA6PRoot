// NA6PCCopyright

#ifndef NA6P_GENHEPMC_H
#define NA6P_GENHEPMC_H

#include "NA6PGenerator.h"

class NA6PGenHepMC : public NA6PGenerator
{
 public:
  using NA6PGenerator::NA6PGenerator;
  
  NA6PGenHepMC() = default;
  NA6PGenHepMC(const std::string& name, const std::string& filename);
  ~NA6PGenHepMC() override = default;
  void generate() override;
  
  void init() override;

  void SetFileName(const std::string& name){mFileName = name;}
  
 protected:
  std::string mFileName = {};
  ulong       mReadEvents = 0;
  
  ClassDefOverride(NA6PGenHepMC, 1); // 
};

#endif
