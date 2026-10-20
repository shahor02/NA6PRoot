// NA6PCCopyright

#ifndef NA6P_GENBEAM_H
#define NA6P_GENBEAM_H

#include "NA6PGenerator.h"

class NA6PGenBeam : public NA6PGenerator
{
 public:
  void init() final;
  void generate() final;
  void generatePrimaryVertex(int maxTrials = 1000000) final;

 private:
  float mZPut = -999;

  ClassDefOverride(NA6PGenBeam, 1); // class to shoot beam particle on the targers
};

#endif
