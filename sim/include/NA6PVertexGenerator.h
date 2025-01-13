// NA6PCCopyright

#ifndef NA6P_VERTEX_GENERATOR_H_
#define NA6P_VERTEX_GENERATOR_H_

#include "NA6PBeamParam.h"

class VertexGenerator
{
 public:
  bool generate(int maxTrials = 100000) const;
};

#endif
