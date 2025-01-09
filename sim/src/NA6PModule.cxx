// NA6PCCopyright
#include "NA6PModule.h"
#include <fairlogger/Logger.h>

void NA6PModule::setActiveID(int i)
{
  if (i < 0 || i >= MaxActiveID) {
    LOGP(fatal, "Module ActiveID must be in the range [0:MaxActiveID), {} requested", i);
  }
  mActiveID = i;
}

int NA6PModule::composeSensorVolID(int id) const
{
  if (mVolIDOffset < 0) {
    LOGP(fatal, "{}: cannot assign sensor ID before module ID is assigned", getName());
  }
  if (id < 0 || id >= MaxVolID - MaxNonSensID) {
    LOGP(fatal, "{}: sensor ID must be in [0:{}) range, {} passed", getName(), MaxVolID - MaxNonSensID, id);
  }
  return mVolIDOffset + id + MaxNonSensID;
}

int NA6PModule::composeNonSensorVolID(int id) const
{
  if (mVolIDOffset < 0) {
    LOGP(fatal, "{}: cannot assign sensor ID before module ID is assigned", getName());
  }
  if (id < 0 || id >= MaxNonSensID) {
    LOGP(fatal, "{}: non-sensor ID must be in [0:{}) range, {} passed", getName(), MaxNonSensID, id);
  }
  return mVolIDOffset + id;
}

bool NA6PModule::stepManager(int vid)
{
  LOGP(warning, "Default stepManager called for {} with volID={}, this should not happen", getName(), vid);
  return false;
}

void NA6PModule::createHitsOutput(const std::string&)
{
  LOGP(warning, "createHitsOutput called for {}, this should not happen", getName());
}

void NA6PModule::closeHitsOutput()
{
  LOGP(warning, "closeHitsOutput called for {}, this should not happen", getName());
}

void NA6PModule::writeHits(const std::vector<int>&)
{
  LOGP(warning, "writeHits called for {}, this should not happen", getName());
}

std::string NA6PModule::addName(const std::string& n)
{
  return fmt::format("{}_{}", n, getName());
}
