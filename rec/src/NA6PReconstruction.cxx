// NA6PCCopyright

#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TSystem.h>
#include "MagneticField.h"
#include "NA6PReconstruction.h"
#include "Propagator.h"

ClassImp(NA6PReconstruction)

  bool NA6PReconstruction::init(const std::string& filename, const std::string& geoname)
{
  // initialize magnetic field
  if (mIsInitialized) {
    LOGP(info, "Reconstruction already intialized");
    return true;
  }
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    Propagator::loadField();
  }
  // load geometry
  if (!gGeoManager && !Propagator::loadGeometry(filename, geoname)) {
    return false;
  }
  mIsInitialized = true;
  return true;
}

void NA6PReconstruction::createClustersOutput()
{
  LOGP(warning, "createClustersOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::clearClusters()
{
  LOGP(warning, "clearClusters called for {}, this should not happen", getName());
}

void NA6PReconstruction::writeClusters()
{
  LOGP(warning, "writeClusters called for {}, this should not happen", getName());
}

void NA6PReconstruction::closeClustersOutput()
{
  LOGP(warning, "closeClustersOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::createVerticesOutput()
{
  LOGP(warning, "createVerticesOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::clearVertices()
{
  LOGP(warning, "clearVertices called for {}, this should not happen", getName());
}

void NA6PReconstruction::writeVertices()
{
  LOGP(warning, "writeVertices called for {}, this should not happen", getName());
}

void NA6PReconstruction::closeVerticesOutput()
{
  LOGP(warning, "closeVerticesOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::createTracksOutput()
{
  LOGP(warning, "createTracksOutput called for {}, this should not happen", getName());
}

void NA6PReconstruction::clearTracks()
{
  LOGP(warning, "clearTracks called for {}, this should not happen", getName());
}

void NA6PReconstruction::writeTracks()
{
  LOGP(warning, "writeTracks called for {}, this should not happen", getName());
}

void NA6PReconstruction::closeTracksOutput()
{
  LOGP(warning, "closeTracksOutput called for {}, this should not happen", getName());
}
