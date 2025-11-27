// NA6PCCopyright

#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TSystem.h>
#include "MagneticField.h"
#include "NA6PReconstruction.h"

ClassImp(NA6PReconstruction)


bool NA6PReconstruction::init(const char* filename, const char* geoname){
  // initialize magnetic field
  if (mIsInitialized) {
    LOGP(info,"Reconstruction alreary intialized");
    return true;
  }
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    auto magField = new MagneticField();
    magField->loadField();
    magField->setAsGlobalField();
  }
  // load geometry
  if (gGeoManager){
    LOGP(info,"Geometry was already loaded");
    mIsInitialized = true;
    return true;
  }
  if (gSystem->Exec(Form("ls -l %s > /dev/null",filename)) != 0){
    LOGP(error,"filename {} does not exist",filename);
    return false;
  }
  TFile* f = TFile::Open(filename);
  gGeoManager = (TGeoManager*)f->Get(geoname);
  if (gGeoManager){
    mIsInitialized = true;
    return true;
  }else{
    LOGP(error,"No geometry with name {} found in file {}",geoname,filename);
    return false;
  }
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

