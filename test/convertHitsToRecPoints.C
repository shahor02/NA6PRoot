#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TRandom3.h>
#include "NA6PVerTelHit.h"
#include "NA6PBaseCluster.h"
#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include "NA6PVerTelReconstruction.h"
#endif

void ConvertHitsToRecPoints(){
  TFile* fh=new TFile("HitsVerTel.root");
  TTree* th=(TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);
  int nEv=th->GetEntriesFast();

  NA6PVerTelReconstruction* recVT = new NA6PVerTelReconstruction();
  recVT->createClustersOutput();
  for(int jEv=0; jEv<nEv; jEv++){
    th->GetEvent(jEv);
    int nHits=vtHits.size();
    printf("Event %d nHits=%d\n",jEv,nHits);
    recVT->clearClusters();
    recVT->hitsToRecPoints(vtHits);
    recVT->writeClusters();
  }
  recVT->closeClustersOutput();
}
