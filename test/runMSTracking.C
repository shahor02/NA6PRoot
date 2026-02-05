#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PMuonSpecCluster.h"
#include "NA6PMuonSpecReconstruction.h"
#endif

void runMSTracking(int firstEv = 0,
		   int lastEv = 99999,
		   const char *dirSimu = "../tst")
{
  
  // kine file needed to get the primary vertex position (temporary)
  TFile* fk=new TFile(Form("%s/MCKine.root",dirSimu));
  TTree* mcTree=(TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fc=new TFile(Form("%s/ClustersMuonSpec.root",dirSimu));
  printf("Open cluster file: %s\n",fc->GetName());
  TTree* tc=(TTree*)fc->Get("clustersMuonSpec");
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;
  tc->SetBranchAddress("MuonSpec", &msClusPtr);
  NA6PMuonSpecReconstruction* msrec = new NA6PMuonSpecReconstruction();
  msrec->init(Form("%s/geometry.root",dirSimu));
  auto msTracker = msrec->getTracker();
  // Example of configuring the tracker for Muon Spectrometer
  //msTracker->setNLayers(6);
  //msTracker->setNumberOfIterations(2);
  //msTracker->setIterationParams(0,0.06,0.1,6.,0.6,0.05,0.05,5.,5.,5.,6);
  int nEv=tc->GetEntries();
  if(lastEv>nEv || lastEv<0) lastEv=nEv;
  if(firstEv<0) firstEv=0;

  TStopwatch timer;
  timer.Start();

  for(int jEv=firstEv; jEv<lastEv; jEv++){
    mcTree->GetEvent(jEv);
    int nPart=mcArr->size();
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for(int jp=0; jp<nPart; jp++){
      auto curPart=mcArr->at(jp);
      if(curPart.IsPrimary()){
        zvert =  curPart.Vz();
        break;
      }
    }
    tc->GetEvent(jEv);
    msrec->setPrimaryVertexPosition(0.,0.,zvert);
    msrec->setClusters(msClus);
    msrec->runTracking();
  }
  msrec->closeTracksOutput();
  
  timer.Stop();
  timer.Print();
}
    
