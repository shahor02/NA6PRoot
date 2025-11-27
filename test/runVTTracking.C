#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PBaseCluster.h"
#include "NA6PVerTelReconstruction.h"
#endif

void runVTTracking(int firstEv = 0,
		   int lastEv = 99999,
		   const char *dirSimu = "Angantyr")
{
  
  // kine file needed to get the primary vertex position (temporary)
  TFile* fk=new TFile(Form("%s/MCKine.root",dirSimu));
  TTree* mcTree=(TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fc=new TFile(Form("%s/ClustersVerTel.root",dirSimu));
  printf("Open cluster file: %s\n",fc->GetName());
  TTree* tc=(TTree*)fc->Get("clustersVerTel");
  std::vector<NA6PBaseCluster> vtClus, *vtClusPtr = &vtClus;
  tc->SetBranchAddress("VerTel", &vtClusPtr);

  NA6PVerTelReconstruction* vtrec = new NA6PVerTelReconstruction();
  vtrec->init(Form("%s/geometry.root",dirSimu));

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
    vtrec->setPrimaryVertexPosition(0.,0.,zvert);
    vtrec->setClusters(vtClus);
    vtrec->runTracking();
  }
  vtrec->closeTracksOutput();
  
  timer.Stop();
  timer.Print();
}
    
