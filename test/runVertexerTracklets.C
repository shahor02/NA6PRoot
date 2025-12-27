#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVectorD.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include "NA6PVerTelCluster.h"
#include "NA6PVertex.h"
#include "NA6PVerTelReconstruction.h"
#endif

void runVertexerTracklets(int firstEv = 0,
                          int lastEv = 999999,
                          const char* dirSimu = "../testN6Proot/bkgNA49")
//   const char *dirSimu = "Angantyr")
{
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  int nEv = mcTree->GetEntries();
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fc=new TFile(Form("%s/ClustersVerTel.root",dirSimu));
  printf("Open cluster file: %s\n",fc->GetName());
  TTree* tc=(TTree*)fc->Get("clustersVerTel");
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
  tc->SetBranchAddress("VerTel", &vtClusPtr);

  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  TH1F* hdz = new TH1F("hdz", "; z_{rec} - z_{gen} (cm); counts", 100, -0.5, 0.5);
  TH1F* hncontr = new TH1F("hncontr", "; N_{contributors}; counts", 100, -0.5, 999.5);
  TH1F* hnvert = new TH1F("hnvert", "; N_{vertices}; counts", 11, -0.5, 10.5);

  NA6PVerTelReconstruction* vtrec = new NA6PVerTelReconstruction();
  vtrec->setRecoParamFile("test2.ini");
  vtrec->initVertexer();
  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    tc->GetEvent(jEv);
    int nPart = mcArr->size();
    float zVertGen = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zVertGen = curPart.Vz();
        break;
      }
    }
    vtrec->setClusters(vtClus);
    vtrec->runVertexerTracklets();
    std::vector<NA6PVertex> zVertices = vtrec->getVertices();
    int nVertices = zVertices.size();
    hnvert->Fill(nVertices);
    int jv = 0;
    for (auto vert : zVertices) {
      float zRec = vert.getZ();
      if (jv == 0) {
        hdz->Fill(zRec - zVertGen);
        hncontr->Fill(vert.getNContributors());
      }
      printf("Vertex %d, z = %f contrib = %d\n", jv++, zRec, vert.getNContributors());
    }
  }
  vtrec->closeVerticesOutput();

  TCanvas* coutp = new TCanvas("coutp", "", 1200, 600);
  coutp->Divide(3, 1);
  coutp->cd(1);
  hnvert->SetLineWidth(2);
  hnvert->Draw();
  coutp->cd(3);
  hdz->SetLineWidth(2);
  hdz->Draw();
  coutp->cd(2);
  hncontr->SetLineWidth(2);
  hncontr->Draw();
}
