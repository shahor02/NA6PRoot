#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include "MagneticField.h"
#include "NA6PVertex.h"
#include "NA6PVertexerTracks.h"
#endif

void runVertexerTracks(const char* dirSimu = ".")
{
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  TFile* ft = new TFile(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TH1F* hncontr = new TH1F("hncontr", ";N_{contributors}", 100, -0.5, 999.5);
  TH1F* hnvert = new TH1F("hnvert", ";Number of reco vertices", 11, -0.5, 10.5);
  TH2F* hxrecgen = new TH2F("hxrecgen", ";x_{gen} (cm);x_{rec} (cm)", 50, -0.05, 0.05, 50, -0.05, 0.05);
  TH2F* hyrecgen = new TH2F("hyrecgen", ";y_{gen} (cm);y_{rec} (cm)", 50, -0.05, 0.05, 50, -0.05, 0.05);
  TH2F* hzrecgen = new TH2F("hzrecgen", ";z_{gen} (cm);z_{rec} (cm)", 50, -4., 2., 50, -4., 2.);
  TH1F* hdx = new TH1F("hdx", ";x_{rec} - x_{gen} (cm)", 100, -0.05, 0.05);
  TH1F* hdy = new TH1F("hdy", ";y_{rec} - y_{gen} (cm)", 100, -0.05, 0.05);
  TH1F* hdz = new TH1F("hdz", ";z_{rec} - z_{gen} (cm)", 100, -0.5, 0.5);

  NA6PVertexerTracks* vertxr = new NA6PVertexerTracks();
  vertxr->setVerbosity(true);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    printf("Event %d particles = %d tracks = %d\n", jEv, nPart, nTracks);

    double xVertGen = 0;
    double yVertGen = 0;
    double zVertGen = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xVertGen = curPart.Vx();
        yVertGen = curPart.Vy();
        zVertGen = curPart.Vz();
        break;
      }
    }
    vertxr->setBeamX(xVertGen);
    vertxr->setBeamY(yVertGen);
    vertxr->createTracksPool(*trArr);
    std::vector<NA6PVertex> vertices;
    vertxr->findVertices(vertices);
    int nVertices = vertices.size();
    hnvert->Fill(nVertices);
    int jv = 0;
    for (auto vert : vertices) {
      double xRec = vert.getX();
      double yRec = vert.getY();
      double zRec = vert.getZ();
      if (jv == 0) {
        hxrecgen->Fill(xVertGen, xRec);
        hyrecgen->Fill(yVertGen, yRec);
        hzrecgen->Fill(zVertGen, zRec);
        hdx->Fill(xRec - xVertGen);
        hdy->Fill(yRec - yVertGen);
        hdz->Fill(zRec - zVertGen);
        hncontr->Fill(vert.getNContributors());
      }
      printf("Vertex %d, z = %f contrib = %d\n", jv++, zRec, vert.getNContributors());
    }
  }

  TCanvas* cv = new TCanvas("cv", "", 1000, 500);
  cv->Divide(2, 1);
  cv->cd(1);
  hnvert->Draw();
  cv->cd(2);
  hncontr->Draw();

  TCanvas* coutp = new TCanvas("coutp", "", 1400, 800);
  coutp->Divide(3, 2);
  coutp->cd(1);
  hxrecgen->Draw("colz");
  coutp->cd(2);
  hyrecgen->Draw("colz");
  coutp->cd(3);
  hzrecgen->Draw("colz");
  coutp->cd(4);
  hdx->Draw();
  coutp->cd(5);
  hdy->Draw();
  coutp->cd(6);
  hdz->Draw();
}
