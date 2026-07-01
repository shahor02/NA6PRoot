#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PMuonSpecModularHit.h"
#include "Propagator.h"
#include "NA6PMatch.h"
#include "NA6PTrack.h"
#include "NA6PTreeStreamRedirector.h"
#include "HistoManager.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#endif

//----- CUTS
float cutEtaMin = 1.5f;
float cutEtaMax = 4.5f;
float cutPtMin = 0.1f;

int minLrVT = 5;
int minLrMS = 6;
//-----

struct MCTrackInfo {
  enum TrStatus : uint8_t { New,
                            Discard,
                            Accept };
  TrStatus status = New;
  int8_t charge = 0;
  int8_t nvtLr = 0;      // number of layers with clusters
  int8_t nmsLr = 0;      // number of layers with clusters
  uint8_t nVTTracks = 0; // number of tracks bearing the label of this MC track
  uint8_t nMSTracks = 0; // number of tracks bearing the label of this MC track
  uint8_t nMatches = 0;  // number of tracks bearing the label of this MC track
  int bestVTTrack = -1;
  int bestMSTrack = -1;
  int bestMatch = -1;

  std::vector<int> vtClID;
  std::vector<int> msClID;
  std::array<float, 3> momMS{}, posMS{};
  bool isReconstructableVT() const { return status == TrStatus::Accept && nvtLr >= minLrVT; }
  bool isReconstructableMS() const { return status == TrStatus::Accept && nmsLr >= minLrMS; }
};

const char* prefHisto[3] = {"VT", "MS", "Match"};
const char* prefAxis[2] = {"1/p", "#eta"};
const char* prefAxisS[2] = {"pInv", "eta"};
const char* parName[5] = {"X", "Y", "Tx", "Ty", "q/Pxz"};
const float parMaxDiff[5] = {0.02, 0.02, 0.005, 0.005, 0.1};

const std::array<int, 3> HMOffsets = {0, 1000, 2000};
const int nBinsPtI = 30;
const int nBinsEta = 30;
const int nBinsDiffKin = 50;

na6p::HistoManager* hm = nullptr;

std::unordered_map<int, MCTrackInfo> mcTruth;
std::unique_ptr<NA6PTreeStreamRedirector> outPerf;

TTree *mcTree = nullptr, *vtClTree = nullptr, *vtTrTree = nullptr, *msClTree = nullptr, *msTrTree = nullptr, *matchTree = nullptr;
TTree* msHitsTree = nullptr;
std::vector<TParticle>* mcArr = nullptr;
std::vector<NA6PVerTelCluster>* vtClVec = nullptr;
std::vector<NA6PTrack>* vtTrVec = nullptr;
std::vector<NA6PMuonSpecCluster>* msClVec = nullptr;
std::vector<NA6PTrack>* msTrVec = nullptr;
std::vector<NA6PMatch>* matchVec = nullptr;
std::vector<NA6PMuonSpecModularHit>* msHitsVec = nullptr;
Long64_t nEntries = 0;

void setupInputs(bool checkVTTracks, bool checkMSTracks, bool checkMatches, const char* dirSimu, const char* dirRec);
void loadEventData(int i);
void buildReconstuctableMCInfo();
void checkMSTracking();
void checkMatching();
void checkVTTracking();
void createOutput(bool checkVTTracks, bool checkMSTracks, bool checkMatches);
void bookHistos(bool checkVTTracks, bool checkMSTracks, bool checkMatches);
void finalizeHistos();
void drawHistos(const na6p::HistoManager* hman, bool same = false);

void runPerfStudy(int firstEv = 0,
                  int lastEv = 99999,
                  bool checkVTTracks = true,
                  bool checkMSTracks = true,
                  bool checkMatches = true,
                  const char* dirSimu = ".",
                  const char* dirRec = ".")
{
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }
  setupInputs(checkVTTracks, checkMSTracks, checkMatches, dirSimu, dirRec);
  bookHistos(checkVTTracks, checkMSTracks, checkMatches);
  if (lastEv < 1) {
    lastEv = nEntries;
  } else if (lastEv > (int)nEntries) {
    lastEv = nEntries;
  }
  outPerf = std::make_unique<NA6PTreeStreamRedirector>("recperf.root");

  for (int ev = firstEv; ev < lastEv; ev++) {
    loadEventData(ev);
    buildReconstuctableMCInfo();
    if (checkVTTracks) {
      checkVTTracking();
    }
    if (checkMSTracks) {
      checkMSTracking();
    }
    if (checkMatches) {
      checkMatching();
    }

    createOutput(checkVTTracks, checkMSTracks, checkMatches);
  }

  outPerf->Close();
  outPerf.reset();
  finalizeHistos();
  drawHistos(hm);
}

void bookHistos(bool checkVTTracks, bool checkMSTracks, bool checkMatches)
{
  hm = new na6p::HistoManager();
  bool doHisto[3] = {checkVTTracks, checkMSTracks, checkMatches};
  float cutAxisMin[3][2] = {{0.f, cutEtaMin}, {0.f, cutEtaMin}, {0.f, cutEtaMin}};   // min 1/p and eta for VT, MS and Matches
  float cutAxisMax[3][2] = {{1.f, cutEtaMax}, {0.4f, cutEtaMax}, {0.4f, cutEtaMax}}; // max 1/p and eta for VT, MS and Matches
  int nbinsAxis[2] = {nBinsPtI, nBinsEta};
  for (int itp = 0; itp < 3; itp++) {
    if (!doHisto[itp]) {
      continue;
    }
    int offset = HMOffsets[itp];
    for (int iax = 0; iax < 2; iax++) {
      // efficiency and purity
      auto h1mc = new TH1F(Form("reconstructable_%s_%s", prefHisto[itp], prefAxis[iax]), Form("reconstructable %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1rec = new TH1F(Form("reconstructed_%s_%s", prefHisto[itp], prefAxis[iax]), Form("reconstructed %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1recCorr = new TH1F(Form("correct_%s_%s", prefHisto[itp], prefAxis[iax]), Form("reconstructed correctly %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1recEff = new TH1F(Form("receff_%s_%s", prefHisto[itp], prefAxis[iax]), Form("rec.eff %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1recPur = new TH1F(Form("purity_%s_%s", prefHisto[itp], prefAxis[iax]), Form("purity, %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      hm->addHisto(h1mc, offset + iax * 100 + 0);
      hm->addHisto(h1rec, offset + iax * 100 + 1);
      hm->addHisto(h1recCorr, offset + iax * 100 + 2);
      hm->addHisto(h1recEff, offset + iax * 100 + 3);
      hm->addHisto(h1recPur, offset + iax * 100 + 4);

      //
      // resolutions
      for (int ip = 0; ip < 5; ip++) {
        auto h2 = new TH2F(Form("difRecMC_%s_vs_%s_%s", parName[ip], prefAxis[iax], prefHisto[itp]),
                           Form("%s Rec-MC %s vs %s;%s;%s data-MC", prefHisto[itp], parName[ip], prefAxis[iax], prefAxis[iax], parName[ip]),
                           nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax], nBinsDiffKin, -parMaxDiff[ip], parMaxDiff[ip]);
        hm->addHisto(h2, offset + iax * 100 + 500 + ip * 10);
      }
    }
  }
  hm->sumw2();
}

void finalizeHistos()
{
  std::unique_ptr<TF1> gs = std::make_unique<TF1>("gs", "gaus", -10., 10.);
  TObjArray harr;
  harr.SetOwner(true);
  for (int itp = 0; itp < 3; itp++) {
    int offset = HMOffsets[itp];
    if (!hm->getHisto(offset + 0 * 100 + 0)) {
      continue;
    }
    for (int iax = 0; iax < 2; iax++) {
      hm->getHisto(offset + iax * 100 + 3)->Divide(hm->getHisto(offset + iax * 100 + 1), hm->getHisto(offset + iax * 100 + 0), 1., 1., "B");
      hm->getHisto(offset + iax * 100 + 4)->Divide(hm->getHisto(offset + iax * 100 + 2), hm->getHisto(offset + iax * 100 + 1), 1., 1., "B");
      //
      // resolutions
      for (int ip = 0; ip < 5; ip++) {
        int hid = offset + iax * 100 + 500 + ip * 10;
        auto h2 = hm->getHisto2F(hid);
        h2->FitSlicesY(gs.get(), 0, -1, 0, "QNR", &harr);
        harr.SetOwner(true);
        TH1* hmean = (TH1*)harr.RemoveAt(1);
        if (hmean) {
          hmean->SetTitle(Form("<%s>", h2->GetTitle()));
          hm->addHisto(hmean, hid + 1);
        }
        TH1* hsig = (TH1*)harr.RemoveAt(2);
        if (hsig) {
          hsig->SetTitle(Form("#sigma(%s)", h2->GetTitle()));
          hm->addHisto(hsig, hid + 2);
        }
        harr.Delete();
      }
    }
  }
}

TCanvas* cnvEff[3] = {};
TCanvas* cnvKin[3][2] = {{0, 0}, {0, 0}, {0, 0}};

void drawHistos(const na6p::HistoManager* hman, bool same)
{
  gStyle->SetTitleW(0.5);
  gStyle->SetOptStat(0);
  const char* opt = same ? "same" : "";
  for (int itp = 0; itp < 3; itp++) {
    int offset = HMOffsets[itp];
    if (!hman->getHisto(offset + 0 * 100 + 0)) {
      continue;
    }
    //--------------- Eff
    if (!cnvEff[itp]) {
      cnvEff[itp] = new TCanvas(Form("cnv%s", prefHisto[itp]), Form("cnv%s", prefHisto[itp]), 600, 800);
      cnvEff[itp]->Divide(2, 2);
    }
    for (int iax = 0; iax < 2; iax++) {
      auto heff = hman->getHisto(offset + iax * 100 + 3);
      auto hpur = hman->getHisto(offset + iax * 100 + 4);
      cnvEff[itp]->cd(1 + 2 * iax);
      heff->Draw(opt);
      gPad->SetGrid();
      cnvEff[itp]->cd(2 + 2 * iax);
      hpur->Draw(opt);
      gPad->SetGrid();
    }
    //--------------- Res
    for (int iax = 0; iax < 2; iax++) {
      if (!cnvKin[itp][iax]) {
        cnvKin[itp][iax] = new TCanvas(Form("cnv%s%s", prefHisto[itp], prefAxisS[iax]), Form("cnv%s%s", prefHisto[itp], prefAxisS[iax]), 600, 800);
        cnvKin[itp][iax]->Divide(2, 5);
        for (int ip = 0; ip < 5; ip++) {
          int hid = offset + iax * 100 + 500 + ip * 10;
          cnvKin[itp][iax]->cd(1 + ip * 2);
          gPad->SetBottomMargin(0.15);
          hman->getHisto(hid + 1)->Draw(opt);
          hman->getHisto(hid + 1)->GetXaxis()->SetLabelSize(0.04);
          hman->getHisto(hid + 1)->GetYaxis()->SetLabelSize(0.04);
          hman->getHisto(hid + 1)->GetXaxis()->SetTitleSize(0.05);
          hman->getHisto(hid + 1)->GetYaxis()->SetTitleSize(0.05);
          gPad->SetGrid();
          cnvKin[itp][iax]->cd(2 + ip * 2);
          gPad->SetBottomMargin(0.15);
          hman->getHisto(hid + 2)->Draw(opt);
          hman->getHisto(hid + 2)->GetXaxis()->SetLabelSize(0.04);
          hman->getHisto(hid + 2)->GetYaxis()->SetLabelSize(0.04);
          hman->getHisto(hid + 2)->GetXaxis()->SetTitleSize(0.05);
          hman->getHisto(hid + 2)->GetYaxis()->SetTitleSize(0.05);
          gPad->SetGrid();
        }
      }
    }
  }
}

void createOutput(bool checkVTTracks, bool checkMSTracks, bool checkMatches)
{
  NA6PTrackParCov dummy;
  dummy.invalidate();
  auto prop = Propagator::Instance();
  for (auto& mcInfoPair : mcTruth) {
    const auto& mcTrackInfo = mcInfoPair.second;
    if (mcTrackInfo.status != MCTrackInfo::Accept || !(mcTrackInfo.isReconstructableVT() || mcTrackInfo.isReconstructableMS()))
      continue;

    const auto& partMC = (*mcArr)[mcInfoPair.first];
    // convert to TrackParam
    float mom[3] = {(float)partMC.Px(), (float)partMC.Py(), (float)partMC.Pz()}, pos[3] = {(float)partMC.Vx(), (float)partMC.Vy(), (float)partMC.Vz()};
    NA6PTrackPar trMC(pos, mom, mcTrackInfo.charge);
    float pInvMC = 1.f / trMC.getP(), etaMC = trMC.getEta();

    if (checkVTTracks && mcTrackInfo.isReconstructableVT()) { // output generated info at the vertex
      hm->getHisto(HMOffsets[0] + 0 * 100 + 0)->Fill(pInvMC); // reconstructable
      hm->getHisto(HMOffsets[0] + 1 * 100 + 0)->Fill(etaMC);
      NA6PTrackParCov trCopy;
      int nVTHits = 0, nClones = 0, partID = -1;
      float chi2 = -1.f;
      if (mcTrackInfo.bestVTTrack >= 0) {
        const auto& vtTr = (*vtTrVec)[mcTrackInfo.bestVTTrack];
        trCopy = vtTr;
        nVTHits = vtTr.getNHits();
        nClones = mcTrackInfo.nVTTracks;
        partID = vtTr.getParticleID();
        chi2 = vtTr.getChi2();
        if (!prop->propagateToZ(trCopy, pos[2])) {
          trCopy.invalidate();
        }
        hm->getHisto(HMOffsets[0] + 0 * 100 + 1)->Fill(pInvMC); // reconstructed
        hm->getHisto(HMOffsets[0] + 1 * 100 + 1)->Fill(etaMC);
        //
        if (partID >= 0) {
          hm->getHisto(HMOffsets[0] + 0 * 100 + 2)->Fill(pInvMC); // reconstructed correctly
          hm->getHisto(HMOffsets[0] + 1 * 100 + 2)->Fill(etaMC);
        }
        if (trCopy.isValid()) { // resolutions
          for (int ip = 0; ip < 5; ip++) {
            hm->getHisto2F(HMOffsets[0] + 0 * 100 + 500 + ip * 10)->Fill(pInvMC, trCopy.getParam(ip) - trMC.getParam(ip)); // reconstructed
            hm->getHisto2F(HMOffsets[0] + 1 * 100 + 500 + ip * 10)->Fill(etaMC, trCopy.getParam(ip) - trMC.getParam(ip));  // reconstructed
          }
        }
      } else {
        trCopy = dummy;
      }
      (*outPerf) << "perfVT" << "evID=" << int(mcTree->GetReadEntry()) << "mcTrID=" << mcInfoPair.first << "pdg=" << partMC.GetPdgCode() << "mcTr=" << trMC
                 << "nvtLr=" << mcTrackInfo.nvtLr << "nmsLr=" << mcTrackInfo.nmsLr
                 << "vtTr=" << trCopy << "nVTCl=" << nVTHits << "chi2vt=" << chi2 << "nClones=" << nClones << "mcID=" << partID << "\n";
    }
    //-----------------------

    if (checkMSTracks && mcTrackInfo.isReconstructableMS()) { // output generated info at the vertex
      NA6PTrackParCov trCopy;
      NA6PTrackPar msMCTr(mcTrackInfo.posMS, mcTrackInfo.momMS, mcTrackInfo.charge);
      hm->getHisto(HMOffsets[1] + 0 * 100 + 0)->Fill(pInvMC); // reconstructable
      hm->getHisto(HMOffsets[1] + 1 * 100 + 0)->Fill(etaMC);

      int nMSHits = 0, nClones = 0, partID = -1;
      float chi2 = -1.f;
      if (mcTrackInfo.bestMSTrack >= 0) {
        const auto& msTr = (*msTrVec)[mcTrackInfo.bestMSTrack];
        trCopy = msTr;
        nMSHits = msTr.getNHits();
        nClones = mcTrackInfo.nMSTracks;
        partID = msTr.getParticleID();
        chi2 = msTr.getChi2();
        if (!prop->propagateToZ(trCopy, msMCTr.getZ())) {
          trCopy.invalidate();
        }
        hm->getHisto(HMOffsets[1] + 0 * 100 + 1)->Fill(pInvMC); // reconstructed
        hm->getHisto(HMOffsets[1] + 1 * 100 + 1)->Fill(etaMC);
        //
        if (partID >= 0) {
          hm->getHisto(HMOffsets[1] + 0 * 100 + 2)->Fill(pInvMC); // reconstructed correctly
          hm->getHisto(HMOffsets[1] + 1 * 100 + 2)->Fill(etaMC);
        }
        if (trCopy.isValid()) { // resolutions
          for (int ip = 0; ip < 5; ip++) {
            hm->getHisto2F(HMOffsets[1] + 0 * 100 + 500 + ip * 10)->Fill(pInvMC, trCopy.getParam(ip) - msMCTr.getParam(ip)); // reconstructed
            hm->getHisto2F(HMOffsets[1] + 1 * 100 + 500 + ip * 10)->Fill(etaMC, trCopy.getParam(ip) - msMCTr.getParam(ip));  // reconstructed
          }
        }
      } else {
        trCopy = dummy;
      }
      (*outPerf) << "perfMS" << "evID=" << int(mcTree->GetReadEntry()) << "mcTrID=" << mcInfoPair.first << "pdg=" << partMC.GetPdgCode()
                 << "mcTrMS0=" << msMCTr << "nvtLr=" << mcTrackInfo.nvtLr << "nmsLr=" << mcTrackInfo.nmsLr
                 << "msTr=" << trCopy << "nMSCl=" << nMSHits << "chi2ms=" << chi2 << "nClones=" << nClones << "mcID=" << partID << "\n";
    }
    //-----------------------
    if (checkMatches && mcTrackInfo.isReconstructableVT() && mcTrackInfo.isReconstructableMS()) { // output generated info at the vertex
      NA6PTrackParCov trCopy;
      int nClones = 0, partID = -1;
      float chi2refit = -1.f, chi2match = -1.f;
      hm->getHisto(HMOffsets[2] + 0 * 100 + 0)->Fill(pInvMC); // matchable
      hm->getHisto(HMOffsets[2] + 1 * 100 + 0)->Fill(etaMC);
      if (mcTrackInfo.bestMatch >= 0) {
        const auto& mtcTr = (*matchVec)[mcTrackInfo.bestMatch];
        trCopy = mtcTr;
        nClones = mcTrackInfo.nMatches;
        partID = mtcTr.getParticleID();
        chi2refit = mtcTr.getChi2Refit();
        chi2match = mtcTr.getChi2Match();
        if (!prop->propagateToZ(trCopy, pos[2])) {
          trCopy.invalidate();
        }
        hm->getHisto(HMOffsets[2] + 0 * 100 + 1)->Fill(pInvMC); // matched
        hm->getHisto(HMOffsets[2] + 1 * 100 + 1)->Fill(etaMC);
        //
        if (partID >= 0) {
          hm->getHisto(HMOffsets[2] + 0 * 100 + 2)->Fill(pInvMC); // matched correctly
          hm->getHisto(HMOffsets[2] + 1 * 100 + 2)->Fill(etaMC);
        }
        if (trCopy.isValid()) { // resolutions
          for (int ip = 0; ip < 5; ip++) {
            hm->getHisto2F(HMOffsets[2] + 0 * 100 + 500 + ip * 10)->Fill(pInvMC, trCopy.getParam(ip) - trMC.getParam(ip));
            hm->getHisto2F(HMOffsets[2] + 1 * 100 + 500 + ip * 10)->Fill(etaMC, trCopy.getParam(ip) - trMC.getParam(ip));
          }
        }
      } else {
        trCopy = dummy;
      }
      (*outPerf) << "perfMatch" << "evID=" << int(mcTree->GetReadEntry()) << "mcTrID=" << mcInfoPair.first << "pdg=" << partMC.GetPdgCode() << "mcTr=" << trMC
                 << "mtcTr=" << trCopy << "chi2refit=" << chi2refit << "chi2match=" << chi2match << "nClones=" << nClones << "mcID=" << partID << "\n";
    }
    //-----------------------
  }
}

void checkVTTracking()
{
  // match tracks in the vtTrVec to MC info
  int ntr = vtTrVec->size();
  for (int itr = 0; itr < ntr; itr++) {
    const auto& tr = (*vtTrVec)[itr];
    int trMCID = tr.getParticleID();
    auto& mcInfo = mcTruth[std::abs(trMCID)];
    mcInfo.nVTTracks++;
    if (mcInfo.bestVTTrack < 0) {
      mcInfo.bestVTTrack = itr;
    } else {
      const auto& trPrev = (*vtTrVec)[mcInfo.bestVTTrack];
      if (trPrev.getNHits() < tr.getNHits() || (trPrev.getNHits() == tr.getNHits() && trPrev.getChi2() > tr.getChi2())) {
        mcInfo.bestVTTrack = itr;
      }
    }
  }
}

void checkMSTracking()
{
  // match tracks in the msTrVec to MC info
  int ntr = msTrVec->size();
  for (int itr = 0; itr < ntr; itr++) {
    const auto& tr = (*msTrVec)[itr];
    int trMCID = tr.getParticleID();
    auto& mcInfo = mcTruth[std::abs(trMCID)];
    mcInfo.nMSTracks++;
    if (mcInfo.bestMSTrack < 0) {
      mcInfo.bestMSTrack = itr;
    } else {
      const auto& trPrev = (*msTrVec)[mcInfo.bestMSTrack];
      if (trPrev.getNHits() < tr.getNHits() || (trPrev.getNHits() == tr.getNHits() && trPrev.getChi2() > tr.getChi2())) {
        mcInfo.bestMSTrack = itr;
      }
    }
  }
}

void checkMatching()
{
  // match tracks in the matchVec to MC info
  int ntr = matchVec->size();
  for (int itr = 0; itr < ntr; itr++) {
    const auto& tr = (*matchVec)[itr];
    int trMCID = tr.getParticleID();
    auto& mcInfo = mcTruth[std::abs(trMCID)];
    mcInfo.nMatches++;
    if (mcInfo.bestMatch < 0) {
      mcInfo.bestMatch = itr;
    } else {
      const auto& trPrev = (*matchVec)[mcInfo.bestMatch];
      if (trPrev.getChi2Match() > tr.getChi2Match()) {
        mcInfo.bestMSTrack = itr;
      }
    }
  }
}

void setupInputs(bool checkVTTracks, bool checkMSTracks, bool checkMatches, const char* dirSimu, const char* dirRec)
{
  // kine file needed to get the primary vertex position (temporary)
  auto fetchTree = [&](const char* fname, const char* tname, TTree*& tree) {
    auto f = TFile::Open(Form("%s/%s", dirSimu, fname));
    if (!f || f->IsZombie() || !(tree = (TTree*)f->Get(tname))) {
      exit(1);
    }
    nEntries = nEntries > 0 ? std::min(tree->GetEntries(), nEntries) : tree->GetEntries();
  };

  fetchTree(Form("%s/MCKine.root", dirSimu), "mckine", mcTree);
  mcTree->SetBranchAddress("tracks", &mcArr);

  if (checkVTTracks || checkMatches) {
    fetchTree(Form("%s/ClustersVerTel.root", dirSimu), "clustersVerTel", vtClTree);
    vtClTree->SetBranchAddress("VerTel", &vtClVec);

    fetchTree(Form("%s/TracksVerTel.root", dirSimu), "tracksVerTel", vtTrTree);
    vtTrTree->SetBranchAddress("VerTel", &vtTrVec);
  }

  if (checkMSTracks || checkMatches) {
    fetchTree(Form("%s/ClustersMuonSpec.root", dirSimu), "clustersMuonSpec", msClTree);
    msClTree->SetBranchAddress("MuonSpec", &msClVec);

    fetchTree(Form("%s/TracksMuonSpec.root", dirSimu), "tracksMuonSpec", msTrTree);
    msTrTree->SetBranchAddress("MuonSpec", &msTrVec);

    fetchTree(Form("%s/HitsMuonSpecModular.root", dirSimu), "hitsMuonSpecModular", msHitsTree);
    msHitsTree->SetBranchAddress("MuonSpecModular", &msHitsVec);
  }

  if (checkMatches) {
    fetchTree(Form("%s/TracksMatching.root", dirSimu), "tracksMatching", matchTree);
    matchTree->SetBranchAddress("Matching", &matchVec);
  }
}

void loadEventData(int i)
{
  if (mcTree)
    mcTree->GetEntry(i);
  if (vtClTree)
    vtClTree->GetEntry(i);
  if (msClTree)
    msClTree->GetEntry(i);
  if (vtTrTree)
    vtTrTree->GetEntry(i);
  if (msTrTree)
    msTrTree->GetEntry(i);
  if (msHitsTree)
    msHitsTree->GetEntry(i);
  if (matchTree)
    matchTree->GetEntry(i);
}

void buildReconstuctableMCInfo()
{
  auto consideMCPart = [&](int partID) {
    auto ent0 = mcTruth.emplace(partID, MCTrackInfo{});
    auto& ent = ent0.first->second;
    if (ent.status == MCTrackInfo::New) {
      const auto& part = (*mcArr)[partID];
      while (true) {
        ent.status = MCTrackInfo::Discard;
        auto ppdg = TDatabasePDG::Instance()->GetParticle(part.GetPdgCode());
        if (!ppdg || !ppdg->Charge()) {
          break;
        }
        auto eta = part.Eta();
        if (eta < cutEtaMin || eta > cutEtaMax) {
          break;
        }
        auto pt = part.Pt();
        if (pt < cutPtMin) {
          break;
        }
        // all checks passed
        // RS	printf("Accept track:%d/ev:%d from Z=%.3f time %.3f\n", partID, int(mcTree->GetReadEntry()), part.Vz(), part.T());
        ent.status = MCTrackInfo::Accept;
        ent.charge = ppdg->Charge() / 3;
        break;
      }
    }
    return ent.status == MCTrackInfo::Discard ? mcTruth.end() : ent0.first; // return iterator
  };

  mcTruth.clear();
  if (vtClVec) {
    for (int i = 0; i < (int)vtClVec->size(); i++) {
      auto mcID = (*vtClVec)[i].getParticleID();
      auto ent = consideMCPart(mcID);
      if (ent != mcTruth.end()) {
        ent->second.vtClID.push_back(i);
      }
    }
  }
  if (msClVec) {
    for (int i = 0; i < (int)msClVec->size(); i++) {
      auto mcID = (*msClVec)[i].getParticleID();
      auto ent = consideMCPart(mcID);
      if (ent != mcTruth.end()) {
        ent->second.msClID.push_back(i);
      }
    }
  }
  // sort clusters in layer increasing order, count layers
  for (auto& mci : mcTruth) {
    {
      auto& vtc = mci.second.vtClID;
      std::sort(vtc.begin(), vtc.end(), [](int i, int j) { return (*vtClVec)[i].getLayer() < (*vtClVec)[j].getLayer(); });
      int prevL = -1;
      for (auto i : vtc) {
        if ((*vtClVec)[i].getLayer() != prevL) {
          mci.second.nvtLr++;
          prevL = (*vtClVec)[i].getLayer();
        }
      }
    }
    {
      auto& msc = mci.second.msClID;
      std::sort(msc.begin(), msc.end(), [](int i, int j) { return (*msClVec)[i].getLayer() < (*msClVec)[j].getLayer(); });
      int prevL = -1;
      for (auto i : msc) {
        if ((*msClVec)[i].getLayer() != prevL) {
          mci.second.nmsLr++;
          prevL = (*msClVec)[i].getLayer();
        }
      }
      if (!msc.empty()) {
        const auto& cl = (*msClVec)[msc.front()];
        auto hit = (*msHitsVec)[cl.getHitID()];
        mci.second.momMS = std::array<float, 3>{(float)hit.getMomIn()[0], (float)hit.getMomIn()[1], (float)hit.getMomIn()[2]};
        mci.second.posMS = std::array<float, 3>{(float)hit.getPosIn()[0], (float)hit.getPosIn()[1], (float)hit.getPosIn()[2]};
      }
    }
  }
}
