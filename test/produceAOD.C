#if !defined(__CINT__) || defined(__MAKECINT__)
#include <vector>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>

#include "NA6PTrack.h"
#include "NA6PVertex.h"
#include "NA6PMatching.h"
#endif

// root -b -q produceAOD.C++g

void produceAOD(
const char* fileNameVerTel = "VerticesVerTel.root",
const char* fileNameTracksVerTel = "TracksVerTel.root",
const char* fileNameTracksMuonSpec = "TracksMuonSpec.root",
const char* fileNameTracksMatching = "TracksMatching.root",
const char* fileNameMC = "MCKine.root",
const char* suffix="")
{
  enum class trackType : uint8_t {
    kVT = 0, // Vertex Telescope standalone
    kMS = 1, // Muon Spectrometer standalone
    kMatched = 2, // matched between VT and MS
  };

  printf("Please make sure that you enabled mtPropagateMatchedTracksToPV=true and vtPropagateTracksToPV=true in reco. params.\n");

  TFile* rootfileVerTel = TFile::Open(fileNameVerTel, "READ");
  TFile* rootfileTracksVerTel = TFile::Open(fileNameTracksVerTel, "READ");
  TFile* rootfileTracksMuonSpec = TFile::Open(fileNameTracksMuonSpec, "READ");
  TFile* rootfileTracksMatching = TFile::Open(fileNameTracksMatching, "READ");
  TFile* rootfileMC = TFile::Open(fileNameMC, "READ");

  TTree* treeVerTel = (TTree*)rootfileVerTel->Get("verticesVerTel");
  std::vector<NA6PVertex>* vecVerTel = nullptr;
  treeVerTel->SetBranchAddress("VerTel", &vecVerTel);

  TTree* treeTracksVerTel = (TTree*)rootfileTracksVerTel->Get("tracksVerTel");
  std::vector<NA6PTrack>* vecTrackVerTel = nullptr;
  treeTracksVerTel->SetBranchAddress("VerTel", &vecTrackVerTel);

  TTree* treeTracksMuonSpec = (TTree*)rootfileTracksMuonSpec->Get("tracksMuonSpec");
  std::vector<NA6PTrack>* vecTrackMuonSpec = nullptr;
  treeTracksMuonSpec->SetBranchAddress("MuonSpec", &vecTrackMuonSpec);

  TTree* treeTracksMatching = (TTree*)rootfileTracksMatching->Get("tracksMatching");
  std::vector<NA6PMatch>* vecTrackMatching = nullptr;
  treeTracksMatching->SetBranchAddress("Matching", &vecTrackMatching);

  float mVX = 0, mVY=0, mVZ=0;
  TTree* treeMC = (TTree*)rootfileMC->Get("mckine");
  std::vector<TParticle>* vecMCParticle = nullptr;
  treeMC->SetBranchAddress("mVX", &mVX);
  treeMC->SetBranchAddress("mVY", &mVY);
  treeMC->SetBranchAddress("mVZ", &mVZ);
  treeMC->SetBranchAddress("tracks", &vecMCParticle);

  TFile *outfile = new TFile(Form("na6pAOD%s.root", suffix), "RECREATE");
  outfile->SetCompressionSettings(505);

  TTree *outTreeCollision = new TTree("collisions", "collision information");
  uint32_t mcCollisionId = 0;
  float posX=-999, posY=-999, posZ=-999, chi2 = -999;
  uint16_t nContributors = 0;
  short vertexType = -1;
  std::vector<int> pvContributorIds;
  outTreeCollision->Branch("posX", &posX, "posX/F");
  outTreeCollision->Branch("posY", &posY, "posY/F");
  outTreeCollision->Branch("posZ", &posZ, "posZ/F");
  outTreeCollision->Branch("chi2", &chi2, "chi2/F");
  outTreeCollision->Branch("nContributors", &nContributors, "nContributors/s");
  outTreeCollision->Branch("pvContributorIds", &pvContributorIds);
  outTreeCollision->Branch("vertexType", &vertexType, "vertexType/S");
  outTreeCollision->Branch("mcCollisionId", &mcCollisionId, "mcCollisionId/i");

  TTree *outTreeCollisionCov = new TTree("collisionsCov", "collision covariance information");
  float covXX = -999, covXY = -999, covXZ= -999, covYY = -999, covYZ = -999, covZZ = -999;
  outTreeCollisionCov->Branch("cXX", &covXX, "cXX/F");
  outTreeCollisionCov->Branch("cYY", &covYY, "cYY/F");
  outTreeCollisionCov->Branch("cZZ", &covZZ, "cZZ/F");
  outTreeCollisionCov->Branch("cXY", &covXY, "cXY/F");
  outTreeCollisionCov->Branch("cYZ", &covYZ, "cYZ/F");
  outTreeCollisionCov->Branch("cXZ", &covXZ, "cXZ/F");

  uint32_t collisionId=0;
  uint8_t trackType=0;
  int8_t sign = 0;
  float chi2Norm = 0;
  uint8_t nClusters = 0;
  float impParX=0, impParY=0;
  float x = 0, y = 0, tx = 0, ty = 0, q2Pxz = 0;
  float cXX=0, cYY=0, cYX=0,
  cTxX = 0, cTxY = 0, cTx2 = 0,
  cTyX = 0, cTyY = 0, cTyTx = 0, cTy2 = 0,
  cQ2PxzX = 0, cQ2PxzY = 0, cQ2PxzTx = 0, cQ2PzTy = 0, cQ2Pxz2 = 0;

  uint32_t mcParticleId=0;
  TTree *outTreeTrackVT = new TTree("tracksVT", "VT track information");
  outTreeTrackVT->Branch("collisionId", &collisionId, "collisionId/i");
  outTreeTrackVT->Branch("x", &x, "x/F");
  outTreeTrackVT->Branch("y", &y, "y/F");
  outTreeTrackVT->Branch("tx", &tx, "tx/F");
  outTreeTrackVT->Branch("ty", &ty, "ty/F");
  outTreeTrackVT->Branch("q2Pxz", &q2Pxz, "q2Pxz/F");
  outTreeTrackVT->Branch("nClusters", &nClusters, "nClusters/b");
  outTreeTrackVT->Branch("chi2Norm", &chi2Norm, "chi2Norm/F");
  outTreeTrackVT->Branch("mcParticleId", &mcParticleId, "mcParticleId/i");

  TTree *outTreeTrackCovVT = new TTree("tracksCovVT", "VT track covariance matrix information"); // call AddFriends at analysis level.
  outTreeTrackCovVT->Branch("cXX", &cXX, "cXX/F");
  outTreeTrackCovVT->Branch("cYY", &cYY, "cYY/F");
  outTreeTrackCovVT->Branch("cYX", &cYX, "cYX/F");
  outTreeTrackCovVT->Branch("cTxX", &cTxX, "cTxX/F");
  outTreeTrackCovVT->Branch("cTxY", &cTxY, "cTxY/F");
  outTreeTrackCovVT->Branch("cTx2", &cTx2, "cTx2/F");
  outTreeTrackCovVT->Branch("cTyX", &cTyX, "cTyX/F");
  outTreeTrackCovVT->Branch("cTyY", &cTyY, "cTyY/F");
  outTreeTrackCovVT->Branch("cTyTx", &cTyTx, "cTyTx/F");
  outTreeTrackCovVT->Branch("cTy2", &cTy2, "cTy2/F");
  outTreeTrackCovVT->Branch("cQ2PxzX", &cQ2PxzX, "cQ2PxzX/F");
  outTreeTrackCovVT->Branch("cQ2PxzY", &cQ2PxzY, "cQ2PxzY/F");
  outTreeTrackCovVT->Branch("cQ2PxzTx", &cQ2PxzTx, "cQ2PxzTx/F");
  outTreeTrackCovVT->Branch("cQ2PzTy", &cQ2PzTy, "cQ2PzTy/F");
  outTreeTrackCovVT->Branch("cQ2Pxz2", &cQ2Pxz2, "cQ2Pxz2/F");

  TTree *outTreeTrackMS = new TTree("tracksMS", "MS track information");
  outTreeTrackMS->Branch("collisionId", &collisionId, "collisionId/i");
  outTreeTrackMS->Branch("x", &x, "x/F");
  outTreeTrackMS->Branch("y", &y, "y/F");
  outTreeTrackMS->Branch("tx", &tx, "tx/F");
  outTreeTrackMS->Branch("ty", &ty, "ty/F");
  outTreeTrackMS->Branch("q2Pxz", &q2Pxz, "q2Pxz/F");
  outTreeTrackMS->Branch("nClusters", &nClusters, "nClusters/b");
  outTreeTrackMS->Branch("chi2Norm", &chi2Norm, "chi2Norm/F");
  outTreeTrackMS->Branch("mcParticleId", &mcParticleId, "mcParticleId/i");

  TTree *outTreeTrackCovMS = new TTree("tracksCovMS", "MS track covariance matrix information"); // call AddFriends at analysis level.
  outTreeTrackCovMS->Branch("cXX", &cXX, "cXX/F");
  outTreeTrackCovMS->Branch("cYY", &cYY, "cYY/F");
  outTreeTrackCovMS->Branch("cYX", &cYX, "cYX/F");
  outTreeTrackCovMS->Branch("cTxX", &cTxX, "cTxX/F");
  outTreeTrackCovMS->Branch("cTxY", &cTxY, "cTxY/F");
  outTreeTrackCovMS->Branch("cTx2", &cTx2, "cTx2/F");
  outTreeTrackCovMS->Branch("cTyX", &cTyX, "cTyX/F");
  outTreeTrackCovMS->Branch("cTyY", &cTyY, "cTyY/F");
  outTreeTrackCovMS->Branch("cTyTx", &cTyTx, "cTyTx/F");
  outTreeTrackCovMS->Branch("cTy2", &cTy2, "cTy2/F");
  outTreeTrackCovMS->Branch("cQ2PxzX", &cQ2PxzX, "cQ2PxzX/F");
  outTreeTrackCovMS->Branch("cQ2PxzY", &cQ2PxzY, "cQ2PxzY/F");
  outTreeTrackCovMS->Branch("cQ2PxzTx", &cQ2PxzTx, "cQ2PxzTx/F");
  outTreeTrackCovMS->Branch("cQ2PzTy", &cQ2PzTy, "cQ2PzTy/F");
  outTreeTrackCovMS->Branch("cQ2Pxz2", &cQ2Pxz2, "cQ2Pxz2/F");

  float chi2Match = 0, chi2Refit = 0;
  uint32_t vtId=0, msId=0;
  uint32_t mcParticleIdVT=0, mcParticleIdMS=0;
  TTree *outTreeTrackMatched = new TTree("tracksMatched", "Matched track information");
  outTreeTrackMatched->Branch("collisionId", &collisionId, "collisionId/i");
  outTreeTrackMatched->Branch("x", &x, "x/F");
  outTreeTrackMatched->Branch("y", &y, "y/F");
  outTreeTrackMatched->Branch("tx", &tx, "tx/F");
  outTreeTrackMatched->Branch("ty", &ty, "ty/F");
  outTreeTrackMatched->Branch("q2Pxz", &q2Pxz, "q2Pxz/F");
  outTreeTrackMatched->Branch("nClusters", &nClusters, "nClusters/b"); // total number of clusters
  outTreeTrackMatched->Branch("chi2Match", &chi2Match, "chi2Match/F");
  outTreeTrackMatched->Branch("chi2Refit", &chi2Refit, "chi2Refit/F");
  outTreeTrackMatched->Branch("vtId", &vtId, "vtId/i");
  outTreeTrackMatched->Branch("msId", &msId, "msId/i");
  outTreeTrackMatched->Branch("mcParticleId", &mcParticleId, "mcParticleId/i");
  outTreeTrackMatched->Branch("mcParticleIdVT", &mcParticleIdVT, "mcParticleIdVT/i");
  outTreeTrackMatched->Branch("mcParticleIdMS", &mcParticleIdMS, "mcParticleIdMS/i");

  TTree *outTreeTrackCovMatched = new TTree("tracksCovMatched", "Matched track covariance matrix information"); // call AddFriends at analysis level.
  outTreeTrackCovMatched->Branch("cXX", &cXX, "cXX/F");
  outTreeTrackCovMatched->Branch("cYY", &cYY, "cYY/F");
  outTreeTrackCovMatched->Branch("cYX", &cYX, "cYX/F");
  outTreeTrackCovMatched->Branch("cTxX", &cTxX, "cTxX/F");
  outTreeTrackCovMatched->Branch("cTxY", &cTxY, "cTxY/F");
  outTreeTrackCovMatched->Branch("cTx2", &cTx2, "cTx2/F");
  outTreeTrackCovMatched->Branch("cTyX", &cTyX, "cTyX/F");
  outTreeTrackCovMatched->Branch("cTyY", &cTyY, "cTyY/F");
  outTreeTrackCovMatched->Branch("cTyTx", &cTyTx, "cTyTx/F");
  outTreeTrackCovMatched->Branch("cTy2", &cTy2, "cTy2/F");
  outTreeTrackCovMatched->Branch("cQ2PxzX", &cQ2PxzX, "cQ2PxzX/F");
  outTreeTrackCovMatched->Branch("cQ2PxzY", &cQ2PxzY, "cQ2PxzY/F");
  outTreeTrackCovMatched->Branch("cQ2PxzTx", &cQ2PxzTx, "cQ2PxzTx/F");
  outTreeTrackCovMatched->Branch("cQ2PzTy", &cQ2PzTy, "cQ2PzTy/F");
  outTreeTrackCovMatched->Branch("cQ2Pxz2", &cQ2Pxz2, "cQ2Pxz2/F");

  float posXMC=-999, posYMC=-999, posZMC=-999;
  TTree *outTreeCollisionMC = new TTree("collisionsMC", "mc collision information");
  outTreeCollisionMC->Branch("posX", &posXMC, "posX/F");
  outTreeCollisionMC->Branch("posY", &posYMC, "posY/F");
  outTreeCollisionMC->Branch("posZ", &posZMC, "posZ/F");

  TTree *outTreeParticleMC = new TTree("particlesMC", "mc particle information");
  float pt=0, eta=0, phi=0, e = 0;
  float vx = 0, vy=0, vz=0;
  int pdgCode = 0, status = 0;
  int firstMotherId = -1, secondMotherId = -1;
  int firstDaughterId = -1, lastDaughterId = -1;
  outTreeParticleMC->Branch("pt", &pt, "pt/F");
  outTreeParticleMC->Branch("eta", &eta, "eta/F");
  outTreeParticleMC->Branch("phi", &phi, "phi/F");
  outTreeParticleMC->Branch("e", &e, "e/F"); // energy
  outTreeParticleMC->Branch("vx", &vx, "vx/F"); // production vertex x
  outTreeParticleMC->Branch("vy", &vy, "vy/F"); // production vertex y
  outTreeParticleMC->Branch("vz", &vz, "vz/F"); // production vertex z
  outTreeParticleMC->Branch("mcCollisionId", &mcCollisionId, "mcCollisionId/i");
  outTreeParticleMC->Branch("pdgCode", &pdgCode, "pdgCode/I");
  outTreeParticleMC->Branch("status", &status, "status/I");
  outTreeParticleMC->Branch("firstMotherId", &firstMotherId, "firstMotherId/I");
  outTreeParticleMC->Branch("secondMotherId", &secondMotherId, "secondMotherId/I");
  outTreeParticleMC->Branch("firstDaughterId", &firstDaughterId, "firstDaughterId/I");
  outTreeParticleMC->Branch("lastDaughterId", &lastDaughterId, "lastDaughterId/I");

  for (int i = 0; i < treeMC->GetEntries(); i++) {
    treeMC->GetEntry(i);
    treeVerTel->GetEntry(i);
    treeTracksVerTel->GetEntry(i);
    treeTracksMuonSpec->GetEntry(i);
    treeTracksMatching->GetEntry(i);
    mcCollisionId = i; // this will be stored in treeCol.

    posXMC = mVX;
    posYMC = mVY;
    posZMC = mVZ;
    outTreeCollisionMC->Fill();

    for (const auto& p : *vecMCParticle) {
      mcCollisionId = i;
      pt = p.Pt();
      eta = p.Eta();
      phi = p.Phi() < 0 ? p.Phi() + 2 * M_PI : p.Phi();
      e = p.Energy();
      vx = p.Vx();
      vy = p.Vy();
      vz = p.Vz();
      pdgCode = p.GetPdgCode();
      status = p.GetStatusCode();
      firstMotherId = p.GetFirstMother();
      secondMotherId = p.GetSecondMother();
      firstDaughterId = p.GetFirstDaughter();
      lastDaughterId = p.GetLastDaughter();
      outTreeParticleMC->Fill();
    } // end of TParticle loop per mc collision

    // printf("mcCollisionId = %d, vecVerTel->size() = %zu, vecTrackVerTel->size() = %zu, vecTrackMuonSpec->size() = %zu, vecTrackMatching->size() = %zu\n", i, vecVerTel->size(), vecTrackVerTel->size(), vecTrackMuonSpec->size(), vecTrackMatching->size());

    for (const auto& vtx : *vecVerTel) {
      posX = vtx.getX();
      posY = vtx.getY();
      posZ = vtx.getZ();
      chi2 = vtx.getChi2();
      nContributors = vtx.getNContributors();
      pvContributorIds = vtx.getTrackIDs();
      vertexType = vtx.getVertexType();
      outTreeCollision->Fill();

      covXX = vtx.getSigmaX2();
      covYY = vtx.getSigmaY2();
      covZZ = vtx.getSigmaZ2();
      covXY = vtx.getSigmaXY();
      covYZ = vtx.getSigmaYZ();
      covXZ = vtx.getSigmaXZ();
      outTreeCollisionCov->Fill();
    } // end of rec. vertex loop

    for (const auto& track : *vecTrackVerTel) {
      collisionId = mcCollisionId; // temporary index. As of 20260709, rec. tracks are not associated to rec. vertex.
      x = track.getX();
      y = track.getY();
      tx = track.getTx();
      ty = track.getTy();
      q2Pxz = track.getQ2Pxz();
      nClusters = track.getNHits();
      chi2Norm = track.getChi2Norm();
      mcParticleId = track.getParticleID();
      outTreeTrackVT->Fill();

      cXX      = track.getSigmaX2();
      cYY      = track.getSigmaY2();
      cYX      = track.getSigmaYX();
      cTxX     = track.getSigmaTxX();
      cTxY     = track.getSigmaTxY();
      cTx2     = track.getSigmaTx2();
      cTyX     = track.getSigmaTyX();
      cTyY     = track.getSigmaTyY();
      cTyTx    = track.getSigmaTyTx();
      cTy2     = track.getSigmaTy2();
      cQ2PxzX  = track.getSigmaQ2PxzX();
      cQ2PxzY  = track.getSigmaQ2PxzY();
      cQ2PxzTx = track.getSigmaQ2PxzTx();
      cQ2PzTy  = track.getSigmaQ2PzTy();
      cQ2Pxz2  = track.getSigmaQ2Pxz2();
      outTreeTrackCovVT->Fill();
    } // end of track loop reconstructed as VTsa.

    for (const auto& track : *vecTrackMuonSpec) {
      collisionId = mcCollisionId; // temporary index. As of 20260709, rec. tracks are not associated to rec. vertex.
      x = track.getX();
      y = track.getY();
      tx = track.getTx();
      ty = track.getTy();
      q2Pxz = track.getQ2Pxz();
      nClusters = track.getNHits();
      chi2Norm = track.getChi2Norm();
      mcParticleId = track.getParticleID();
      outTreeTrackMS->Fill();

      cXX      = track.getSigmaX2();
      cYY      = track.getSigmaY2();
      cYX      = track.getSigmaYX();
      cTxX     = track.getSigmaTxX();
      cTxY     = track.getSigmaTxY();
      cTx2     = track.getSigmaTx2();
      cTyX     = track.getSigmaTyX();
      cTyY     = track.getSigmaTyY();
      cTyTx    = track.getSigmaTyTx();
      cTy2     = track.getSigmaTy2();
      cQ2PxzX  = track.getSigmaQ2PxzX();
      cQ2PxzY  = track.getSigmaQ2PxzY();
      cQ2PxzTx = track.getSigmaQ2PxzTx();
      cQ2PzTy  = track.getSigmaQ2PzTy();
      cQ2Pxz2  = track.getSigmaQ2Pxz2();
      outTreeTrackCovMS->Fill();
    } // end of track loop reconstructed as MSsa.

    for (const auto& track : *vecTrackMatching) {
      collisionId = mcCollisionId; // temporary index. As of 20260709, rec. tracks are not associated to rec. vertex.
      x = track.getX();
      y = track.getY();
      tx = track.getTx();
      ty = track.getTy();
      q2Pxz = track.getQ2Pxz();
      nClusters = track.getNClusters();
      chi2Match = track.getChi2Match();
      chi2Refit = track.getChi2Refit();
      vtId = track.getIndexVT();
      msId = track.getIndexMS();
      mcParticleId = track.getParticleID();
      mcParticleIdVT = track.getParticleIDVT();
      mcParticleIdMS = track.getParticleIDMS();
      outTreeTrackMatched->Fill();

      cXX      = track.getSigmaX2();
      cYY      = track.getSigmaY2();
      cYX      = track.getSigmaYX();
      cTxX     = track.getSigmaTxX();
      cTxY     = track.getSigmaTxY();
      cTx2     = track.getSigmaTx2();
      cTyX     = track.getSigmaTyX();
      cTyY     = track.getSigmaTyY();
      cTyTx    = track.getSigmaTyTx();
      cTy2     = track.getSigmaTy2();
      cQ2PxzX  = track.getSigmaQ2PxzX();
      cQ2PxzY  = track.getSigmaQ2PxzY();
      cQ2PxzTx = track.getSigmaQ2PxzTx();
      cQ2PzTy  = track.getSigmaQ2PzTy();
      cQ2Pxz2  = track.getSigmaQ2Pxz2();
      outTreeTrackCovMatched->Fill();
    } // end of track loop after matching.

  } // end of mc collision loop

  outfile->WriteTObject(outTreeCollision);
  outfile->WriteTObject(outTreeCollisionCov);
  outfile->WriteTObject(outTreeTrackVT);
  outfile->WriteTObject(outTreeTrackCovVT);
  outfile->WriteTObject(outTreeTrackMS);
  outfile->WriteTObject(outTreeTrackCovMS);
  outfile->WriteTObject(outTreeTrackMatched);
  outfile->WriteTObject(outTreeTrackCovMatched);
  outfile->WriteTObject(outTreeCollisionMC);
  outfile->WriteTObject(outTreeParticleMC);
  outfile->Close();

  rootfileVerTel->Close();
  rootfileTracksVerTel->Close();
  rootfileTracksMuonSpec->Close();
  rootfileTracksMatching->Close();
  rootfileMC->Close();
}

