#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PMCGenHeader.h"
#include "NA6PMCEventHeader.h"
#include "ConfigurableParam.h"
#include "NA6PLayoutParam.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <fairlogger/Logger.h>
#include <cstdio>
#include <bitset>

#endif

TChain* loadUserChain(const char* inpData, const char* chName, bool recursive = false);
void addRootFileDirectory(TChain* chain, const char* flName, const char* chName, bool recursive = false);

uint64_t getACTSgeomID(uint64_t layer, uint64_t sensitive)
{
  // layer numbered as 1, 2, 3 ...
  // sensitive numbered as 1,2,3,4 for VT, 1 for MS
  uint64_t volume = 1;
  uint64_t boundar = 0;
  uint64_t approach = 0;
  uint64_t channel = 0;
  uint64_t id = 0;
  layer *= 2;
  id += volume << 56;
  id += boundar << 48;
  id += layer << 36;
  id += approach << 28;
  id += sensitive << 8;
  id += channel;
  return id;
}

void convert2ACTS(int iev,
                  const std::vector<TParticle>& mcParts,
                  const std::vector<NA6PVerTelHit>& vtHits,
                  const std::vector<NA6PMuonSpecHit>& msHits);

void convert2ACTS(const std::string& dirname = "./",
                  const std::string& fnameMCKin = "MCKine.root",
                  const std::string& fnameVTHits = "HitsVerTel.root",
                  const std::string& fnameMCHits = "HitsMuonSpec.root",
                  const std::string& confOpts = "")
{
  std::string dirnameL = dirname;
  if (dirnameL.empty()) {
    dirnameL = ".";
  } else if (dirnameL.back() != '/') {
    dirnameL += '/';
  }
  if (!confOpts.empty()) {
    na6p::conf::ConfigurableParam::updateFromString(confOpts);
  }
  NA6PMCEventHeader mcHeader, *mcHeaderPtr = &mcHeader;
  std::vector<TParticle> mcParts, *mcPartsPtr = &mcParts;
  std::vector<NA6PMuonSpecHit> msHits, *msHitsPtr = &msHits;
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;

  TChain* treeKin = loadUserChain(fmt::format("{}{}", dirnameL, fnameMCKin).c_str(), "mckine");
  TChain* treeVTH = loadUserChain(fmt::format("{}{}", dirnameL, fnameVTHits).c_str(), "hitsVerTel");
  TChain* treeMSH = loadUserChain(fmt::format("{}{}", dirnameL, fnameMCHits).c_str(), "hitsMuonSpec");
  int nent = 0;
  if (!treeKin || !treeVTH || !treeMSH || !(nent = treeKin->GetEntries()) || treeVTH->GetEntries() != nent || treeMSH->GetEntries() != nent) {
    LOGP(fatal, "Problem in opening input file or empty input");
  }
  treeKin->SetBranchAddress("header", &mcHeaderPtr);
  treeKin->SetBranchAddress("tracks", &mcPartsPtr);
  treeVTH->SetBranchAddress("VerTel", &vtHitsPtr);
  treeMSH->SetBranchAddress("MuonSpec", &msHitsPtr);
  for (int iev = 0; iev < nent; iev++) {
    treeKin->GetEntry(iev);
    treeVTH->GetEntry(iev);
    treeMSH->GetEntry(iev);
    LOGP(info, "Event#{} Vtx:[{:.2f},{:.2f},{:.2f}] {} Tracks ({} primaries), Hits: VT: {} MS: {}", iev, mcHeader.getVX(), mcHeader.getVY(), mcHeader.getVZ(), mcHeader.getNTracks(), mcHeader.getNPrimaries(), vtHits.size(), msHits.size());

    convert2ACTS(iev, mcParts, vtHits, msHits);
  }
}

TChain* loadUserChain(const char* inpData, const char* chName, bool recursive)
{
  if (!inpData) {
    return nullptr;
  }
  TChain* chain = new TChain(chName);
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    addRootFileDirectory(chain, inpData, chName, recursive);
  } else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      LOGP(error, "Failed on input filename {}", inpData);
      return nullptr;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while (!flName.IsNull()) {
      flName = flName.Strip(TString::kBoth, ' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {
        flName.ReadLine(inpf);
        continue;
      }
      flName = flName.Strip(TString::kBoth, ',');
      flName = flName.Strip(TString::kBoth, '"');
      addRootFileDirectory(chain, flName.Data(), chName, recursive);
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n < 1) {
    LOGP(error, "Obtained chain is empty");
    return nullptr;
  } else
    LOGP(info, "Opened {} chain with {} entries", chName, n);
  return chain;
}

void addRootFileDirectory(TChain* chain, const char* flName, const char* chName, bool recursive)
{
  if (!recursive) {
    printf("Adding %s\n", flName);
    chain->AddFile(flName);
  } else {
    TFile* f = TFile::Open(flName);
    if (!f || f->IsZombie()) {
      LOGP(error, "Failed to open {}", flName);
      delete f;
      return;
    }
    auto* lstD = f->GetListOfKeys();
    TIter iter(lstD->MakeIterator());
    while (TObject* obj = iter()) {
      TKey* theKey = (TKey*)obj;
      TString kname = theKey->GetClassName();
      if (!kname.BeginsWith("TDirectory")) {
        continue;
      }
      if (f->Get(Form("%s/%s", theKey->GetName(), chName))) {
        LOGP(info, "Adding {}/{}", flName, theKey->GetName());
        chain->AddFile(flName, TTree::kMaxEntries, Form("%s/%s", theKey->GetName(), chName));
      }
    }
    delete f;
  }
}

void convert2ACTS(int iev,
                  const std::vector<TParticle>& mcParts,
                  const std::vector<NA6PVerTelHit>& vtHits,
                  const std::vector<NA6PMuonSpecHit>& msHits)
{
  std::string eventFileName = fmt::format("event{:09}-particles.csv", iev);
  std::string primFileName = fmt::format("event{:09}-primaries.csv", iev);
  std::string hitsFileName = fmt::format("event{:09}-hits.csv", iev);
  std::FILE* fpcsvKin = std::fopen(eventFileName.c_str(), "w");
  std::FILE* fpcsvPrim = std::fopen(primFileName.c_str(), "w");

  fprintf(fpcsvKin, "particle_id,particle_type,process,vx,vy,vz,vt,px,py,pz,m,q\n");
  fprintf(fpcsvPrim, "particle_id,particle_type,process,vx,vy,vz,vt,px,py,pz,m,q\n");
  int decCount = 0;
  std::vector<uint64_t> doneBC;
  doneBC.resize(mcParts.size());
  std::vector<uint64_t> secV;
  secV.resize(mcParts.size());
  printf("Loop on %d particles\n", (int)mcParts.size());
  for (int ipartTop = 0; ipartTop < (int)mcParts.size(); ipartTop++) {
    // clang-format off
    auto getBarCode = [&iev, &ipartTop](int gen, int sub, int isv)
    {
      uint64_t b =
	( (uint64_t(iev+1)    & ((0x1<<12)-1)) << ( 12 + 16 + 8 + 16 ) ) |
	( (uint64_t(isv)      & ((0x1<<12)-1)) << (      16 + 8 + 16 ) ) |
	( (uint64_t(ipartTop) & ((0x1<<16)-1)) << (           8 + 16 ) ) |
	( (uint64_t(gen)      & ((0x1<<8)-1))  << (               16 ) ) |
	uint64_t(sub);
      return b;
    };
    // clang-format on

    auto storePart = [iev, ipartTop, fpcsvKin, fpcsvPrim, &doneBC, &secV, &decCount, &mcParts, getBarCode](int ip, int gen, int sub, int imoth, const auto& self) -> void {
      if (doneBC[ip]) {
        //	printf("Particle %d already stored\n",ip);
        return;
      }
      const auto& part = mcParts[ip];
      const auto pdgPart = part.GetPDG();
      double fcmtomm = 10.f;
      if (!pdgPart) {
        LOGP(info, "Skipping unknown particle {}", part.GetPdgCode());
        return;
      }
      doneBC[ip] = getBarCode(gen, sub, imoth);

      fprintf(fpcsvKin, "%lu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", doneBC[ip], std::abs(part.GetPdgCode()), 0,
              part.Vx() * fcmtomm, part.Vy() * fcmtomm, part.Vz() * fcmtomm, part.T(),
              part.Px(), part.Py(), part.Pz(), part.GetCalcMass(), pdgPart->Charge() / 3);
      if (part.IsPrimary()) {
        fprintf(fpcsvPrim, "%lu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", doneBC[ip], std::abs(part.GetPdgCode()), 0,
                part.Vx() * fcmtomm, part.Vy() * fcmtomm, part.Vz() * fcmtomm, part.T(),
                part.Px(), part.Py(), part.Pz(), part.GetCalcMass(), pdgPart->Charge() / 3);
      }
      int dtF = part.GetFirstDaughter();
      int dtL = part.GetLastDaughter();
      if (dtF > -1) {
        gen++;
        decCount++;
        secV[ip] = decCount;
        int subG = 0;
        for (int i = dtF; i <= dtL; i++) {
          self(i, gen, subG++, secV[ip], self);
        }
      }
    };
    storePart(ipartTop, 0, 0, 0, storePart);
  }
  std::fclose(fpcsvKin);
  std::fclose(fpcsvPrim);

  std::FILE* fpcsvHits = std::fopen(hitsFileName.c_str(), "w");
  fprintf(fpcsvHits, "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");

  auto getGeomID = [](int sensID) {
    const auto& param = NA6PLayoutParam::Instance();
    uint64_t b = 0;
    return b;
  };
  double fcmtomm = 10.f;
  for (const auto& hit : vtHits) {
    int det = hit.getDetectorID();
    double m = mcParts[hit.getTrackID()].GetCalcMass();
    double en = std::sqrt(hit.getPX() * hit.getPX() + hit.getPY() * hit.getPY() + hit.getPZ() * hit.getPZ() + m * m);
    auto particle_id = doneBC[hit.getTrackID()];
    if (particle_id == 0) {
      int idP = hit.getTrackID();
      const auto& part = mcParts[idP];
      const auto pdgPart = part.GetPDG();
      printf("Hit w/o particle: %lu %d pdgPart=%x, pdgcode=%d\n", particle_id, hit.getTrackID(), pdgPart, part.GetPdgCode());
    }
    int layer = det / 4 + 1; // layer from 1 to 5
    int sens = det % 4 + 1;  // sens from 1 to 4
    uint64_t geometry_id = getACTSgeomID(layer, sens);
    fprintf(fpcsvHits, "%lu,%lu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", particle_id, geometry_id,
            hit.getX() * fcmtomm, hit.getY() * fcmtomm, hit.getZ() * fcmtomm, hit.getTime(),
            hit.getPX(), hit.getPY(), hit.getPZ(), en,
            0., 0., 0., 0., 0);
  }
  for (const auto& hit : msHits) {
    int det = hit.getDetectorID();
    double m = mcParts[hit.getTrackID()].GetCalcMass();
    double en = std::sqrt(hit.getPX() * hit.getPX() + hit.getPY() * hit.getPY() + hit.getPZ() * hit.getPZ() + m * m);
    auto particle_id = doneBC[hit.getTrackID()];
    int layer = det + NA6PLayoutParam::Instance().nVerTelPlanes + 1;
    int sens = 1;
    uint64_t geometry_id = getACTSgeomID(layer, sens);
    fprintf(fpcsvHits, "%lu,%lu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", particle_id, geometry_id,
            hit.getX() * fcmtomm, hit.getY() * fcmtomm, hit.getZ() * fcmtomm, hit.getTime(),
            hit.getPX(), hit.getPY(), hit.getPZ(), en,
            0., 0., 0., 0., 0);
  }
  std::fclose(fpcsvHits);
}
