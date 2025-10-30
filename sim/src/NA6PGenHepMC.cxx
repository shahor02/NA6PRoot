// NA6PCCopyright

#include "NA6PGenHepMC.h"
#include "NA6PMCStack.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/FourVector.h"
#include <TDatabasePDG.h>
#include <TMath.h>
#include <cmath>

NA6PGenHepMC::NA6PGenHepMC(const std::string& name, const std::string& filname, bool storeDecayed) : NA6PGenerator(name), mFileName(filname), mStoreDecayedPrimaries(storeDecayed)
{
  LOGP(info, "create NA6PGenHepMC {} {}, decayed primaries will {}be stored", name, filname, storeDecayed ? "" : "not ");
}

long NA6PGenHepMC::canGenerateMaxEvents() const
{
  if (mNEvInTree < 0 || isInitDone()) {
    LOGP(fatal, "Initialization was not done");
  }
  return mNEvInTree;
}

void NA6PGenHepMC::init()
{
  LOGP(info, "NA6PGenHepMC init");
  if (isInitDone()) {
    return;
  }
  try {
    auto fl = TFile::Open(mFileName.c_str());
    assert(fl);
    TTree* tree = (TTree*)fl->Get(HEPTreeName.c_str());
    assert(tree);
    mNEvInTree = tree->GetEntries();
    delete tree;
    fl->Close();
  } catch (...) {
    LOGP(fatal, "Failed to extract the HEPMC tree {} from file {}", HEPTreeName, mFileName);
  }

  mHEPRootFileReader = std::make_unique<HepMC3::ReaderRootTree>(mFileName);
  if (mHEPRootFileReader->failed()) {
    LOGP(fatal, "No HepMC input file found {}", mFileName);
  }
  NA6PGenerator::init();
}

void NA6PGenHepMC::generate()
{
  if (mReadEvents >= mNEvInTree) {
    LOGP(error, "Something is wrong: event {} beyond the available range [0:{}] is requested, wrapping to 0", mReadEvents, mNEvInTree - 1);
    mReadEvents = 0;
  }
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  auto mcHead = getStack()->getEventHeader();
  double xv = mcHead->getVX();
  double yv = mcHead->getVY();
  double zv = mcHead->getVZ();
  std::unordered_map<int, int> partID2StoreID;

  HepMC3::GenEvent evt;
  mHEPRootFileReader->read_event(evt); // Read one event
  if (mHEPRootFileReader->failed()) {
    LOGP(info, " End of file reached . Exit .\n");
    return;
  }

  int nparticlesgood = 0;
  int pdgCode, status;
  float phi, pt, pZ, en, sn, cs;
  float vx, vy, vz, tof; // tof in mm/c
  int tobetracked = 0;

  std::vector<int> partons = {1, 2, 3, 4, 5, 6, 21}; // to remove quarks and gluons
  int cnt = 0;
  // Loop on particles
  for (const auto p : evt.particles()) {
    cnt++;
    auto isParton = [](int pabs) { return (pabs >= 1 & pabs <= 6) || pabs == 21; };
    auto isDiquark = [](int pabs) { return pabs > 1000 && pabs < 10000 && ((pabs / 10) % 10) == 0; };
    auto isColorState = [](int pabs) { return pabs > 1e6; };

    // Upstream and downstream vertices
    auto prod_vtx = p->production_vertex();
    auto end_vtx = p->end_vertex();
    auto pidabs = std::abs(p->pid());

    if (prod_vtx) {
      if (isParton(pidabs) ||  // Remove quarks and gluons
          isDiquark(pidabs) || // Remove diquarks
          isColorState(pidabs) // Remove particles like octet states for charmonia and similar
      ) {
        continue;
      }
      // Intermediate particles have both a production vertex and an end vertex
      if (end_vtx) {
        if (!mStoreDecayedPrimaries) {
          continue;
        }
        tobetracked = 0; // decayed particles should not be tracked by GEANT4
      } else {
        tobetracked = 1;
      }
      partID2StoreID[p->id()] = nparticlesgood++;

      int mothertobestored = -1;
      for (const auto& mother : prod_vtx->particles_in()) {
        auto mothEntry = partID2StoreID.find(mother->id());
        if (mothEntry != partID2StoreID.end()) {
          mothertobestored = mothEntry->second;
          break;
        }
      }

      phi = p->momentum().phi();
      pt = p->momentum().perp();
      pZ = p->momentum().pz();
      en = p->momentum().e();
      status = p->status();
      sn = std::sin(phi);
      cs = std::cos(phi);
      pdgCode = p->pid();
      HepMC3::FourVector pos = prod_vtx->position();
      vx = pos.x() + xv; // shift to account for generated position of the primary vertex
      vy = pos.y() + yv;
      vz = pos.z() + zv;
      tof = pos.t();
      int dummy = 0;

      getStack()->PushTrack(tobetracked, mothertobestored, pdgCode, pt * cs, pt * sn, pZ, en, vx, vy, vz, tof, 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., status);
    } else {
      LOGP(info, "particle ID# {} with no production vertex\n", p->id());
      continue;
    }
  }

  auto info = fmt::format("{}_x{}", getName(), nparticlesgood);                              // test
  mcHead->getGenHeaders().emplace_back(nparticlesgood, 0, mcHead->getNPrimaries(), 0, info); // test
  mcHead->incNPrimaries(nparticlesgood);
  ++mReadEvents;
}
