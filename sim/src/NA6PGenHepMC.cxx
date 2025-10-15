// NA6PCCopyright

#include "NA6PGenHepMC.h"
#include "NA6PMCStack.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/ReaderRootTree.h"
#include <TDatabasePDG.h>
#include <TMath.h>
#include <cmath>


NA6PGenHepMC::NA6PGenHepMC(const std::string& name, const std::string& filname) : NA6PGenerator(name), mFileName(filname), mReadEvents(0)
{}



void NA6PGenHepMC::init()
{
  if (isInitDone()) {
    return;
  }
  NA6PGenerator::init();
}

void NA6PGenHepMC::generate()
{

  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  auto mcHead = getStack()->getEventHeader();
  double xv = mcHead->getVX();
  double yv = mcHead->getVY();
  double zv = mcHead->getVZ();

  HepMC3::ReaderRootTree mInput(mFileName);
  if (mInput.failed()) {
    LOGP(fatal, "No HepMC input file found {}", mFileName);
  }
  
  mInput.skip(mReadEvents);
  HepMC3::GenEvent evt;
  mInput.read_event(evt); //Read one event
  if (mInput.failed()) 
    {
      LOGP (info, " End of file reached . Exit .\n"); 
      return;
    }

  int nparticles=0;
  int nparticlesgood=0;
  int nparticleskept=0; //test
  
  int pdgCode, status;
  float phi, pt, pZ, en, sn, cs;
  float vx, vy, vz, tof; // tof in mm/c
  bool tobetracked;
  
  std::vector<int> partons = {1, 2, 3, 4, 5, 6, 21}; // to remove quarks and gluons
  std::vector<int> mothersid; // vector to store id of particles after purging kinematics from gluons, quarks, diquarks and exotica
  
  // Loop on particles
  for (const auto p : evt.particles()) {

    auto isParton = [](int pabs) { return (pabs>=1&pabs<=6) || pabs==21; };
    auto isDiquark = [](int pabs) { return pabs > 1000 && pabs < 10000 && ((pabs/10)%10) == 0; };
    auto isColorState = [](int pabs) { return pabs>1e6; };
    
    // Upstream and downstream verticese
    auto prod_vtx = p->production_vertex();
    auto end_vtx = p->end_vertex();
    auto pidabs = std::abs(p->pid());
    
    if (prod_vtx) {
      // Intermediate particles have both a production vertex and an end vertex
      if (end_vtx) {
	if ( isParton(pidabs) ||   // Remove quarks and gluons
	     isDiquark(pidabs) ||  // Remove diquarks
	     isColorState(pidabs)  // Remove particles like octet states for charmonia and similar
	     ) {
	  mothersid.push_back(-1);
	  continue;
	} else {
	  mothersid.push_back(nparticlesgood); // Add the particle to the list of potential mothers
	  nparticlesgood++;
	  tobetracked = kFALSE; // not to be tracked by GEANT4
	}
      } else { // Final state particles do not have an end vertex
	tobetracked = kTRUE; // to be tracked by GEANT4
	mothersid.push_back(nparticlesgood);
	nparticlesgood++;
      }
     
      int mothertobestored=-1;
      for (const auto& mother : prod_vtx->particles_in()) {
	auto motherpidabs = std::abs(mother->pid());
	if( !isParton(motherpidabs) && !isDiquark(motherpidabs) && !isColorState(motherpidabs)) { // RSCHECK here was std::abs(p->pid())<1e6, I changed this to !isColorState(motherpidabs)
	  mothertobestored = mothersid[mother->id()]; //find index of the mother renormalized to the list of particles that are pushed to the stack
	  // RS? break;
	}
      }
      
      phi = p->momentum().phi();
      pt  = p->momentum().perp();
      pZ  = p->momentum().pz();
      en  = p->momentum().e();
      status = p->status();
      sn = std::sin(phi); 
      cs = std::cos(phi);
      pdgCode = p->pid();
      HepMC3::FourVector pos=prod_vtx->position();
      vx = pos.x() + xv; // shift to account for generated position of the primary vertex
      vy = pos.y() + yv;
      vz = pos.z() + zv;
      tof = pos.t();           
      int dummy = 0

      // Temporary, to avoid crash
      //      if(mothertobestored == -1) { 
      //        tobetracked = 1; //test
      nparticleskept++; //test
      int toBeDone = tobetracked ? 1 : 0;
      if (mothertobestored >= 0) {}
      getStack()->PushTrack(tobetracked,mothertobestored,pdgCode, pt*cs, pt*sn, pZ, en, vx, vy, vz, tof, 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., status); 
	//      } 
//Temporary, to avoid crash
      
    } else {
      LOGP(info, "particle {} with no production vertex\n",p->id()-1);
      mothersid.push_back(-1);
      continue;
    }
    nparticles++;
  }
  // register generator header in the MCHeader
//  static std::string info = fmt::format("{}_x{}", getName(), nparticlesgood);
//  mcHead->getGenHeaders().emplace_back(nparticlesgood, 0, mcHead->getNPrimaries(), 0, info);
//  mcHead->incNPrimaries(nparticlesgood);

// Temporary, to avoid crash
  static std::string info = fmt::format("{}_x{}", getName(), nparticleskept); //test
  mcHead->getGenHeaders().emplace_back(nparticleskept, 0, mcHead->getNPrimaries(), 0, info); //test
  mcHead->incNPrimaries(nparticleskept); //test
// Temporary, to avoid crash

  evt.clear();
  mothersid.clear();
  mInput.close();
  ++mReadEvents;
}


