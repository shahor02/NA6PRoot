// NA6PCCopyright

/*
Based on:
//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2014 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file E03/src/Ex03MCStack.cxx
/// \brief Implementation of the Ex03MCStack class
///
/// Geant4 ExampleN03 adapted to Virtual Monte Carlo
///
/// \date 06/03/2002
/// \author I. Hrivnacova; IPN, Orsay
*/

#include <Riostream.h>
#include <TParticle.h>
#include <fairlogger/Logger.h>
#include "NA6PMCStack.h"

using namespace std;

//_____________________________________________________________________________
NA6PMCStack::NA6PMCStack(int size) : mParticles(0), mCurrentTrack(-1), mNPrimary(0)
{
  /// Standard constructor
  /// \param size  The stack size

  mParticles = new TClonesArray("TParticle", size);
}

//_____________________________________________________________________________
NA6PMCStack::NA6PMCStack() : mParticles(0), mCurrentTrack(-1), mNPrimary(0)
{
  /// Default constructor
}

//_____________________________________________________________________________
NA6PMCStack::~NA6PMCStack()
{
  /// Destructor

  if (mParticles) {
    mParticles->Delete();
  }
  delete mParticles;
}

// private methods

// public methods

//_____________________________________________________________________________
void NA6PMCStack::PushTrack(int toBeDone, int parent, int pdg,
                            double px, double py, double pz, double e,
                            double vx, double vy, double vz, double tof,
                            double polx, double poly, double polz,
                            TMCProcess mech, int& ntr, double weight, int is)
{
  /// Create a new particle and push into stack;
  /// adds it to the particles array (mParticles) and if not done to the
  /// stack (mStack).
  /// Use TParticle::fMother[1] to store Track ID.
  /// \param toBeDone  1 if particles should go to tracking, 0 otherwise
  /// \param parent    number of the parent track, -1 if track is primary
  /// \param pdg       PDG encoding
  /// \param px        particle momentum - x component [GeV/c]
  /// \param py        particle momentum - y component [GeV/c]
  /// \param pz        particle momentum - z component [GeV/c]
  /// \param e         total energy [GeV]
  /// \param vx        position - x component [cm]
  /// \param vy        position - y component  [cm]
  /// \param vz        position - z component  [cm]
  /// \param tof       time of flight [s]
  /// \param polx      polarization - x component
  /// \param poly      polarization - y component
  /// \param polz      polarization - z component
  /// \param mech      creator process VMC code
  /// \param ntr       track number (is filled by the stack
  /// \param weight    particle weight
  /// \param is        generation status code

  const int kFirstDaughter = -1;
  const int kLastDaughter = -1;

  TClonesArray& particlesRef = *mParticles;
  int trackId = GetNtrack();
  TParticle* particle = new (particlesRef[trackId]) TParticle(pdg, is, parent, trackId,
                                                              kFirstDaughter, kLastDaughter,
                                                              px, py, pz, e, vx, vy, vz, tof);

  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(mech);

  if (parent < 0) {
    mNPrimary++;
  }

  if (toBeDone) {
    mStack.push(particle);
  }

  ntr = GetNtrack() - 1;
  if (mVerbosity > 1) {
    LOGP(info, "PushTrack({},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{},{:.2f},{:.2f},{:.2f},{},{},{},{})",
         toBeDone, parent, pdg, px, py, pz, e, vx, vy, vz, tof, polx, poly, polz, (int)mech, ntr, weight, is);
  }
}

//_____________________________________________________________________________
TParticle* NA6PMCStack::PopNextTrack(int& itrack)
{
  /// Get next particle for tracking from the stack.
  /// \return       The popped particle object
  /// \param track  The index of the popped track

  itrack = -1;
  if (mStack.empty()) {
    return nullptr;
  }

  TParticle* particle = mStack.top();
  mStack.pop();

  if (!particle) {
    return nullptr;
  }
  mCurrentTrack = particle->GetSecondMother();
  itrack = mCurrentTrack;
  if (mVerbosity > 1) {
    LOGP(info, "PopNextTrack {} | NParticles:{}, NPrimaries:{} StackSize:{}", itrack, GetNtrack(), GetNprimary(), mStack.size());
    particle->Print();
  }
  return particle;
}

//_____________________________________________________________________________
TParticle* NA6PMCStack::PopPrimaryForTracking(int i)
{
  /// Return \em i -th particle in mParticles.
  /// \return   The popped primary particle object
  /// \param i  The index of primary particle to be popped
  if (mVerbosity > 1) {
    LOGP(info, "PopPrimaryForTracking {} | NParticles:{}, NPrimaries:{} StackSize:{}", i, GetNtrack(), GetNprimary(), mStack.size());
  }
  if (i < 0 || i >= mNPrimary) {
    LOGP(fatal, "PopPrimaryForTracking: index {} out of range {}", i, mNPrimary);
  }
  return (TParticle*)mParticles->At(i);
}

//_____________________________________________________________________________
void NA6PMCStack::Print(Option_t* /*option*/) const
{
  /// Print info for all particles.
  LOGP(info, "NA6PMCStack Info: NParticles:{}, NPrimaries:{}", GetNtrack(), GetNprimary());

  for (int i = 0; i < GetNtrack(); i++)
    GetParticle(i)->Print();
}

//_____________________________________________________________________________
void NA6PMCStack::clear()
{
  /// Delete contained particles, reset particles array and stack.

  mCurrentTrack = -1;
  mNPrimary = 0;
  mParticles->Clear();
  mMCHeader.clear();
  mPVGenerated = false;
  while (!mStack.empty()) {
    mStack.pop();
  }
}

//_____________________________________________________________________________
void NA6PMCStack::SetCurrentTrack(int track)
{
  /// Set the current track number to a given value.
  /// \param  track The current track number
  if (mVerbosity > 1) {
    LOGP(info, "SetCurrentTrack {} | NParticles:{}, NPrimaries:{} StackSize:{}", track, GetNtrack(), GetNprimary(), mStack.size());
  }
  mCurrentTrack = track;
}
