// NA6PCCopyright

/*
Based on:
//------------------------------------------------
// The Virtual Monte Carlo examples
// Copyright (C) 2007 - 2015 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file  E03/include/Ex03MCStack.h
/// \brief Definition of the Ex03MCStack class
///
/// Geant4 ExampleN03 adapted to Virtual Monte Carlo
///
/// \author I. Hrivnacova; IPN, Orsay
*/

#ifndef NA6P_MCSTACK_H
#define NA6P_MCSTACK_H

#include <TVirtualMCStack.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <stack>

class TParticle;
class TClonesArray;

class NA6PMCStack : public TVirtualMCStack
{
 public:
  NA6PMCStack(int size);
  NA6PMCStack();
  ~NA6PMCStack() override;

  // methods
  void PushTrack(int toBeDone, int parent, int pdg, double px, double py, double pz, double e,
                 double vx, double vy, double vz, double tof, double polx, double poly, double polz,
                 TMCProcess mech, int& ntr, double weight, int is) override;
  TParticle* PopNextTrack(int& track) override;
  TParticle* PopPrimaryForTracking(int i) override;
  void Print(Option_t* option = "") const override;
  void clear();

  // set methods
  void SetCurrentTrack(int track) override;

  // get methods
  int GetNtrack() const override { return mParticles->GetEntriesFast(); }
  int GetNprimary() const override { return mNPrimary; }
  TParticle* GetCurrentTrack() const override { return GetParticle(mCurrentTrack); }
  int GetCurrentTrackNumber() const override { return mCurrentTrack; }
  int GetCurrentParentTrackNumber() const override
  {
    TParticle* current = GetCurrentTrack();
    return current ? current->GetFirstMother() : -1;
  }
  TParticle* GetParticle(int id) const { return (TParticle*)mParticles->At(id); }

  void addHit(int detBit) { GetCurrentTrack()->SetBit(detBit); }
  void setVerbosity(int v) { mVerbosity = v; }
  auto getVerbosity() const { return mVerbosity; }

 private:
  // data members
  std::stack<TParticle*> mStack{};    //!< The stack of particles (transient)
  TClonesArray* mParticles = nullptr; ///< The array of particle (persistent)
  int mCurrentTrack = -1;             ///< The current track number
  int mNPrimary = 0;                  ///< The number of primaries
  int mVerbosity = 0;                 //!

  ClassDefOverride(NA6PMCStack, 1); // NA6PMCStack
};

#endif
