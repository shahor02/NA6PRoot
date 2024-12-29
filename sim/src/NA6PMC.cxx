// NA6PCCopyright

#include "NA6PMC.h"
#include "NA6PDetector.h"
#include "NA6PModule.h"
#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
#include "StringUtils.h"
#include "MiscUtils.h"
#include "KeyValParam.h"

#include <TGeoManager.h>
#include <TVirtualMC.h>
#include <TMCManager.h>
#include <TVirtualMCStack.h>
#include <TGeoVolume.h>
#include <TMCProcess.h>
#include <TInterpreter.h>
#include <TFile.h>
#include <TTree.h>
#include <fairlogger/Logger.h>
#include <filesystem>

using Str=na6p::utils::Str;

NA6PMC::NA6PMC(const char* name, const char* title) : TVirtualMCApplication(name, title)
{
  mDet = std::make_unique<NA6PDetector>();
  mStack = std::make_unique<NA6PMCStack>(200);
  LOGP(info, "NA6PMC application initialized.");
}

NA6PMC::~NA6PMC()
{
  closeKineOutput();
  mDet->closeHitsOutput();
  LOGP(info, "NA6PMC application terminated.");
}

void NA6PMC::ConstructGeometry()
{
  LOGP(info, "Constructing geometry");
  mDet->createGeometry();
  TVirtualMC::GetMC()->SetRootGeometry();
}

void NA6PMC::InitGeometry() {
  LOGP(info, "Init geometry...");
}

void NA6PMC::ConstructOpGeometry()
{
  LOGP(info, "Constructing optional geometry...");
}

void NA6PMC::AddParticles()
{
  if (mVerbosity>0) {
    LOGP(info, "Adding particles to the simulation...");
  }
  // Add particle definitions here
}

void NA6PMC::GeneratePrimaries()
{
  // Implement primary particle generation here
  mGenerator->generate();
  if (mVerbosity>0) {
    LOGP(info, "Generated {} primary particles with {}", mStack->GetNprimary(), mGenerator->getName());
    if (mVerbosity>2) {
      mStack->Print();
    }
  }
}

void NA6PMC::BeginEvent()
{
  if (mVerbosity>0) {
    LOGP(info, "Beginning event {}", mEvCount);
  }
  mStack->clear();
  clearHits();
  mMCHeader.clear();
  mMCTracks.clear();
}

void NA6PMC::FinishEvent()
{
  if (mVerbosity>0) {    
    std::string rep = fmt::format(" Tracked {} particles ({} primaries) | Hits:", mStack->GetNtrack(), mStack->GetNprimary());
    for (int im=0; im<mDet->getNActiveModules(); im++) {
      const auto* det = mDet->getActiveModule(im);
      rep += fmt::format(" {}:{}", det->getName(), det->getNHits());
    }
    LOGP(info, "Finishing event {} | {}", mEvCount, rep);
  }
  selectTracksToSave();
  writeKine();
  mDet->writeHits(mRemap);
  mEvCount++;
}

void NA6PMC::BeginPrimary()
{
  if (mVerbosity > 1) {
    LOGP(info, "Beginning primary particle tracking");
  }
}

void NA6PMC::FinishPrimary()
{
  if (mVerbosity > 1) {
    LOGP(info, "Finishing primary particle tracking");
  }
}

void NA6PMC::Stepping()
{
  // Implement per-step processing logic here
  if (!TVirtualMC::GetMC()->TrackCharge()) {
    return;
  }
  int volCopy;
  TVirtualMC::GetMC()->CurrentVolID(volCopy);
  if (NA6PModule::isSensor(volCopy)) {
    mDet->getModule(NA6PModule::volID2ModuleID(volCopy))->stepManager(volCopy);
  } else {
    // non-sensitive volumes treatment
  }
}

bool NA6PMC::setupGenerator(const std::string& s)
{
  TInterpreter::EErrorCode errCode = TInterpreter::EErrorCode::kNoError;
  NA6PGenerator* gen = nullptr;
  if (Str::endsWith(s, ".C") || Str::endsWith(s, ".C+") || Str::endsWith(s, ".C++") || Str::endsWith(s, ".C+g") || Str::endsWith(s, ".cxx++g")) {
    TInterpreter::Instance()->LoadMacro(s.c_str(), &errCode);
    if (errCode != TInterpreter::EErrorCode::kNoError) {
      LOGP(fatal, "Loading of generator macro {} leads to error TInterpreter::EErrorCode={}", s, int(errCode));
    }
    std::filesystem::path mpth(s);
    std::string confMacroFun = fmt::format("{}()",mpth.stem().string());    
    gen = dynamic_cast<NA6PGenerator*>((TObject*)TInterpreter::Instance()->ProcessLineSynch(confMacroFun.c_str(), &errCode));
    if (errCode != TInterpreter::EErrorCode::kNoError || !gen) {
      LOGP(fatal, "Execution of generator function {} from macro {} leads to error TInterpreter::EErrorCode={} and returned pointer {}",
	   confMacroFun, s, int(errCode), (void*)gen);
    }
    gen->setStack(mStack.get());
    if (mRandomSeed) {
      gen->setRandomSeed(mRandomSeed);
    }
    gen->init();
    mGenerator.reset(gen);
    LOGP(info, "Initialized generator {} via macro {}", mGenerator->getName(), s);
  }

  return mGenerator != nullptr;
}

void NA6PMC::setRandomSeed(Long64_t r)
{
  if (r<0) {
    std::chrono::system_clock::time_point t0, t = std::chrono::high_resolution_clock::now();
    mRandomSeed = std::chrono::duration_cast<std::chrono::nanoseconds>(t-t0).count();
  } else {
    mRandomSeed = r;  
  }
  if (mRandomSeed) {
    gRandom->SetSeed(mRandomSeed);
    LOGP(info, "Set gRandom seed to {}", mRandomSeed);
  }
}

void NA6PMC::init()
{
  if (mVerbosity < 2) {
    MiscUtils::silenceStdOut("for BuildPhysics()");
  }  

  TVirtualMC::GetMC()->Init();

  TVirtualMC::GetMC()->BuildPhysics();
  if (mVerbosity < 2) {
    MiscUtils::reviveStdOut();
  }  
  
  TVirtualMC::GetMC()->SetStack(mStack.get());

  mDet->setVerbosity(mVerbosity);
  mStack->setVerbosity(mVerbosity);
  if (mGenerator) {
    mGenerator->setVerbosity(mVerbosity);
  }

  auto dir = na6p::utils::Str::rectifyDirectory(na6p::conf::KeyValParam::Instance().output_dir);
  if (dir.empty()) {
    dir = na6p::utils::Str::rectifyDirectory(".");
  }
  createKineOutput(dir);
  mDet->createHitsOutput(dir);
}
  
void NA6PMC::clearHits()
{
  if (mVerbosity>0) {
    LOGP(info, "clearHits");
  }
  for (int i=0;i<mDet->getNModules();i++) {
    auto det = mDet->getModule(i);
    det->clearHits();
  }
}

void NA6PMC::createKineOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}MCKine.root",outDir);
  mKineFile = TFile::Open(nm.c_str(), "recreate");
  mKineTree = new TTree("mckine", "MC kinematics");
  mKineTree->Branch("header", &mMCHeaderPtr);
  mKineTree->Branch("tracks", &mMCTracksPtr);
  LOGP(info, "Will store MC Kinematics in {}", nm);
}

void NA6PMC::closeKineOutput()
{
  if (mKineTree && mKineFile) {
    mKineFile->cd();
    mKineTree->Write();
    delete mKineTree;
    mKineTree = 0;
    mKineFile->Close();
    delete mKineFile;
    mKineFile = 0;
  }
}

void NA6PMC::writeKine()
{  
  if (mKineTree) {
    mKineTree->Fill();    
  }
}

void NA6PMC::selectTracksToSave()
{
  // store track if primary or if it has hits (then store all ancestors)
  int ntrIni = mStack->GetNtrack(), nPrimIni = mStack->GetNprimary();
  mRemap.clear();
  mRemap.resize(mStack->GetNtrack(), -1);
  int nkeep = 0;
  for (int i=0;i<nPrimIni;i++) { // generator primaries are stored as is
    mRemap[i] = nkeep++;
  }
  for (int i=nPrimIni;i<ntrIni;i++) {
    const auto* part = mStack->GetParticle(i);
    if (NA6PModule::testActiveIDBits(*part)) { // has hits
      mRemap[i] = 0;
      int mothID;
      LOGP(debug, "will save contributor {}, mother {}({})", i, part->GetFirstMother(), part->GetFirstMother()>=0 ? mRemap[part->GetFirstMother()] : -1);
      while ((mothID=part->GetFirstMother())>=0 && mRemap[mothID] < 0) {
	mRemap[mothID] = 0;
	part = mStack->GetParticle(mothID);
	LOGP(debug, "will save mother {}, its mother {}({})", mothID, part->GetFirstMother(), part->GetFirstMother()>=0 ? mRemap[part->GetFirstMother()] : -1);
      }
    }
  }
  // enumerate mRemapped indices
  for (int i=nPrimIni;i<ntrIni;i++) {
    if (mRemap[i] == 0) {
      mRemap[i] = nkeep++;
    }
  }
  // update headers
  mMCHeader.setNTracks(nkeep);
  mMCHeader.setNPrimaries(nPrimIni);
  mMCHeader.setEventID(mEvCount);
  
  // record selected tracks, mRemapping mother/daughter indices
  for (int i=0;i<ntrIni;i++) {
    if (mRemap[i] >= 0) {
      auto* part = mStack->GetParticle(i);
      mMCTracks.push_back(*part);      
      auto& partN = mMCTracks.back();
      partN.SetLastMother(-1);
      int mothID = partN.GetFirstMother();
      if (mothID >= 0) {
	if (mRemap[mothID]<0) {
	  LOGP(debug, "Mother {} of {} should have been spared but it is not!", mothID, i);
	}
	partN.SetFirstMother(mRemap[mothID]);
	auto& moth = mMCTracks[mRemap[mothID]]; // set daughters
	if (moth.GetFirstDaughter()<0) {
	  moth.SetFirstDaughter(mRemap[i]);
	}
	moth.SetLastDaughter(mRemap[i]);
      }
    }
  }
  
  LOGP(info, "Will save {} tracks out of {}", mMCTracks.size(), mRemap.size());
}
