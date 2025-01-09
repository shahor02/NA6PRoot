#include "TGeoManager.h"
#include "TVirtualMC.h"
#include "TVirtualMCApplication.h"
#include "TMCManager.h"
#include "TVirtualMCStack.h"
#include "TGeoVolume.h"
#include "TMCProcess.h"
#include "TMCGeometry.h"
#include <iostream>

// Forward declarations
void CreateGeometry(); // Assuming this is already defined to create TGeo geometry

class NA6PSim : public TVirtualMCApplication
{
 public:
  NA6PSim(const char* name, const char* title) : TVirtualMCApplication(name, title) {}

  void InitMC() override
  {
    std::cout << "Initializing Virtual Monte Carlo..." << std::endl;
    TVirtualMC::GetMC()->Init();
  }

  void ConstructGeometry() override
  {
    std::cout << "Constructing geometry..." << std::endl;
    CreateGeometry();
    gGeoManager->CloseGeometry();
    TVirtualMC::GetMC()->SetRootGeometry();
  }

  void RunMC(Int_t nEvents)
  {
    for (Int_t i = 0; i < nEvents; ++i) {
      std::cout << "Processing event: " << i + 1 << std::endl;
      TVirtualMC::GetMC()->ProcessRun(i);
    }
  }

  void FinishRun()
  {
    std::cout << "Finishing simulation run..." << std::endl;
    TVirtualMC::GetMC()->FinishRun();
  }

  void FinishEvent() override
  {
    std::cout << "Finishing event..." << std::endl;
    TVirtualMC::GetMC()->FinishEvent();
  }

  void BeginEvent() override
  {
    std::cout << "Starting new event..." << std::endl;
    TVirtualMC::GetMC()->BeginEvent();
  }

  void StepManager() override
  {
    // Implement step-specific processing if required
    std::cout << "StepManager processing step..." << std::endl;
    /*
      void StepManager() override {
    // Get the current volume name
    const char* volumeName = TVirtualMC::GetMC()->CurrentVolName();

    // Get the track ID and step information
    Int_t trackID = TVirtualMC::GetMC()->GetStack()->GetCurrentTrackNumber();
    Double_t pos[3], mom[3];
    Double_t step, energy, time;

    TVirtualMC::GetMC()->TrackPosition(pos[0], pos[1], pos[2]);
    TVirtualMC::GetMC()->TrackMomentum(mom[0], mom[1], mom[2]);
    step = TVirtualMC::GetMC()->CurrentStep();
    energy = TVirtualMC::GetMC()->Etot();
    time = TVirtualMC::GetMC()->TrackTime();

    // Perform volume-dependent action
    if (std::string(volumeName) == "Volume1") {
        std::cout << "In Volume1: Track " << trackID << " at position ("
                  << pos[0] << ", " << pos[1] << ", " << pos[2] << ")"
                  << " with step length " << step << " cm." << std::endl;
    } else if (std::string(volumeName) == "Volume2") {
        std::cout << "In Volume2: Track " << trackID << " has total energy "
                  << energy << " GeV and momentum ("
                  << mom[0] << ", " << mom[1] << ", " << mom[2] << ")." << std::endl;
    } else {
        std::cout << "In volume " << volumeName << ": No specific action defined." << std::endl;
    }
    }
     */
  }
};

int main(int argc, char** argv)
{
  // Initialize ROOT application
  TApplication app("VMCApp", &argc, argv);

  // Create Virtual Monte Carlo Manager
  TMCManager* mcManager = TMCManager::Instance();
  mcManager->SetApplication(new NA6PSim("VMCApp", "My Virtual Monte Carlo Application"));

  // Initialize Monte Carlo engine
  mcManager->Init();

  // Configure and process simulation
  mcManager->SetMC("TGeant4"); // Example: Using Geant4
  mcManager->ProcessRun(10);   // Process 10 events

  // Cleanup and finalize
  mcManager->FinishRun();

  return 0;
}
