#include <fairlogger/Logger.h>
#include <TVirtualMC.h>
#include "NA6PMCStack.h"
#include "NA6PMC.h"
#include "NA6PSimMisc.h"

int testHooks(int arg, bool inout)
{
  if (arg == UserHook::ADDParticles) {
    LOGP(info, "Calling user hook from {} of the AddParticles method", inout ? "entry" : "exit");
    {
      LOGP(info, "Restricting J/psi decays to dimuon channel");
      int Jpsi_decay_mode[6][3] = {{0}};
      float Jpsi_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      Jpsi_decay_mode[0][0] = 13;  // mu-
      Jpsi_decay_mode[0][1] = -13; // mu+
      TVirtualMC::GetMC()->SetDecayMode(443, Jpsi_BR, Jpsi_decay_mode);
      
      LOGP(info, "Restricting Phi decays to dimuon channel");
      int Phi_decay_mode[6][3] = {{0}};
      float Phi_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      Phi_decay_mode[0][0] = 13;  // mu-
      Phi_decay_mode[0][1] = -13; // mu+
      TVirtualMC::GetMC()->SetDecayMode(333, Phi_BR, Phi_decay_mode);

      LOGP(info, "Restricting omega2body decays to dimuon channel");
      int Omega_decay_mode[6][3] = {{0}};
      float Omega_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
      Omega_decay_mode[0][0] = 13;  // mu-
      Omega_decay_mode[0][1] = -13; // mu+
      TVirtualMC::GetMC()->SetDecayMode(223, Omega_BR, Omega_decay_mode);
    }

  } else if (arg == UserHook::SelectParticles) {
    if (!inout) {
      return 0; // do nothing
    }
    LOGP(info, "Calling user hook from selectTracksToSave method entrance");
    auto mc = (NA6PMC*)TVirtualMC::GetMC();
    auto stack = mc->getMCStack();
    int ntr = stack->GetNtrack(), nPrim = stack->GetNprimary();
    for (int i = nPrim; i < ntr; i++) {
      auto* part = stack->GetParticle(i);
      int idMoth = part->GetFirstMother();
      int mothPdg = -1;
      while (idMoth >= 0) {
        auto* currMoth = stack->GetParticle(idMoth);
        int absPdg = std::abs(currMoth->GetPdgCode());
        if (absPdg == 11 || absPdg == 22 || absPdg == 211 || absPdg == 130 || absPdg == 321 || absPdg == 2212 || absPdg == 2112) { // stop if a "stable" particle is found in the ancestors
          break;
        }
        if ((absPdg > 400 && absPdg < 600) || (absPdg > 4000 && absPdg < 6000) || absPdg == 100443 || absPdg == 20443) {
          part->SetBit(UserHook::KeepParticleBit); // force saving particle i (and its ancestors)
          break;
        }
        idMoth = currMoth->GetFirstMother();
      }
    }
  } else {
    LOGP(error, "Unknown hook ID {}", arg);
    return -1;
  }
  return 0;
}
