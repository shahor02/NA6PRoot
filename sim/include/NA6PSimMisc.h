// NA6PCCopyright
#ifndef NA6P_SIM_MISC_H_
#define NA6P_SIM_MISC_H_

// misc definitions for simulation

struct UserHook {
  static constexpr int ADDParticles = 0;                  // 1st arg passed to hook method at the entrance and exit of NA6PMC::AddParticles
  static constexpr int SelectParticles = 1;               // 1st arg passed to hook method from the NA6PMC::selectTracksToSave
  static constexpr int KeepParticleBit = 0x1 << (9 + 14); // user may set this bit to the particle in the user hook method call with SelectParticles

  ClassDefNV(UserHook, 1);
};

#endif
