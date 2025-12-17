set(THIS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(GenPYTHIA GenPYTHIA.cxx)

# Keep it consistent with the rest of the project (do not force global flags here).
target_compile_features(GenPYTHIA PRIVATE cxx_std_17)

# Ensure CMake can find custom FindPythia8.cmake shipped in <repo>/dependencies
# (Safe even if the top-level already set CMAKE_MODULE_PATH.)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/dependencies")

# --- Dependencies ------------------------------------------------------------

# ROOT
# Prefer imported targets. Components list is conservative; adjust if your GenPYTHIA needs more.
#find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist)

# HepMC3 (usually provides a config package with imported targets)
# Prefer CONFIG; fall back to MODULE if your environment requires it.
#find_package(HepMC3 REQUIRED CONFIG)

# Pythia8 via <repo>/dependencies/FindPythia8.cmake
#find_package(Pythia8 REQUIRED)

# --- Linking ----------------------------------------------------------------
# ROOT
target_link_libraries(GenPYTHIA PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)

# HepMC3 main library
if(TARGET HepMC3::HepMC3)
  target_link_libraries(GenPYTHIA PRIVATE HepMC3::HepMC3)
else()
  # Fallback if the package exports a different target name
  message(FATAL_ERROR "HepMC3 imported target HepMC3::HepMC3 was not found. Check your HepMC3 CMake package.")
endif()

# HepMC3 ROOT IO (equivalent to -lHepMC3rootIO)
if(TARGET HepMC3::rootIO)
  target_link_libraries(GenPYTHIA PRIVATE HepMC3::rootIO)
else()
  # Fallback: link by library name if no imported target exists
  find_library(HEPMC3_ROOTIO_LIB HepMC3rootIO)
  if(HEPMC3_ROOTIO_LIB)
    target_link_libraries(GenPYTHIA PRIVATE ${HEPMC3_ROOTIO_LIB})
  else()
    message(FATAL_ERROR "HepMC3rootIO requested but neither HepMC3::rootIO nor libHepMC3rootIO was found.")
  endif()
endif()

# Pythia8 (support both imported-target and variable-based FindPythia8.cmake styles)
if(TARGET Pythia8::Pythia8)
  target_link_libraries(GenPYTHIA PRIVATE Pythia8::Pythia8)
elseif(DEFINED PYTHIA8_LIBRARIES)
  target_link_libraries(GenPYTHIA PRIVATE ${PYTHIA8_LIBRARIES})
elseif(DEFINED PYTHIA8_LIBRARY)
  target_link_libraries(GenPYTHIA PRIVATE ${PYTHIA8_LIBRARY})
else()
  message(FATAL_ERROR "Pythia8 found, but no target/variables (Pythia8::Pythia8 or PYTHIA8_LIBRARIES) were provided.")
endif()

# --- Install (so it’s available with `make install` like the rest of the repo) ---
install(TARGETS GenPYTHIA RUNTIME DESTINATION bin)
