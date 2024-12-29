// NA6PCCopyright

#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include <TApplication.h>
#include <TMCManager.h>
#include <TGeoGlobalMagField.h>
#include <boost/program_options.hpp>
#include "MagneticField.h"
#include "NA6PMC.h"
#include "TG4RunConfiguration.h"
#include "TGeant4.h"

int main(int argc, char** argv)
{
  const std::string layoutIni = "na6pLayout.ini";
  fair::Logger::OnFatal([]() { throw std::runtime_error("Fatal error"); });

  namespace bpo = boost::program_options;
  bpo::variables_map vm;
  bpo::options_description opt_general("Usage:\n  " + std::string(argv[0]));
  bpo::options_description opt_hidden("");
  bpo::options_description opt_all;
  bpo::positional_options_description opt_pos;

  try {
    auto add_option = opt_general.add_options();
    add_option("help,h", "Print this help message");
    add_option("verbosity,v", bpo::value<uint32_t>()->default_value(0), "verbosity level [0 = no output]");
    add_option("configKeyValues", bpo::value<std::string>()->default_value(""), "comma-separated configKeyValues");
    add_option("load-ini", bpo::value<std::string>()->default_value(""), "load configurables from ini file (if defined), overridden by configKeyValues");
    add_option("disable-write-ini", bpo::value<bool>()->default_value(false)->implicit_value(true), "do not write ini file");
    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);

    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }
  auto flini = vm["load-ini"].as<std::string>();
  na6p::conf::ConfigurableParam::updateFromString(vm["configKeyValues"].as<std::string>());
  if (!flini.empty()) {
    na6p::conf::ConfigurableParam::updateFromFile(flini,"",true);
  }
  LOGP(info, "Printing all configs");
  na6p::conf::ConfigurableParam::printAllKeyValuePairs();

  auto mc = new NA6PMC("NA6PMCApp", "NA6P Virtual Monte Carlo Application");
  auto runConfig = new TG4RunConfiguration("geomRoot", "FTFP_BERT");
  auto geant4 = new TGeant4("TGeant4", "Geant4 Monte Carlo Engine", runConfig, argc, argv);
  auto magField = new MagneticField();
  magField->loadFlukaField();
  magField->setAsGlobalField();
  geant4->SetMagField( TGeoGlobalMagField::Instance()->GetField() );
  //mc->InitMC("");
  const int nEvents = 10; // Number of events to simulate
  for (int i = 0; i < nEvents; ++i) {
    std::cout << "Processing event " << i + 1 << "..." << std::endl;
    
    // Begin the event
    mc->BeginEvent();
    
    // Process the event using Geant4
    geant4->ProcessRun(i);
    
    // Finish the event
    mc->FinishEvent();
  }
  if (!vm["disable-write-ini"].as<bool>()) {
    na6p::conf::ConfigurableParam::writeINI(layoutIni);
    LOGP(info, "Stored configurable params to {}/{}", na6p::conf::ConfigurableParam::getOutputDir(), layoutIni);
  }
  //  ~/aliroot/sw/SOURCES/GEANT4_VMC/v5-3/v5-3/examples/Gflash
  delete geant4;
  delete mc;
  return 0;
}
