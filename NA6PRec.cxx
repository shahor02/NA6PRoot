// NA6PCCopyright

#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include <TApplication.h>
#include <TMCManager.h>
#include <TGeoGlobalMagField.h>
#include <boost/program_options.hpp>
#include <filesystem>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PVerTelHit.h"
#include "MagneticField.h"
#include "StringUtils.h"
#include "NA6PVerTelReconstruction.h"

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
    add_option("verbosity,v", bpo::value<int>()->default_value(0), "verbosity level [0 = no output]");
    add_option("configKeyValues", bpo::value<std::string>()->default_value(""), "comma-separated configKeyValues");
    add_option("load-ini", bpo::value<std::string>()->default_value(""), "load configurables from ini file (if defined), overridden by configKeyValues");
    add_option("load-recoparam", bpo::value<std::string>()->default_value(""), "load reco parameters from ini file (if defined), overridden by configKeyValues");
    add_option("firstevent,f", bpo::value<int32_t>()->default_value(0), "first event");
    add_option("lastevent,l", bpo::value<int32_t>()->default_value(-1), "last event");
    add_option("doHitsToRecPoints,cl", bpo::value<bool>()->default_value(true), "run hits->clusters");
    add_option("doTrackletVertex,vert", bpo::value<bool>()->default_value(true), "run tracklet vertexer");
    add_option("doVTTracking,vt", bpo::value<bool>()->default_value(true), "run VT tracker");
    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);

    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl
              << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }
  auto flini = vm["load-ini"].as<std::string>();
  na6p::conf::ConfigurableParam::updateFromString(vm["configKeyValues"].as<std::string>());
  if (!flini.empty()) {
    na6p::conf::ConfigurableParam::updateFromFile(flini, "", true);
  }
  auto flrp = vm["load-recoparam"].as<std::string>();
  if (!flrp.empty()) {
    na6p::conf::ConfigurableParam::updateFromFile(flrp, "", true);
  }
  LOGP(info, "Printing all configs");
  na6p::conf::ConfigurableParam::printAllKeyValuePairs();
  int firstEv = vm["firstevent"].as<int32_t>();
  int lastEv = vm["lastevent"].as<int32_t>();

  // mag field definition
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  NA6PVerTelReconstruction* vtrec = new NA6PVerTelReconstruction();
  vtrec->setGeometryFile("geometry.root");

  if (vm["doHitsToRecPoints"].as<bool>()) {

    TFile* fhVT = TFile::Open("HitsVerTel.root");
    if (!fhVT || fhVT->IsZombie()) {
      LOGP(error, "Cannot open HitsVerTel.root");
      if (fhVT)
        delete fhVT;
      return -1;
    }
    TTree* thVT = (TTree*)fhVT->Get("hitsVerTel");
    if (!thVT) {
      LOGP(error, "Cannot find tree 'hitsVerTel' in HitsVerTel.root");
      fhVT->Close();
      delete fhVT;
      return -1;
    }
    std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
    thVT->SetBranchAddress("VerTel", &vtHitsPtr);
    int nEvVT = thVT->GetEntriesFast();

    vtrec->createClustersOutput();
    for (int jEv = 0; jEv < nEvVT; jEv++) {
      thVT->GetEvent(jEv);
      int nHits = vtHits.size();
      LOGP(info, "VerTel Event {} nHits= {}", jEv, nHits);
      vtrec->clearClusters();
      vtrec->hitsToRecPoints(vtHits);
      vtrec->writeClusters();
    }
    vtrec->closeClustersOutput();
    fhVT->Close();
    delete fhVT;
  } else {
    LOGP(info, "Hits -> Recpoints disabled from input options");
  }

  if (vm["doVTTracking"].as<bool>() || vm["doTrackletVertex"].as<bool>()) {

    TFile* fk = TFile::Open("MCKine.root");
    if (!fk || fk->IsZombie()) {
      LOGP(error, "Cannot open MCKine.root");
      if (fk)
        delete fk;
      return -1;
    }
    TTree* mcTree = (TTree*)fk->Get("mckine");
    if (!mcTree) {
      LOGP(error, "Cannot find tree 'mckine' in MCKine.root");
      fk->Close();
      delete fk;
      return -1;
    }
    std::vector<TParticle>* mcArr = nullptr;
    mcTree->SetBranchAddress("tracks", &mcArr);

    TFile* fc = TFile::Open("ClustersVerTel.root");
    if (!fc || fc->IsZombie()) {
      LOGP(error, "Cannot open ClustersVerTel.root");
      if (fc)
        delete fc;
      return -1;
    }
    TTree* tc = (TTree*)fc->Get("clustersVerTel");
    if (!tc) {
      LOGP(error, "Cannot find tree 'clustersVerTel' in ClustersVerTel.root");
      fc->Close();
      delete fc;
      return -1;
    }
    std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
    tc->SetBranchAddress("VerTel", &vtClusPtr);

    int nEv = tc->GetEntries();
    if (lastEv > nEv || lastEv < 0)
      lastEv = nEv;
    if (firstEv < 0)
      firstEv = 0;

    if (vm["doTrackletVertex"].as<bool>())
      vtrec->initVertexer();
    if (vm["doVTTracking"].as<bool>())
      vtrec->initTracker();

    TStopwatch timer;
    timer.Start();

    for (int jEv = firstEv; jEv < lastEv; jEv++) {
      LOGP(info, "Process event {}", jEv);
      vtrec->clearEvent();
      mcTree->GetEvent(jEv);
      int nPart = mcArr->size();
      float zvert = 0;
      // get primary vertex position from the Kine Tree
      for (int jp = 0; jp < nPart; jp++) {
        auto curPart = mcArr->at(jp);
        if (curPart.IsPrimary()) {
          zvert = curPart.Vz();
          break;
        }
      }
      tc->GetEvent(jEv);
      NA6PVertex pvert;
      pvert.setXYZ(0., 0., zvert);
      vtrec->setPrimaryVertex(&pvert);
      vtrec->setClusters(vtClus);
      if (vm["doTrackletVertex"].as<bool>())
        vtrec->runVertexerTracklets();
      if (vm["doVTTracking"].as<bool>())
        vtrec->runTracking();
    }
    vtrec->closeVerticesOutput();
    vtrec->closeTracksOutput();
    timer.Stop();
    LOGP(info, "------ VT reconstruction time ------");
    timer.Print();
    LOGP(info, "------------------------------");
    fk->Close();
    delete fk;
    fc->Close();
    delete fc;
  }
  delete vtrec;

  return 0;
}
