#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#include "HepMC3/WriterRootTree.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TFile.h"

#include <filesystem>

bool isDiquark(int);

using namespace Pythia8;
namespace fs = std::filesystem;

int main(int argc, char* argv[])
{

  // Default for external parameters
  int nEvent = 100;                   // number of events
  float energy = 150.;                // projectile energy per nucleon
  float ymin = -0.5, ymax = 0.5;      // rapidity range for output
  int specparticle = 213;             // Additional histos for a specific particle outside the "common" ones
  std::string inpFile = "pp_MB.cmnd"; // external command file

  // Reading external parameters
  bool help = false;
  for (int ig = 0; ig < argc; ++ig) {
    TString str = argv[ig];
    if (str.Contains("--help") || str.Contains("--help")) {
      help = true;
    } else if (str.Contains("--events")) {
      sscanf(argv[++ig], "%d", &nEvent);
    } else if (str.Contains("--energy")) {
      sscanf(argv[++ig], "%f", &energy);
    } else if (str.Contains("--ymin")) {
      sscanf(argv[++ig], "%f", &ymin);
    } else if (str.Contains("--ymax")) {
      sscanf(argv[++ig], "%f", &ymax);
    } else if (str.Contains("--particle")) {
      sscanf(argv[++ig], "%d", &specparticle);
    } else if (str.Contains("--cmnd")) {
      inpFile = argv[++ig];
    }
  }
  if (help) {
    printf("Available options:\n");
  }
  printf(
    "--events [%d]\n"
    "--energy [%.1f]\n"
    "--ymin   [%.2f]\n"
    "--ymax   [%.2f]\n"
    "--particle [%d]\n"
    "--cmnd [%s]\n",
    nEvent, energy, ymin, ymax, specparticle, inpFile.c_str());
  if (help) {
    return 0;
  }

  Pythia pythia;
  if (!fs::exists(inpFile)) {
    std::cerr << "Error: command file '" << inpFile << "' not found.\n";
    return 1;
  }

  pythia.readFile(inpFile);

  char en[30];
  pythia.readString(en);

  // General histograms
  TH1D* hidall = new TH1D("hidall", "particle id all", 10000, 0., 10000.);
  TH1D* hmall = new TH1D("hmall", "mass all", 100, 0., 10.);
  TH1D* hyall = new TH1D("hyall", "rapidity all", 100, -2.5, 2.5);
  TH1D* hptall = new TH1D("hptall", "transverse momentum all", 100, 0., 10.);
  TH2D* hptvsyall = new TH2D("hptvsyall", "transverse momentum vs rapidity all", 100, -2.5, 2.5, 100, 0., 10.);

  // Histograms for specific particle indicated in the input list
  TH1D* hmspec = new TH1D(Form("mass id%d", specparticle), Form("mass id%d", specparticle), 100, 0., 10.);
  TH1D* hyspec = new TH1D(Form("rapidity id%d", specparticle), Form("rapidity id%d", specparticle), 100, -2.5, 2.5);
  TH1D* hptspec = new TH1D(Form("transverse momentum_id%d", specparticle), Form("transverse momentum_id%d", specparticle), 100, 0., 10.);
  TH2D* hptvsyspec = new TH2D(Form("transverse momentum vs rapidity id%d", specparticle), Form("transverse momentum vs rapidity id%d", specparticle), 100, -2.5, 2.5, 100, 0., 10.);

  // List of "common" hadrons for tuning purposes
  const std::vector<int> hadid = {211, 321, 2212, 3122, 3222, 3112, 3322, 3312, 3334, 130, 310}; // pi, kch, p, lambda, sigma (2 states), xi (2 states), omega, ksh, kl
  // and related histograms
  std::vector<TH1I> hmult;
  for (int i = 0; i < (int)hadid.size(); i++) {
    hmult.emplace_back(Form("hmult_id%d", hadid[i]), Form("hmult_id%d", hadid[i]), 1000, 0., 1000.);
  }
  std::vector<TH1D> hpt;
  std::vector<TH1D> hmt;
  std::vector<TH1D> hy;
  for (int i = 0; i < (int)hadid.size(); i++) {
    hpt.emplace_back(Form("hpt_id%d", hadid[i]), Form("hpt_id%d", hadid[i]), 100, 0., 10.);
    hmt.emplace_back(Form("hmt_id%d", hadid[i]), Form("hmt_id%d", hadid[i]), 100, 0., 10.);
    hy.emplace_back(Form("hy_id%d", hadid[i]), Form("hy_id%d", hadid[i]), 100, -2.5, 2.5);
  }

  std::vector<int> partons = {1, 2, 3, 4, 5, 6, 21}; // to remove quarks and gluons

  pythia.init();

  float ycm = 0.5 * 0.5 * TMath::Log((energy + TMath::Sqrt(energy * energy - 0.938 * 0.938)) / (energy - TMath::Sqrt(energy * energy - 0.938 * 0.938)));

  // Initialize HepMC3 ROOT output
  HepMC3::WriterRootTree writer("GenPYTHIA_HepMC3.root");
  HepMC3::Pythia8ToHepMC3 toHepMC3; // converter object

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    if (!pythia.next())
      continue;

    // print heavy-ion info
    const Pythia8::HIInfo* hi = pythia.info.hiInfo;
    if (hi) {
      printf("Impact parameter = %3.2f Ncoll = %d Npart = %d\n", hi->b(), hi->nCollTot(), hi->nPartProj() + hi->nPartTarg());
    }

    // Write event to HepMC3 tree
    HepMC3::GenEvent hepmcevt;
    toHepMC3.fill_next_event(pythia, &hepmcevt); // convert Pythia to HepMC3
    writer.write_event(hepmcevt);                // write directly to ROOT

    std::vector<int> hadmult((int)hadid.size(), 0);

    // Loop for event analysis
    for (int i = 0; i < pythia.event.size(); ++i) {

      if (std::find(partons.begin(), partons.end(), std::abs(pythia.event[i].id())) != partons.end())
        continue; // it is a parton, skip
      if (isDiquark(pythia.event[i].id()))
        continue; // it is a diquark, skip
      if (pythia.event[i].id() == 90)
        continue; // it is the initial system, skip
      if (!pythia.event[i].isFinal())
        continue; // skip non-final particles

      TLorentzVector cckin;
      cckin.SetPxPyPzE(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
      float mcc = cckin.M();
      float ycc = cckin.Rapidity();
      float ptcc = cckin.Pt();
      float mtcc = std::sqrt(ptcc * ptcc + mcc * mcc);
      bool ysel = (ycc - ycm) > ymin && (ycc - ycm) < ymax;

      // Extract y, pt, mt and multiplicity for "standard" hadrons and fill histograms
      for (int j = 0; j < (int)hadid.size(); j++) {
        auto it = std::find(hadid.begin(), hadid.end(), std::abs(pythia.event[i].id()));
        if (it != hadid.end()) {
          size_t index = std::distance(hadid.begin(), it);
          hy[index].Fill(ycc - ycm);
          if (ysel) {
            hadmult[index]++;
            hpt[index].Fill(ptcc);
            hmt[index].Fill(mtcc - mcc);
          }
          break;
        }
      }

      // Fill general histograms
      hidall->Fill(pythia.event[i].id());
      hmall->Fill(mcc);
      hyall->Fill(ycc - ycm);
      hptall->Fill(ptcc);
      hptvsyall->Fill(ycc - ycm, ptcc);
      // Fill histograms for a specific particle
      if (std::abs(pythia.event[i].id()) == specparticle) {
        hmspec->Fill(mcc);
        hyspec->Fill(ycc - ycm);
        hptspec->Fill(ptcc);
        hptvsyspec->Fill(ycc - ycm, ptcc);
      }

      //      printf("Event %d particle %d rapidity %f id %d name %s ysel =%d\n",iEvent,i,ycc-ycm,pythia.event[i].id(),pythia.event[i].name().c_str(),ysel)
    }

    for (int k = 0; k < (int)hadid.size(); k++) {
      hmult[k].Fill(hadmult[k]);
    }
  }
  writer.close();
  pythia.stat();

  // Combine infos for various sigma and Xi baryons
  TH1D* hptsigmaplus = &hpt[4];
  TH1D* hptsigmaminus = &hpt[5];
  TH1D* hptxi0 = &hpt[6];
  TH1D* hptximinus = &hpt[7];
  hptsigmaplus->Add(hptsigmaminus);
  hptxi0->Add(hptximinus);

  // Calculate Tslope for main hadrons
  std::vector<float> Tslope, errTslope;
  TF1* fmt = new TF1("fmt", "[0]*pow(x,[2])*exp(-x/[1])", 0., 10.);
  for (int nh = 0; nh < (int)hadid.size(); nh++) {
    fmt->SetParameters(1000., 0.2, 1.);
    hmt[nh].Fit(fmt);
    Tslope.push_back(fmt->GetParameter(1));
    errTslope.push_back(fmt->GetParError(1));
  }

  char fname[40];
  sprintf(fname, "PYTHIA8_controlhistos.root");
  TFile* fcc = new TFile(fname, "RECREATE");
  hidall->Write();
  hmall->Write();
  hyall->Write();
  hptall->Write();
  hptvsyall->Write();
  hmspec->Write();
  hyspec->Write();
  hptspec->Write();
  hptvsyspec->Write();

  for (int nh = 0; nh < (int)hadid.size(); nh++) {
    hmult[nh].Write();
    hpt[nh].Write();
    hmt[nh].Write();
    hy[nh].Write();
  }

  fcc->Close();

  // Summary printout
  float deltay = ymax - ymin;
  printf("Mean pT in %3.2f<y<%3.2f\n pi = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f, T = %4.3f +/- %4.3f\n kch = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f, T = %4.3f +/- %4.3f \n p = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f, T = %4.3f +/- %4.3f \n lambda = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f, T = %4.3f +/- %4.3f \n sigma = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f \n xi = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f \n omega = %3.2f +/- %3.2f GeV/c, dN/dy = %3.2f \n\n",
         ymin, ymax,
         hpt[0].GetMean(), hpt[0].GetMeanError(), hmult[0].GetMean() / deltay, Tslope[0], errTslope[0],
         hpt[1].GetMean(), hpt[1].GetMeanError(), hmult[1].GetMean() / deltay, Tslope[1], errTslope[1],
         hpt[2].GetMean(), hpt[2].GetMeanError(), hmult[2].GetMean() / deltay, Tslope[2], errTslope[2],
         hpt[3].GetMean(), hpt[3].GetMeanError(), hmult[3].GetMean() / deltay, Tslope[3], errTslope[3],
         hptsigmaplus->GetMean(), hptsigmaplus->GetMeanError(), (hmult[4].GetMean() + hmult[5].GetMean()) / deltay,
         hptxi0->GetMean(), hptxi0->GetMeanError(), (hmult[6].GetMean() + hmult[7].GetMean()) / deltay,
         hpt[8].GetMean(), hpt[8].GetMeanError(), hmult[8].GetMean() / deltay);
}

bool isDiquark(int ipdg)
{
  return (ipdg > 1000 && ipdg < 10000 && (ipdg / 10) % 10 == 0);
}
