#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoNavigator.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TRandom3.h>
#include "NA6PVerTelHit.h"
#include "NA6PTrack.h"
#include "NA6PBaseCluster.h"
#include "NA6PFastTrackFitter.h"
#include "MagneticField.h"
#include "NA6PTreeStreamRedirector.h"
#include "ConfigurableParam.h"
#include "NA6PLayoutParam.h"
#endif


// simple macro with an example of how to use NA6PFastFitter to fit tracks starting from their hits in the VT
// uses the MC truth information to group the hits produced by the same particle

void SuperposHistos(TH1* h1, TH1* h2, TH1* h3){
  h1->SetLineColor(1);
  h1->SetLineWidth(2);
  if(h2->GetMaximum()>h1->GetMaximum()) h1->SetMaximum(1.02*h2->GetMaximum());
  if(h3->GetMaximum()>h1->GetMaximum()) h1->SetMaximum(1.02*h3->GetMaximum());
  h1->SetMinimum(0);
  h1->Draw();
  h2->SetLineColor(2);
  h2->SetLineWidth(2);
  h2->Draw("sames");
  h3->SetLineColor(kGreen+1);
  h3->SetLineWidth(2);
  h3->Draw("sames");
  gPad->Update();
  TPaveStats* st1=(TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  if(st1){
    st1->SetY1NDC(0.72);
    st1->SetY2NDC(0.92);
    st1->SetTextColor(1);
  }
  TPaveStats* st2=(TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  if(st2){
    st2->SetY1NDC(0.51);
    st2->SetY2NDC(0.71);
    st2->SetTextColor(2);
  }
  TPaveStats* st3=(TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
  if(st3){
    st3->SetY1NDC(0.30);
    st3->SetY2NDC(0.50);
    st3->SetTextColor(kGreen+1);
  }
  gPad->Modified();
}


void FillMeanAndRms(TH2F* hImpParVsP, TH1F* hImpParMean, TH1F* hImpParRms, TH1F* hImpParSig){
  if (!hImpParVsP || !hImpParMean || !hImpParRms || !hImpParSig) {
    printf("One of pointers is not set: hImpParVsP=%p hImpParMean=%p hImpParRms=%p hImpParSig=%p\n", hImpParVsP, hImpParMean, hImpParRms, hImpParSig);
    return;
  }
  int nMomBins = hImpParVsP->GetXaxis()->GetNbins();
  hImpParMean->SetStats(0);
  hImpParRms->SetStats(0);
  hImpParSig->SetStats(0);
  for(int jp = 1; jp <= nMomBins; jp++){
    TH1D* htmp1 = hImpParVsP->ProjectionY("htmp1",jp,jp);
    double mean = htmp1->GetMean();
    double emean = htmp1->GetMeanError();
    double rms = htmp1->GetRMS();
    double erms = htmp1->GetRMSError();
    hImpParMean->SetBinContent(jp,mean);
    hImpParMean->SetBinError(jp,emean);
    hImpParRms->SetBinContent(jp,rms);
    hImpParRms->SetBinError(jp,erms);
    double gw = rms;
    double egw = erms;
    if(htmp1->GetEntries()>20. && htmp1->Fit("gaus","Q","",-3.*rms,3.*rms)) {
      TF1* fg = (TF1*)htmp1->GetListOfFunctions()->FindObject("gaus");
      if (!fg) continue;
      gw = fg->GetParameter(2);
      egw = fg->GetParError(2);
    }
    hImpParSig->SetBinContent(jp,gw);
    hImpParSig->SetBinError(jp,egw);
    delete htmp1;
  }
}

void runFastFitOnHits(int firstEv = 0,
		      int lastEv = 99999999,
		      double clures = 5.e-4,
		      const char *dirSimu = "../testN6Proot/pions/50evwithtarget",
		      int minITShits=5) {
  
  // na6p::conf::ConfigurableParam::updateFromFile(Form("%s/na6pLayout.ini",dirSimu),"",true);
  // const auto& param = NA6PLayoutParam::Instance();
  
  TFile* fk=new TFile(Form("%s/MCKine.root",dirSimu));
  TTree* mcTree=(TTree*)fk->Get("mckine");
  int nEv=mcTree->GetEntries();
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh=new TFile(Form("%s/HitsVerTel.root",dirSimu));
  TTree* th=(TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  int nMomBins = 20;
  double maxP = 10.;
  TH1F* hEtaGen = new TH1F("hEtaGen",";#eta;counts",20,0.,5.);
  TH1F* hEtaReco = new TH1F("hEtaReco",";#eta;counts",20,0.,5.);
  TH1F* hDeltaPx = new TH1F("hDeltaPx",";p_{x}^{rec}-p_{x}^{gen} (GeV/c);counts",100,-0.2,0.2);
  TH1F* hDeltaPy = new TH1F("hDeltaPy",";p_{y}^{rec}-p_{y}^{gen} (GeV/c);counts",100,-0.2,0.2);
  TH1F* hDeltaPz = new TH1F("hDeltaPz",";p_{z}^{rec}-p_{z}^{gen} (GeV/c);counts",100,-0.5,0.5);
  TH1F* hDeltaP = new TH1F("hDeltaP",";p^{rec}-p^{gen} (GeV/c);counts",100,-0.5,0.5);
  TH1F* hDeltaPhi = new TH1F("hDeltaPhi",";#varphi^{rec}-#varphi^{gen};counts",100,-M_PI/4.,M_PI/4.);
  TH1F* hDeltaEta = new TH1F("hDeltaEta",";#eta^{rec}-#eta^{gen};counts",100,-0.5,0.5);
  TH1F* hImpParX = new TH1F("hImpParX",";Track Imp. Par. X (#mum)};counts",100,-500,500);
  TH1F* hImpParY = new TH1F("hImpParY",";Track Imp. Par. Y (#mum)};counts",100,-500,500);
  TH2F* hImpParXVsP = new TH2F("hImpParXVsP",";p (GeV/c);Track Imp. Par. X (#mum)};counts",nMomBins,0.,maxP,100,-500.,500.);
  TH2F* hImpParYVsP = new TH2F("hImpParYVsP",";p (GeV/c);Track Imp. Par. Y (#mum)};counts",nMomBins,0.,maxP,100,-500.,500.);
  TH2F* hDeltaPxVsP = new TH2F("hDeltaPxVsP",";p (GeV/c);p_{x}^{rec}-p_{x}^{gen} (GeV/c);counts",nMomBins,0.,maxP,500,-0.5,0.5);
  TH2F* hDeltaPyVsP = new TH2F("hDeltaPyVsP",";p (GeV/c);p_{y}^{rec}-p_{y}^{gen} (GeV/c);counts",nMomBins,0.,maxP,500,-0.5,0.5);
  TH2F* hDeltaPzVsP = new TH2F("hDeltaPzVsP",";p (GeV/c);p_{z}^{rec}-p_{z}^{gen} (GeV/c);counts",nMomBins,0.,maxP,500,-0.5,0.5);
  TH2F* hDeltaPVsP = new TH2F("hDeltaPVsP",";p (GeV/c);p^{rec}-p^{gen} (GeV/c);counts",nMomBins,0.,maxP,100,-0.5,0.5);
  TH2F* hRelDeltaPxVsP = new TH2F("hRelDeltaPxVsP",";p (GeV/c);(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen};counts",nMomBins,0.,maxP,100,-0.5,0.5);
  TH2F* hRelDeltaPyVsP = new TH2F("hRelDeltaPyVsP",";p (GeV/c);(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen};counts",nMomBins,0.,maxP,100,-0.5,0.5);
  TH2F* hRelDeltaPzVsP = new TH2F("hRelDeltaPzVsP",";p (GeV/c);(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};counts",nMomBins,0.,maxP,100,-0.5,0.5);
  TH2F* hRelDeltaPVsP = new TH2F("hRelDeltaPVsP",";p (GeV/c);(p^{rec}-p^{gen})/p_{gen};counts",nMomBins,0.,maxP,100,-0.5,0.5);
  TH1F* hIsGood = new TH1F("hIsGood","",2,-0.5,1.5);
  hIsGood->GetXaxis()->SetBinLabel(1,"bad");
  hIsGood->GetXaxis()->SetBinLabel(2,"good");
  TH1F* hNclu = new TH1F("hNclu",";n_{ITSclus};counts",7,-0.5,6.5);
  TH1F* hChi2NDF = new TH1F("hChi2NDF",";#chi^{2}/nDOF;counts",100,0.,10.);

  if(lastEv>nEv || lastEv<0) lastEv=nEv;
  if(firstEv<0) firstEv=0;
  
  NA6PTreeStreamRedirector outTr("fitCheck.root");

  NA6PFastTrackFitter* fitter = new NA6PFastTrackFitter();
  fitter->loadGeometry(Form("%s/geometry.root",dirSimu));
  fitter->setNLayersVT(5);
  fitter->setMaxChi2Cl(10.);
  fitter->setPropagateToPrimaryVertex(true);
  // fitter->setSeedFromTwoOutermostHits();
  // fitter->setCharge(2);
  // fitter->setParticleHypothesis(2212);
  // fitter->disableMaterialCorrections();
  for(int jEv=firstEv;jEv<lastEv; jEv++){
    mcTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart=mcArr->size();
    int nHits=vtHits.size();
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for(int jp=0; jp<nPart; jp++){
      auto curPart=mcArr->at(jp);
      if(curPart.IsPrimary()){
	zvert =  curPart.Vz();
	break;
      }
    }
    fitter->setPrimaryVertexZ(zvert);
    for(int jp=0; jp<nPart; jp++){
      const auto& curPart=mcArr->at(jp);
      double pxPart=curPart.Px();
      double pyPart=curPart.Py();
      double pzPart=curPart.Pz();
      double momPart=curPart.P();
      double phiPart=curPart.Phi();
      double thetaPart=std::acos(pzPart/momPart);
      double etaPart=-std::log(std::tan(thetaPart/2.));
      std::array<std::unique_ptr<NA6PBaseCluster>, 5> clusters{ nullptr, nullptr, nullptr, nullptr, nullptr };
      double xclu[5],yclu[5],zclu[5];
      if(curPart.IsPrimary()){
	hEtaGen->Fill(etaPart);
	int maskHits = 0;
	fitter->cleanupAndStartFit();
	for(const auto& hit : vtHits) {
	  int idPart=hit.getTrackID();
	  if(idPart==jp){
	    int nDet=hit.getDetectorID();
	    int nLay=nDet/4;
	    maskHits += (1<<nLay);
	    if(clures > 0){
	      xclu[nLay] = gRandom->Gaus(hit.getX(),clures); 
	      yclu[nLay] = gRandom->Gaus(hit.getY(),clures); 
	      zclu[nLay] = hit.getZ(); 
	    }else{
	      xclu[nLay] = hit.getX(); 
	      yclu[nLay] = hit.getY(); 
	      zclu[nLay] = hit.getZ(); 
	    }
	    NA6PBaseCluster* clu = new NA6PBaseCluster(xclu[nLay], yclu[nLay], zclu[nLay], idPart);
	    clu->setDetectorID(nDet);
	    if(clures>0)  clu->setErr(clures*clures,0.,clures*clures);
	    else clu->setErr(1.e-6,0.,1e-6);
	    fitter->addClusterVT(nLay,clu);
	  }
	}
	// track only particles with 5 hits in the VT
	if (maskHits != 31) continue;

	
	// Uncomment the next lines to use the MC truth as seed for the track
	//	double xyz0[3]={xclu[4],yclu[4],zclu[4]};	
	//	double pxyz0[3]={curPart.Px(),curPart.Py(),curPart.Pz()};
	//	printf("Initialize seed at %f %f %f p = %f %f %f\n",xyz0[0],xyz0[1],xyz0[2],pxyz0[0],pxyz0[1],pxyz0[2]);
	//	int sign = curPart.GetPdgCode() > 0 ? 1 : -1;
	//	fitter->setSeed(xyz0,pxyz0,sign);

	NA6PTrack* currTr = fitter->fitTrackPointsVT();
	//	fitter->propagateToZ(currTr,zvert); 
	double chiFit = -1.;
	if(currTr){
	  hIsGood->Fill(1);
	  int nClusters = currTr->getNVTHits();
	  hNclu->Fill(nClusters);
	  if(nClusters<5) continue;
	  chiFit = currTr->getChi2();
	  //	  printf("chi2 = %f  chi2/ndf = %f\n",chiFit,currTr->GetNormChi2());
	  double pxyz[3];
	  currTr->getPXYZ(pxyz);
	  double pxtr=pxyz[0];
	  double pytr=pxyz[1];
	  double pztr=pxyz[2];
	  double momtr=currTr->getP();
	  double phitr=std::atan2(pytr,pxtr);
	  double thetatr=std::acos(pztr/momtr);
	  double etatr=-std::log(std::tan(thetatr/2.));
	  double impparX = currTr->getXLab();
	  double impparY = currTr->getYLab();
	  hEtaReco->Fill(etatr);
	  hDeltaPx->Fill(pxtr-pxPart);
	  hDeltaPy->Fill(pytr-pyPart);
	  hDeltaPz->Fill(pztr-pzPart);
	  hDeltaP->Fill(momtr-momPart);
	  hDeltaPxVsP->Fill(momtr,pxtr-pxPart);
	  hDeltaPyVsP->Fill(momtr,pytr-pyPart);
	  hDeltaPzVsP->Fill(momtr,pztr-pzPart);
	  hDeltaPVsP->Fill(momtr,momtr-momPart);
	  hRelDeltaPxVsP->Fill(momtr,(pxtr-pxPart)/pxPart);
	  hRelDeltaPyVsP->Fill(momtr,(pytr-pyPart)/pyPart);
	  hRelDeltaPzVsP->Fill(momtr,(pztr-pzPart)/pzPart);
	  hRelDeltaPVsP->Fill(momtr,(momtr-momPart)/momPart);
	  hDeltaPhi->Fill(phitr-phiPart);
	  hDeltaEta->Fill(etatr-etaPart);
	  hImpParX->Fill(impparX*1e4);
	  hImpParY->Fill(impparY*1e4);
	  hImpParXVsP->Fill(momtr,impparX*1e4);
	  hImpParYVsP->Fill(momtr,impparY*1e4);
	  hChi2NDF->Fill(chiFit);
	}else{
	  hIsGood->Fill(0); 
	}	
	TParticle trFit;
	if (currTr) {
	  double pos[3],mom[3];
	  currTr->getPXYZ(mom);
	  TLorentzVector v;
	  v.SetXYZM(mom[0], mom[1], mom[2], 0.14);
	  trFit.SetMomentum(v);
	  trFit.SetProductionVertex(currTr->getXLab(), currTr->getYLab(), currTr->getZLab(), 0);
	}
	outTr << "out" << "mcTr=" << curPart << "fitTr=" << trFit << "fitChi2=" << chiFit << "\n";
      }
    }
  }
  outTr.Close();

  TH1F* hImpParXMean = new TH1F("hImpParXMean",";p (GeV/c);<Imp Par X> (#mum)",nMomBins,0.,maxP);
  TH1F* hImpParYMean = new TH1F("hImpParYMean",";p (GeV/c);<Imp Par Y> (#mum)",nMomBins,0.,maxP);
  TH1F* hImpParXRms = new TH1F("hImpParXRms",";p (GeV/c);rms (Imp Par X) (#mum)",nMomBins,0.,maxP);
  TH1F* hImpParYRms = new TH1F("hImpParYRms",";p (GeV/c);rms (Imp Par Y) (#mum)",nMomBins,0.,maxP);
  TH1F* hImpParXSig = new TH1F("hImpParXSig",";p (GeV/c);#sigma(Imp Par X) (#mum)",nMomBins,0.,maxP);
  TH1F* hImpParYSig = new TH1F("hImpParYSig",";p (GeV/c);#sigma(Imp Par Y) (#mum)",nMomBins,0.,maxP);
  FillMeanAndRms(hImpParXVsP,hImpParXMean,hImpParXRms,hImpParXSig);
  FillMeanAndRms(hImpParYVsP,hImpParYMean,hImpParYRms,hImpParYSig);

  TH1F* hDeltaPMean = new TH1F("hDeltaPMean",";p (GeV/c);<p_{rec}-p_{gen}> (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPRms = new TH1F("hDeltaPRms",";p (GeV/c);rms (p_{rec}-p_{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPSig = new TH1F("hDeltaPSig",";p (GeV/c);#sigma(p_{rec}-p_{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hRelDeltaPMean = new TH1F("hRelDeltaPMean",";p (GeV/c);<(p_{rec}-p_{gen})/p_{gen}>",nMomBins,0.,maxP);
  TH1F* hRelDeltaPRms = new TH1F("hRelDeltaPRms",";p (GeV/c);rms (p_{rec}-p_{gen})/p_{gen}",nMomBins,0.,maxP);
  TH1F* hRelDeltaPSig = new TH1F("hRelDeltaPSig",";p (GeV/c);#sigma(p_{rec}-p_{gen})/p_{gen} (GeV/c)",nMomBins,0.,maxP);
  FillMeanAndRms(hDeltaPVsP,hDeltaPMean,hDeltaPRms,hDeltaPSig);
  FillMeanAndRms(hRelDeltaPVsP,hRelDeltaPMean,hRelDeltaPRms,hRelDeltaPSig);
  
  TH1F* hDeltaPxMean = new TH1F("hDeltaPxMean",";p (GeV/c);<p_{x}^{rec}-p_{x}^{gen}> (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPxRms = new TH1F("hDeltaPxRms",";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPxSig = new TH1F("hDeltaPxSig",";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hRelDeltaPxMean = new TH1F("hRelDeltaPxMean",";p (GeV/c);<(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}>",nMomBins,0.,maxP);
  TH1F* hRelDeltaPxRms = new TH1F("hRelDeltaPxRms",";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}",nMomBins,0.,maxP);
  TH1F* hRelDeltaPxSig = new TH1F("hRelDeltaPxSig",";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen} (GeV/c)",nMomBins,0.,maxP);
  FillMeanAndRms(hDeltaPxVsP,hDeltaPxMean,hDeltaPxRms,hDeltaPxSig);
  FillMeanAndRms(hRelDeltaPxVsP,hRelDeltaPxMean,hRelDeltaPxRms,hRelDeltaPxSig);

  TH1F* hDeltaPyMean = new TH1F("hDeltaPyMean",";p (GeV/c);<p_{y}^{rec}-p_{x}^{gen}> (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPyRms = new TH1F("hDeltaPyRms",";p (GeV/c);rms (p_{y}^{rec}-p_{x}^{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPySig = new TH1F("hDeltaPySig",";p (GeV/c);#sigma(p_{y}^{rec}-p_{x}^{gen}) (GeV/c)",nMomBins,0.,maxP);
   TH1F* hRelDeltaPyMean = new TH1F("hRelDeltaPyMean",";p (GeV/c);<(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}>",nMomBins,0.,maxP);
  TH1F* hRelDeltaPyRms = new TH1F("hRelDeltaPyRms",";p (GeV/c);rms (p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}",nMomBins,0.,maxP);
  TH1F* hRelDeltaPySig = new TH1F("hRelDeltaPySig",";p (GeV/c);#sigma(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen} (GeV/c)",nMomBins,0.,maxP);
  FillMeanAndRms(hDeltaPyVsP,hDeltaPyMean,hDeltaPyRms,hDeltaPySig);
  FillMeanAndRms(hRelDeltaPyVsP,hRelDeltaPyMean,hRelDeltaPyRms,hRelDeltaPySig);

  TH1F* hDeltaPzMean = new TH1F("hDeltaPzMean",";p (GeV/c);<p_{z}^{rec}-p_{z}^{gen}> (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPzRms = new TH1F("hDeltaPzRms",";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hDeltaPzSig = new TH1F("hDeltaPzSig",";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen}) (GeV/c)",nMomBins,0.,maxP);
  TH1F* hRelDeltaPzMean = new TH1F("hRelDeltaPzMean",";p (GeV/c);<(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}>",nMomBins,0.,maxP);
  TH1F* hRelDeltaPzRms = new TH1F("hRelDeltaPzRms",";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}",nMomBins,0.,maxP);
  TH1F* hRelDeltaPzSig = new TH1F("hRelDeltaPzSig",";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen} (GeV/c)",nMomBins,0.,maxP);
  FillMeanAndRms(hDeltaPzVsP,hDeltaPzMean,hDeltaPzRms,hDeltaPzSig);
  FillMeanAndRms(hRelDeltaPzVsP,hRelDeltaPzMean,hRelDeltaPzRms,hRelDeltaPzSig);
  
  TCanvas* ct = new TCanvas("ct","",1000,800);
  ct->Divide(2,2);
  ct->cd(1);
  hIsGood->Draw();
  ct->cd(2);
  hNclu->Draw();
  ct->cd(3);
  hChi2NDF->Draw();

  TCanvas* cip = new TCanvas("cip","",1200,800);
  cip->Divide(2,2);
  cip->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXMean->SetMarkerStyle(25);
  hImpParXMean->SetMarkerColor(1);
  hImpParXMean->SetLineColor(1);
  hImpParXMean->SetMinimum(-20);
  hImpParXMean->SetMaximum(20);
  hImpParXMean->Draw("P");
  cip->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXSig->SetMarkerStyle(25);
  hImpParXSig->SetMarkerColor(1);
  hImpParXSig->SetLineColor(1);
  hImpParXSig->SetMinimum(0);
  hImpParXSig->SetMaximum(120);
  hImpParXSig->Draw("P");
  cip->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYMean->SetMarkerStyle(25);
  hImpParYMean->SetMarkerColor(1);
  hImpParYMean->SetLineColor(1);
  hImpParYMean->SetMinimum(-20);
  hImpParYMean->SetMaximum(20);
  hImpParYMean->Draw("P");
  cip->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYSig->SetMarkerStyle(25);
  hImpParYSig->SetMarkerColor(1);
  hImpParYSig->SetLineColor(1);
  hImpParYSig->SetMinimum(0);
  hImpParYSig->SetMaximum(120);
  hImpParYSig->Draw("P");

  TCanvas* cmom = new TCanvas("cmom","",1400,900);
  cmom->Divide(3,3);
  cmom->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxMean->SetMarkerStyle(25);
  hDeltaPxMean->SetMarkerColor(1);
  hDeltaPxMean->SetLineColor(1);
  hDeltaPxMean->SetMinimum(-0.05);
  hDeltaPxMean->SetMaximum(0.05);
  hDeltaPxMean->Draw("P");
  cmom->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxSig->SetMarkerStyle(25);
  hDeltaPxSig->SetMarkerColor(1);
  hDeltaPxSig->SetLineColor(1);
  hDeltaPxSig->SetMinimum(0);
  hDeltaPxSig->SetMaximum(0.05);
  hDeltaPxSig->Draw("P");
  cmom->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPxSig->SetMarkerStyle(25);
  hRelDeltaPxSig->SetMarkerColor(1);
  hRelDeltaPxSig->SetLineColor(1);
  hRelDeltaPxSig->SetMinimum(0);
  hRelDeltaPxSig->SetMaximum(0.1);
  hRelDeltaPxSig->Draw("P");
  cmom->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPyMean->SetMarkerStyle(25);
  hDeltaPyMean->SetMarkerColor(1);
  hDeltaPyMean->SetLineColor(1);
  hDeltaPyMean->SetMinimum(-0.05);
  hDeltaPyMean->SetMaximum(0.05);
  hDeltaPyMean->Draw("P");
  cmom->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPySig->SetMarkerStyle(25);
  hDeltaPySig->SetMarkerColor(1);
  hDeltaPySig->SetLineColor(1);
  hDeltaPySig->SetMinimum(0);
  hDeltaPySig->SetMaximum(0.05);
  hDeltaPySig->Draw("P");
  cmom->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPySig->SetMarkerStyle(25);
  hRelDeltaPySig->SetMarkerColor(1);
  hRelDeltaPySig->SetLineColor(1);
  hRelDeltaPySig->SetMinimum(0);
  hRelDeltaPySig->SetMaximum(0.1);
  hRelDeltaPySig->Draw("P");
  cmom->cd(7);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzMean->SetMarkerStyle(25);
  hDeltaPzMean->SetMarkerColor(1);
  hDeltaPzMean->SetLineColor(1);
  hDeltaPzMean->SetMinimum(-0.05);
  hDeltaPzMean->SetMaximum(0.05);
  hDeltaPzMean->Draw("P");
  cmom->cd(8);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzSig->SetMarkerStyle(25);
  hDeltaPzSig->SetMarkerColor(1);
  hDeltaPzSig->SetLineColor(1);
  hDeltaPzSig->SetMinimum(0);
  hDeltaPzSig->SetMaximum(0.2);
  hDeltaPzSig->Draw("P");
  cmom->cd(9);
  hRelDeltaPzSig->SetMarkerStyle(25);
  hRelDeltaPzSig->SetMarkerColor(1);
  hRelDeltaPzSig->SetLineColor(1);
  hRelDeltaPzSig->SetMinimum(0);
  hRelDeltaPzSig->SetMaximum(0.05);
  hRelDeltaPzSig->Draw("P");

  TCanvas* ce = new TCanvas("ce","",1200,600);
  ce->Divide(2,1);
  ce->cd(1);
  hEtaGen->SetLineColor(kGray);
  hEtaGen->SetLineWidth(2);
  hEtaGen->Draw();
  hEtaReco->SetLineColor(kGreen+1);
  hEtaReco->SetLineWidth(2);
  hEtaReco->Draw("same");
  ce->cd(2);
  TH1F* hEffRec = (TH1F*)hEtaReco->Clone("hEffRec");
  hEffRec->Divide(hEtaReco,hEtaGen,1.,1.,"B");
  hEffRec->SetMaximum(1.05);
  hEffRec->Draw("");

  TFile* outFFit = new TFile("TrackingFastFit.root","recreate");
  hImpParXVsP->Write();
  hImpParYVsP->Write();
  hDeltaPxVsP->Write();
  hDeltaPyVsP->Write();
  hDeltaPzVsP->Write();
  hDeltaPVsP->Write();
  hRelDeltaPxVsP->Write();
  hRelDeltaPyVsP->Write();
  hRelDeltaPzVsP->Write();
  hRelDeltaPVsP->Write();
  hEtaGen->Write();
  hEtaReco->Write();
  outFFit->Close();

}
