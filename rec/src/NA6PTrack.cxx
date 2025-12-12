// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TGeoGlobalMagField.h>
#include "NA6PBaseCluster.h"
#include "NA6PTrack.h"

ClassImp(NA6PTrack)

//_______________________________________________________________________
NA6PTrack::NA6PTrack() :
  mMass{0.140},
  mChi2{0.0},
  mNLayers{6},
  mClusterMap{0},
  mNClusters{0},
  mNClustersVT{0},
  mNClustersMS{0},
  mNClustersTR{0},
  mClusterIndices{},
  mClusterPartID{},
  mParticleID{-2},
  mExtTrack{}
{
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
}

//_______________________________________________________________________
NA6PTrack::NA6PTrack(const double* xyz, const double* pxyz, int sign, double errLoose) :
  mMass{0.140},
  mChi2{0.0},
  mNLayers{6},
  mClusterMap{0},
  mNClusters{0},
  mNClustersVT{0},
  mNClustersMS{0},
  mNClustersTR{0},
  mClusterIndices{},
  mClusterPartID{},
  mParticleID{-2},
  mExtTrack{}
{
  // initialize arrays
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  // initialize the track parameters and covariance
  init(xyz, pxyz, sign, errLoose);
}

//_______________________________________________________________________
void NA6PTrack::reset() 
{
  mMass = 0.14; 
  mChi2 = 0; 
  mClusterMap = 0;
  mExtTrack.Reset();
  resetCovariance();
  mParticleID = -2;
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  mNClusters = mNClusters = mNClustersMS = mNClustersTR = 0;
}
  
//_______________________________________________________________________
Bool_t NA6PTrack::init(const double *xyz, const double *pxyz, int sign, double errLoose)
{
  // Init with track position/momentum in usual Lab frame
  // If errLoose>0 then scale initially small errors by this amount
  double xyzL[3],pxyzL[3];
  lab2trk(xyz,xyzL);
  lab2trk(pxyz,pxyzL);
  double cov[21] = {1.e-6,    // assign small errors first
		    0.   ,1.e-6, 
		    0.   ,0.   ,1.e-6, 
		    0.   ,0.   ,0.   ,1.e-4,
		    0.   ,0.   ,0.   ,0.    ,1.e-4, 
		    0.   ,0.   ,0.   ,0.    ,0.   ,1.e-3};

  mExtTrack.set(xyzL,pxyzL,cov,sign);
  if (errLoose>0) mExtTrack.ResetCovariance(errLoose);
  mExtTrack.Rotate(pxyz[2] > 0 ? 0.:TMath::Pi());
  return true;
}

//_______________________________________________________________________
Bool_t NA6PTrack::propagateToZBxByBz(double z, double maxDZ, double xOverX0, double xTimesRho)
{
  // propagate the track to position Z in uniform material with xOverX0 rad lgt and xTimesRho lgt*density
  //
  double zCurr = getZLab();
  double dz = z - zCurr;
  if (TMath::Abs(dz)<kAlmost0) return true;
  int nz = TMath::Abs(dz)/maxDZ + 1;
  double zstep = dz/nz;
  double xyz[3],bxyz[3],bxyzFwd[3];
  for (int iz=0;iz<nz;iz++) {
    getXYZ(xyz);              // coordinates in Lab frame
    TGeoGlobalMagField::Instance()->Field(xyz,bxyz);
    lab2trk(bxyz,bxyzFwd);  // field in Fwd frame
    zCurr += zstep;
    if (!propagateToZBxByBz(zCurr,bxyzFwd)) return false;
    if (TMath::Abs(xTimesRho)>1e-6 && 
	!correctForMeanMaterial(xOverX0/nz, xTimesRho/nz)) return false;
  }
  return true;
  //
}

//_______________________________________________________________________
Bool_t NA6PTrack::propagateToDCA(NA6PTrack* partner)
{
  // propagate the track to position Z of closest approach to partner track
  //
  double xyz[3],bxyz[3],bxyzFwd[3];
  getXYZ(xyz);              // coordinates in Lab frame
  TGeoGlobalMagField::Instance()->Field(xyz,bxyz);
  lab2trk(bxyz,bxyzFwd);  // field in Fwd frame
  double zthis=0, zpartner=0;
  double dca=mExtTrack.getDCA(&partner->mExtTrack,bxyzFwd[2],zthis,zpartner);

  if (!propagateToZBxByBz(zthis) || !partner->propagateToZBxByBz(zpartner)) return false;
  return true;
  //
}


//_______________________________________________________________________
double NA6PTrack::getSigmaP2() const
{
  // error^2 on total momentum, P = sqrt(1+tgl^2)/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  double tgl  = mExtTrack.getTgl();
  double tglE  = mExtTrack.getSigmaTgl2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvtgE = mExtTrack.getSigma1PtTgl();
  //
  double tp12 = TMath::Sqrt(1.+tgl*tgl);
  double dt =  tgl/pinv/tp12;  // dP/dtgl
  double dc = -tp12/pinv/pinv; // dP/dC
  double err2 = dt*dt*tglE +dc*dc*pinvE + dt*dc*pinvtgE;
  return err2;
  //
}

//_______________________________________________________________________
double NA6PTrack::getSigmaPX2() const
{
  // error^2 on Px in Lav frame (Z in the tracking frame, = tgl/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  double tgl  = mExtTrack.getTgl();
  double tglE  = mExtTrack.getSigmaTgl2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvtgE = mExtTrack.getSigma1PtTgl();
  //
  double dt =  1./pinv;      // dP/dtgl
  double dc = -tgl/pinv/pinv; // dP/dC
  double err2 = dt*dt*tglE +dc*dc*pinvE + dt*dc*pinvtgE;
  return err2;
  //
}

//_______________________________________________________________________
double NA6PTrack::getSigmaPY2() const
{
  // error^2 on Py in Lab frame (Y in the tracking frame, = sinp/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  double cosAlp=TMath::Cos(mExtTrack.getAlpha()), sinAlp=TMath::Sin(mExtTrack.getAlpha());
  double snp  = mExtTrack.getSnp();
  double csp  = TMath::Sqrt((1. - snp)*(1. + snp));
  double snpE  = mExtTrack.getSigmaSnp2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvsnpE = mExtTrack.getSigma1PtSnp();
  //
  double ds =  (cosAlp-sinAlp*snp/csp)/pinv; // dP/dsnp
  double dc = -(snp*cosAlp+csp*sinAlp)/pinv/pinv; // dP/dC
  double err2 = ds*ds*snpE +dc*dc*pinvE + ds*dc*pinvsnpE;
  return err2;

  //
}

//_______________________________________________________________________
double NA6PTrack::getSigmaPZ2() const
{
  // error^2 on Pz in Lab frame (X in the tracking frame, = (sqrt(1-snp^2)*cosAlp-snp*sinAlp)/(1/pt)
  double pinv = mExtTrack.getSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  double cosAlp=TMath::Cos(mExtTrack.getAlpha()), sinAlp=TMath::Sin(mExtTrack.getAlpha());
  double snp  = mExtTrack.getSnp();
  double snpE  = mExtTrack.getSigmaSnp2();
  double pinvE = mExtTrack.getSigma1Pt2();
  double pinvsnpE = mExtTrack.getSigma1PtSnp();
  //
  double csp = TMath::Sqrt( (1.-snp)*(1.+snp) );
  double ds = -snp/pinv*(cosAlp/csp + sinAlp);        // dP/dsnp
  double dc = -(csp*cosAlp-snp*sinAlp)/pinv/pinv;     // dP/dC
  double err2 = ds*ds*snpE +dc*dc*pinvE + ds*dc*pinvsnpE;
  return err2;
  //
}

//_______________________________________________________________________
void NA6PTrack::resetCovariance(float err)
{
  // reset cov matrix
  double *trCov  = mExtTrack.getCovarianceForUpdate();
  const double *trPars = mExtTrack.getParameter();
  const double kLargeErr2Coord = 50*50;
  const double kLargeErr2Dir = 0.6*0.6;
  const double kLargeErr2PtI = 2.;
  for (int ic=15;ic--;) trCov[ic] = 0.;
  trCov[kY2]   = trCov[kZ2]   = err<0 ? kLargeErr2Coord : err*err; 
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
  mExtTrack.checkCovariance();
}

//_______________________________________________________________________
std::string NA6PTrack::asString() const
{
  double pxyz[3];
  getPXYZ(pxyz);
  return fmt::format("Track: Nclusters:{} NVTclusters:{}  chi2:{} pos:{:.4f},{:.4f},{:.4f} mom:{:.3f},{:.3f},{:.3f}",
                     mNClusters,mNClusters,mChi2,getXLab(),getYLab(),getZLab(),pxyz[0],pxyz[1],pxyz[2]);
}

//_______________________________________________________________________
void NA6PTrack::print() const
{
  LOGP(info, "{}", asString());
}

//_______________________________________
void NA6PTrack::addCluster(const NA6PBaseCluster* clu, int cluIndex, double chi2) {

  mNClusters++;
  mChi2 += chi2;
  // int nDet = clu->getDetectorID();
  int trackID = clu->getParticleID();
  
  int nLay = clu->getLayer();
  mClusterPartID[nLay] = trackID;
  mClusterIndices[nLay] = cluIndex;
  mClusterMap |= (1<<nLay);

  // if (nDet < mNLayers * 4) {
  // }
  // else if (nDet >= 20 && nDet < 24) mNClustersMS++;
  // else if (nDet >= 24) mNClustersTR++;
}

//____________________________________________
void NA6PTrack::imposeKinematics(const double* xyzLab,const double* cosinesLab,
				   double en, double mass, int charge) 
{
  // RS: note: we assume p=e here to avoid problem with e+ e- interpreted as muon
  double p = en*en - mass*mass;
  if (p<=0) {
    printf("Anomalous kinematics: E:%e M:%e",en,mass);
    exit(1);
  }
  p = TMath::Sqrt(p);
  double pxyz[3] = {p*cosinesLab[0],p*cosinesLab[1],p*cosinesLab[2]};
  init(xyzLab,pxyz,charge,1.e4);
  setMass(mass);
  resetCovariance();// reset cov.matrix
}
