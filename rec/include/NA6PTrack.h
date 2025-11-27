// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_TRACK_H
#define NA6P_TRACK_H

#include <string>
#include <Rtypes.h>
#include "ExtTrackPar.h"

// Basic track class

class NA6PBaseCluster;

class NA6PTrack
{
 public:
  
  enum {kNDOF=5, kMaxVTLr=16};
  enum {kY2=0,kZ2=2,kSnp2=5,kTgl2=9,kPtI2=14};
  enum {kY,kZ,kSnp,kTgl,kPtI};

  NA6PTrack();
  NA6PTrack(const double *xyz, const double *pxyz, int sign, double errLoose=-1);
  NA6PTrack(const NA6PTrack&) = default;
  NA6PTrack& operator=(const NA6PTrack&) = default;
  virtual ~NA6PTrack() {}

  bool    init(const double *xyz, const double *pxyz, int sign, double errLoose=-1);
  void    imposeKinematics(const double* xyzLab,const double* cosinesLab, double en, double mass, int charge);
  void    reset(); 
  void    resetCovariance(float err=-1.);

  double        getMass()                  const {return mMass;}
  int           getNVTLayers()             const {return mNVTLayers;}
  double        getChi2()                  const {return mChi2;}
  double        getChi2VT()                const {return mChi2VT;}
  const ExtTrackPar& getTrackExtParam()    const {return mExtTrack; }
  const double* getCovariance()            const {return mExtTrack.getCovariance(); }
  int          getNVTHits()                const {return mNClustersVT;}
  int          getNMSHits()                const {return mNClustersMS;}
  int          getNTRHits()                const {return mNClustersTR;}
  int          getNHits()                  const {return mNClusters;}
  uint32_t     getClusterMap()             const {return mClusterMap;}
  bool         hasClusterOnVTLayer(int lr) const {return (lr<mNVTLayers) ?  mClusterMap&(1<<lr) : false;}
  int          getVTClusterIndex(int lr)   const {return lr < kMaxVTLr ? mClusterIndices[lr] : -1; }
  int          getParticleLabel(int lr)    const {return lr < kMaxVTLr ? mClusterPartID[lr] : -2; }
  int          getParticleID()             const {return mParticleID;}
  int          getCAIteration()            const {return mCAIteration;}
  
  double       getAlpha()                  const {return mExtTrack.getAlpha();}
  double       getCharge()                 const {return mExtTrack.Charge();}
  bool         negDir()                    const {return std::abs(mExtTrack.getAlpha()) > M_PI/2.;}
  double       getXLab()                   const {return mExtTrack.getY();} 
  double       getYLab()                   const {return mExtTrack.getZ();} 
  double       getZLab()                   const {return negDir() ? -mExtTrack.getX():mExtTrack.getX();}
  double       getR()                      const {double x=getXLab(),y=getYLab(),r=x*x+y*y; return r>0? std::sqrt(r):0;}
  double       getXTF()                    const {return mExtTrack.getX();}
  double       getYTF()                    const {return mExtTrack.getY();}
  double       getZTF()                    const {return mExtTrack.getZ();}
  void         getXYZ(double *xyz)         const;
  void         getPXYZ(double *pxyz)       const;  
  double       getP()                      const {return mExtTrack.getP();}
  double       getSigmaX2()                const {return mExtTrack.getSigmaY2();}
  double       getSigmaY2()                const {return mExtTrack.getSigmaZ2();}
  double       getSigmaXY()                const {return mExtTrack.getSigmaZY();}  
  double       getSigmaP2()                const;
  double       getSigmaPX2()               const;
  double       getSigmaPY2()               const;
  double       getSigmaPZ2()               const;
  double       getNormChi2()               const {return mNClusters<3 ? 0 :  mChi2 / ( (mNClusters<<1)-kNDOF);}
  double       getNormChi2VT()             const {return mNClustersVT<3 ? 0 :  mChi2VT / ( (mNClustersVT<<1)-kNDOF);}
  double       getPredictedChi2(double* p, double* cov) const {return mExtTrack.getPredictedChi2(p,cov);}
  
  void   setMass(double m)         {mMass = m;}
  void   setNVTLayers(int n)       {mNVTLayers = n;}
  void   setChi2(double chi2)      {mChi2 = chi2;}
  void   setChi2VT(double chi2)    {mChi2VT = chi2;}
  void   setParticleLabel(int idx, int lr)  { if (lr < kMaxVTLr) mClusterPartID[lr] = idx;}
  void   setVTClusterIndex(int idx, int lr) { if (lr < kMaxVTLr) mClusterIndices[lr] = idx;}
  void   setParticleID(int idx)    {mParticleID = idx;}
  void   setCAIteration(int iter)  {mCAIteration = iter;}

  static void lab2trk(const double *vLab, double *vTrk); 
  static void trk2lab(const double *vTrk, double *vLab); 

  void   addCluster(const NA6PBaseCluster* clu, int cluIndex, double chi2);
  bool   correctForMeanMaterial(double xOverX0, double xTimesRho, bool anglecorr = false){
    return mExtTrack.correctForMeanMaterial(xOverX0, xTimesRho, mMass, anglecorr);
  }
  bool   propagateToZBxByBz(double z, double maxDZ=1.0, double xOverX0=0., double xTimesRho=0.);
  bool   propagateToZBxByBz(double z, const double *bxyz){return mExtTrack.propagateToBxByBz(z,bxyz);}
  bool   propagateToDCA(NA6PTrack* partner);
  bool   update(double p[2],double cov[3]) {return mExtTrack.Update(p,cov);}

  virtual void  print() const;
  std::string asString() const;

  
 protected:
  double   mMass = 0.140;                        // particle mass
  double   mChi2 = 0.f;                          // total chi2
  double   mChi2VT = 0.f;                        // total chi2 VT
  int      mNVTLayers = 5;                       // number of vT laters
  uint32_t mClusterMap = 0;                      // pattern of clusters per layer
  int      mNClusters = 0;                       // total hits
  int      mNClustersVT = 0;                     // total VT hits
  int      mNClustersMS = 0;                     // total MS hits
  int      mNClustersTR = 0;                     // total TR hits
  std::array<int, kMaxVTLr> mClusterIndices{};   // cluster indices
  std::array<int, kMaxVTLr> mClusterPartID{};    // particle ID (per cluster)
  int      mParticleID = -1;                     // particle ID (MC truth)
  int      mCAIteration = -1;                    //! CA iteration (for debug)
  ExtTrackPar mExtTrack;                         // track params
  
  ClassDefNV(NA6PTrack,1)
};

//_______________________________________________________________________
inline void  NA6PTrack::lab2trk(const double *vLab, double *vTrk)
{
  // convert alice coordinates to modified
  vTrk[0] = vLab[2];
  vTrk[1] = vLab[0];
  vTrk[2] = vLab[1];
}

//_______________________________________________________________________
inline void  NA6PTrack::trk2lab(const double *vTrk, double *vLab)
{
  // convert modified coordinates to Lab ones
  vLab[0] = vTrk[1];
  vLab[1] = vTrk[2];
  vLab[2] = vTrk[0];
}

//_______________________________________________________________________
inline void NA6PTrack::getXYZ(double *xyz) const
{
  // track position in Lab coordinates
  double xyzTrk[3];
  mExtTrack.getXYZ(xyzTrk);
  trk2lab(xyzTrk, xyz);
}

//_______________________________________________________________________
inline void NA6PTrack::getPXYZ(double *pxyz) const
{
  // track position in Lab coordinates
  double pxyzTrk[3];
  mExtTrack.getPxPyPz(pxyzTrk);
  trk2lab(pxyzTrk, pxyz);
}


#endif

