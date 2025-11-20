#ifndef EXTTRACKPAR_H
#define EXTTRACKPAR_H
/*****************************************************************************
 *              "External" track parametrisation class                       *
 *                                                                           *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuthal angle  *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *                                                                           *
 * The parameters are estimated at an exact position x in a local coord.     *
 * system rotated by angle alpha with respect to the global coord.system.    *
 *        Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                     *
 *****************************************************************************/

#include "TPolyMarker3D.h"
#include "TMath.h"

#ifndef ALIEXTERNALTRACKPARAM_H

const Double_t kAlmost0=Double_t(FLT_MIN);
const Double_t kAlmost1=1. - kAlmost0;
const Double_t kVeryBig=1./kAlmost0;
const Double_t kMostProbablePt=0.35;
const Double_t kB2C=-0.299792458e-3;
const Double_t kAlmost0Field=1.e-13;

const Double_t kC0max=100*100, // SigmaY<=100cm
               kC2max=100*100, // SigmaZ<=100cm
               kC5max=1*1,     // SigmaSin<=1
               kC9max=1*1,     // SigmaTan<=1
               kC14max=100*100; // Sigma1/Pt<=100 1/GeV
#endif

class ExtTrackPar: public TObject {
 public:
  ExtTrackPar();
  ExtTrackPar(const ExtTrackPar &);
  ExtTrackPar& operator=(const ExtTrackPar & trkPar);
  ExtTrackPar(Double_t x, Double_t alpha, 
			const Double_t param[5], const Double_t covar[15]);
  ExtTrackPar(Double_t xyz[3],Double_t pxpypz[3],
			Double_t cv[21],Short_t sign);

  virtual ~ExtTrackPar(){}
  
  template <typename T>
  void set(T x, T alpha, const T param[5], const T covar[15]) {
    //  Sets the parameters
    if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
    else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();
    fX=x; fAlpha=alpha;
    for (Int_t i = 0; i < 5; i++)  fP[i] = param[i];
    for (Int_t i = 0; i < 15; i++) fC[i] = covar[i];

    checkCovariance();

  }

  void setParamOnly(double x, double alpha, const double param[5]) {
    //  Sets the parameters, neglect cov matrix
    if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
    else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();
    fX=x; fAlpha=alpha;
    for (Int_t i = 0; i < 5; i++)  fP[i] = param[i];
  }

  void set(Double_t xyz[3],Double_t pxpypz[3],Double_t cv[21],Short_t sign);

  static void setMostProbablePt(Double_t pt) { fgMostProbablePt=pt; }
  static Double_t getMostProbablePt() { return fgMostProbablePt; }

  void Reset();
  void ResetCovariance(Double_t s2);
  void AddCovariance(const Double_t cov[15]);

  const Double_t *getParameter() const {return fP;}
  const Double_t *getCovariance() const {return fC;}
  Double_t *getCovarianceForUpdate() {return fC;} // to be used when need to update the cov mat
  virtual  Bool_t IsStartedTimeIntegral() const {return kFALSE;}
  virtual  void   AddTimeStep(Double_t ) {} // dummy method, real stuff is done in AliKalmanTrack
  Double_t getAlpha() const {return fAlpha;}
  Double_t getX() const {return fX;}
  Double_t getY()    const {return fP[0];}
  Double_t getZ()    const {return fP[1];}
  Double_t getSnp()  const {return fP[2];}
  virtual Double_t getTgl()  const {return fP[3];}
  Double_t getSigned1Pt()  const {return fP[4];}

  Double_t getSigmaY2() const {return fC[0];}
  Double_t getSigmaZY() const {return fC[1];}
  Double_t getSigmaZ2() const {return fC[2];}
  Double_t getSigmaSnpY() const {return fC[3];}
  Double_t getSigmaSnpZ() const {return fC[4];}
  Double_t getSigmaSnp2() const {return fC[5];}
  Double_t getSigmaTglY() const {return fC[6];}
  Double_t getSigmaTglZ() const {return fC[7];}
  Double_t getSigmaTglSnp() const {return fC[8];}
  Double_t getSigmaTgl2() const {return fC[9];}
  Double_t getSigma1PtY() const {return fC[10];}
  Double_t getSigma1PtZ() const {return fC[11];}
  Double_t getSigma1PtSnp() const {return fC[12];}
  Double_t getSigma1PtTgl() const {return fC[13];}
  Double_t getSigma1Pt2() const {return fC[14];}

  // additional functions for AliVParticle
  Double_t Px() const;
  Double_t Py() const;
  Double_t Pz() const { return Pt()*getTgl(); }
  Double_t Pt() const { return TMath::Abs(getSignedPt()); }
  Double_t P() const { return getP(); }
  Bool_t   PxPyPz(Double_t p[3]) const { return getPxPyPz(p); }
  
  Double_t Xv() const;
  Double_t Yv() const;
  Double_t Zv() const {return getZ();}
  Bool_t   XvYvZv(Double_t x[3]) const { return getXYZ(x); }

  Double_t OneOverPt() const { return 1./Pt(); }
  Double_t Phi() const;
  Double_t PhiPos() const;
  Double_t Theta() const;
  virtual Double_t E() const;
  virtual Double_t M() const;
  Double_t Eta() const;
  virtual Double_t Y() const;
  virtual Short_t  Charge() const { return (Short_t)getSign(); }
  virtual const Double_t *PID() const { return 0x0; }

  virtual Int_t    getID() const { return -999; }
  virtual UChar_t  getITSClusterMap() const {return 0; }
  virtual ULong_t  getStatus() const { return 0; }

  Double_t getSign() const {return (fP[4]>0) ? 1 : -1;}
  Double_t getP() const;
  Double_t getSignedPt() const {
    return (TMath::Abs(fP[4])>kAlmost0) ? 1./fP[4]:TMath::Sign(kVeryBig,fP[4]);
  }
  Double_t get1P() const;
  virtual Double_t getC(Double_t b) const {return fP[4]*b*kB2C;}
  void getDZ(Double_t x,Double_t y,Double_t z,Double_t b,Float_t dz[2]) const; 
  Double_t getD(Double_t xv, Double_t yv, Double_t b) const; 
  Double_t getLinearD(Double_t xv, Double_t yv) const; 

  Bool_t correctForMeanMaterial(Double_t xOverX0, Double_t xTimesRho, 
        Double_t mass,  Bool_t anglecorr=kFALSE,
	Double_t (*f)(Double_t)=ExtTrackPar::BetheBlochSolid);

  Bool_t correctForMeanMaterialdEdx(Double_t xOverX0, Double_t xTimesRho, 
	Double_t mass, Double_t dEdx, Bool_t anglecorr=kFALSE);

  Bool_t correctForMeanMaterialZA(Double_t xOverX0, Double_t xTimesRho, 
                                  Double_t mass,
                                  Double_t zOverA=0.49848,
                                  Double_t density=2.33,
                                  Double_t exEnergy=173e-9,
                                  Double_t jp1=0.20,
                                  Double_t jp2=3.00,
                                  Bool_t anglecorr=kFALSE
  );

  //
  // Bethe-Bloch formula parameterizations
  //
  static Double_t BetheBlochAleph(Double_t bg,
                                  Double_t kp1=0.76176e-1,
                                  Double_t kp2=10.632,
                                  Double_t kp3=0.13279e-4,
                                  Double_t kp4=1.8631,
                                  Double_t kp5=1.9479
				  );
  static Double_t BetheBlochGeant(Double_t bg,
                                  Double_t kp0=2.33,
                                  Double_t kp1=0.20,
                                  Double_t kp2=3.00,
                                  Double_t kp3=173e-9,
                                  Double_t kp4=0.49848
				  );
    
  static Double_t BetheBlochSolid(Double_t bg);
  static Double_t BetheBlochGas(Double_t bg);

  Double_t getPredictedChi2(const Double_t p[2],const Double_t cov[3]) const;

  Double_t 
    getPredictedChi2(const Double_t p[3],const Double_t covyz[3],const Double_t covxyz[3]) const;

  Double_t getPredictedChi2(const ExtTrackPar *t) const;

  Bool_t 
    propagateTo(Double_t p[3],Double_t covyz[3],Double_t covxyz[3],Double_t b);

  Double_t *getResiduals(Double_t *p,Double_t *cov,Bool_t updated=kTRUE) const;
  Bool_t Update(const Double_t p[2],const Double_t cov[3]);
  Bool_t Rotate(Double_t alpha);
  Bool_t RotateParamOnly(Double_t alpha);
  Bool_t Invert();
  Bool_t propagateTo(Double_t x, Double_t b);
  Bool_t propagateParamOnlyTo(Double_t xk, Double_t b);
  Bool_t propagate(Double_t alpha, Double_t x, Double_t b);
  Bool_t propagateBxByBz(Double_t alpha, Double_t x, Double_t b[3]);
  Bool_t propagateParamOnlyBxByBzTo(Double_t xk, const Double_t b[3]);
  void   propagate(Double_t len,Double_t x[3],Double_t p[3],Double_t bz) const;
  Bool_t Intersect(Double_t pnt[3], Double_t norm[3], Double_t bz) const;

  static void g3helx3(Double_t qfield, Double_t step, Double_t vect[7]); 
  Bool_t propagateToBxByBz(Double_t x, const Double_t b[3]);

  void getHelixParameters(Double_t h[6], Double_t b) const;
  Double_t getDCA(const ExtTrackPar *p, Double_t b,
    Double_t &xthis,Double_t &xp) const;
  Double_t propagateToDCA(ExtTrackPar *p, Double_t b);
  
  void getDirection(Double_t d[3]) const;
  Bool_t getPxPyPz(Double_t *p) const;  
  Bool_t getXYZ(Double_t *p) const;
  Bool_t getCovarianceXYZPxPyPz(Double_t cv[21]) const;
  Bool_t getPxPyPzAt(Double_t x, Double_t b, Double_t p[3]) const;
  Bool_t getXYZAt(Double_t x, Double_t b, Double_t r[3]) const;
  Double_t getParameterAtRadius(Double_t r, Double_t bz, Int_t parType) const;

  Bool_t getYAt(Double_t x,  Double_t b,  Double_t &y) const;
  Bool_t getZAt(Double_t x,  Double_t b,  Double_t &z) const;
  Double_t getYAtFast(Double_t x, Double_t b) const {double y=0; return getYAt(x,b,y) ? y : -99999;}
  Double_t getZAtFast(Double_t x, Double_t b) const {double z=0; return getZAt(x,b,z) ? z : -99999;}
  void Print(Option_t* option = "") const;
  Double_t getSnpAt(Double_t x,Double_t b) const;
  Bool_t getXatLabR(Double_t r,Double_t &x, Double_t bz, Int_t dir=0) const;
  Bool_t getXYZatR(Double_t xr,Double_t bz, Double_t *xyz=0, Double_t* alpSect=0) const;

  //Deprecated
  Bool_t correctForMaterial(Double_t d, Double_t x0, Double_t mass,
	 Double_t (*f)(Double_t)=ExtTrackPar::BetheBlochSolid);

  Bool_t getDistance(ExtTrackPar *param2, Double_t x, Double_t dist[3], Double_t b);
  static Int_t getIndex(Int_t i, Int_t j);
  Int_t getLabel() const {return -1;} 
  Int_t PdgCode()  const {return 0;}

  //
  // visualization (M. Ivanov)
  //
  virtual void FillPolymarker(TPolyMarker3D *pol, Float_t magf, Float_t minR, Float_t maxR, Float_t stepR);
  virtual void DrawTrack(Float_t magF, Float_t minR, Float_t maxR, Float_t stepR);

  virtual Bool_t Translate(Double_t *vTrasl,Double_t *covV);

  void checkCovariance();

  static Bool_t  getUseLogTermMS()                {return fgUseLogTermMS;} 
  static void    setUseLogTermMS(Bool_t v=kTRUE)  {fgUseLogTermMS = v;} 

  //---------------------------------------------------------------------------
  //--the calibration interface--
  //--to be used in online calibration/QA
  //--should also be implemented in ESD so it works offline as well
  //-----------
  virtual Int_t getTrackParam         ( ExtTrackPar & ) const {return 0;}
  virtual Int_t getTrackParamRefitted ( ExtTrackPar & ) const {return 0;}
  virtual Int_t getTrackParamIp       ( ExtTrackPar & ) const {return 0;}
  virtual Int_t getTrackParamTPCInner ( ExtTrackPar & ) const {return 0;}
  virtual Int_t getTrackParamOp       ( ExtTrackPar & ) const {return 0;}
  virtual Int_t getTrackParamCp       ( ExtTrackPar & ) const {return 0;}
  virtual Int_t getTrackParamITSOut   ( ExtTrackPar & ) const {return 0;}
  //
  // coordinate system conversions
  Bool_t   Local2GlobalMomentum(Double_t p[3], Double_t alpha) const;
  Bool_t   Local2GlobalPosition(Double_t r[3], Double_t alpha) const;
  Bool_t   Global2LocalMomentum(Double_t p[3], Short_t charge, Double_t &alpha) const;
  Bool_t   Global2LocalPosition(Double_t r[3], Double_t alpha) const;

  //
 protected:

/*  protected: */
 private:
  Double_t &Par(Int_t i) {return fP[i];}
  Double_t &Cov(Int_t i) {return fC[i];}
 protected:
  Double32_t           fX;     // X coordinate for the point of parametrisation
  Double32_t           fAlpha; // Local <-->global coor.system rotation angle
  Double32_t           fP[5];  // The track parameters
  Double32_t           fC[15]; // The track parameter covariance matrix

  static Double32_t    fgMostProbablePt; // "Most probable" pt
                                         // (to be used if Bz=0)
  static Bool_t        fgUseLogTermMS;   // use log term in Mult.Stattering evaluation
  ClassDef(ExtTrackPar, 8)
};

inline void ExtTrackPar::ResetCovariance(Double_t s2) {
  //
  // Reset the covarince matrix to "something big"
  //

  s2 = TMath::Abs(s2);
  Double_t fC0=fC[0]*s2,
           fC2=fC[2]*s2,
           fC5=fC[5]*s2,
           fC9=fC[9]*s2,
           fC14=fC[14]*s2;
 
  if (fC0>kC0max)  fC0 = kC0max;
  if (fC2>kC2max)  fC2 = kC2max;
  if (fC5>kC5max)  fC5 = kC5max;
  if (fC9>kC9max)  fC9 = kC9max;
  if (fC14>kC14max)  fC14 = kC14max;


    fC[0] = fC0;
    fC[1] = 0.;  fC[2] = fC2;
    fC[3] = 0.;  fC[4] = 0.;  fC[5] = fC5;
    fC[6] = 0.;  fC[7] = 0.;  fC[8] = 0.;  fC[9] = fC9;
    fC[10]= 0.;  fC[11]= 0.;  fC[12]= 0.;  fC[13]= 0.;  fC[14] = fC14;
}

#endif
