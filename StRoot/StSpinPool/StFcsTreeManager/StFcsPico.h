/*
  The purpose of this header is to represent the basic "pico" level information about FCS information from simulation and data. The other reason to keep all these classes separate is so that you can load this library into ROOT without any other STAR libraries.
  
  @[June 13, 2023](David Kapukchyan) > Copied from StFcsPicoTree.h and only left the relevant "pico" classes
 */

#ifndef StFcsPico_H
#define StFcsPico_H

//ROOT headers
#include "TObject.h"

//A dummy class to hold just basic FCS hit information
class StFcsPicoHit : public TObject{
public:
  StFcsPicoHit() {}
  
  UShort_t mDetId = 0;
  UShort_t mChId = 0;
  UInt_t mAdcSum = 0;
  Float_t mEnergy = 0;

  //STAR coordinates
  Float_t mXstar = 0;
  Float_t mYstar = 0;
  Float_t mZstar = 0;
  
  ClassDef(StFcsPicoHit,1);
};

//A dummy class to hold just basic FCS cluster information
class StFcsPicoCluster : public TObject{
public:
  StFcsPicoCluster(){}
  
  Int_t mId = -1;  //Cluster Id
  UShort_t mDetId = 0;
  Int_t mCategory=0;
  Int_t mNTowers=0;
  Int_t mNNeighbor=0;
  Int_t mNPoints=0;
  Float_t mEnergy = 0;
  Float_t mX=0;         //This is local x
  Float_t mY=0;         //This is local y
  Float_t mSigmaMin=0;
  Float_t mSigmaMax=0;
  Float_t mTheta=0;
  Float_t mChi2Ndf1Photon=0;
  Float_t mChi2Ndf2Phoron=0;

  //Lorentz 4 momentum of cluster
  Double_t mLorentzX = 0;
  Double_t mLorentzY = 0;
  Double_t mLorentzZ = 0;
  Double_t mLorentzE = 0;

  ClassDef(StFcsPicoCluster,1);
};

//A dummy clas to hold just the basic FCS point information
class StFcsPicoPoint : public TObject{
public:
  StFcsPicoPoint(){}

  Short_t mNS = -1;
  Float_t mEnergy = 0;
  Float_t mXlocal=0;
  Float_t mYlocal=0;
  Int_t mNParentClusterPhotons=0;

  //STAR xyx
  Double_t mXstar = 0;
  Double_t mYstar = 0;
  Double_t mZstar = 0;

  //Lorentz 4 momentum of point
  Double_t mLorentzX = 0;
  Double_t mLorentzY = 0;
  Double_t mLorentzZ = 0;
  Double_t mLorentzE = 0;  

  ClassDef(StFcsPicoPoint,1);
};

//A dummy  class to hold just the basic G2t track information
class StPicoG2tTrack : public TObject
{
public:
  StPicoG2tTrack() {}

  long id = 0;
  long ge_pid  = 0;
  double mPx = 0;
  double mPy = 0;
  double mPz = 0;
  double mE  = 0;
  double mEta = 0;
  double phi()  { return atan2(mPy,mPx); }
  double pt()   { return sqrt( mPx*mPx + mPy*mPy ); }
  double ptot() { return sqrt( mPx*mPx + mPy*mPy + mPz*mPz ); }
  double theta(){ return 2.0*atan(exp(-1.0*mEta)); }
  double mXProj = 0;
  double mYProj = 0;
  double mZProj = 0;

  ClassDef(StPicoG2tTrack,1);
};



#endif
