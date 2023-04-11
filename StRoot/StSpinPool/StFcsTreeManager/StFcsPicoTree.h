/*
  The purpose of this class is to represent our own *TTree* that will contain more basic "pico" level information about FCS data from the MuDSTs
  
  @[April 3, 2023](David Kapukchyan) > First instance
 */

#ifndef StFcsPicoTree_H
#define StFcsPicoTree_H

//ROOT headers
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

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

class StFcsPicoTree : public TTree
{
  
 public:
  StFcsPicoTree();
  StFcsPicoTree(const char* name, const char* title, Int_t splitlevel=99);
  StFcsPicoTree(const char* name, TFile* file);
  virtual ~StFcsPicoTree();
  //StFcsPicoTree(const StFcsPicoTree& tm) = delete;
  //StFcsPicoTree& operator=(const StFcsPicoTree& tm) = delete;

  void LoadFile(TFile* file, const char* treename);    //!< Grab the TClonesArray objects from a differently named StFcsPicoTree and load them into *this* tree  

  TClonesArray* GetHits()    const { return mHits; }
  TClonesArray* GetClusters()const { return mClusters; }
  TClonesArray* GetPoints()  const { return mPoints; }

  StFcsPicoHit* ConstructedHit(Int_t ihit);            //!< If ihit exists in #mHits then return it otherwise create one, i.e. calls mHits->ConstructedAt(ihit);
  StFcsPicoCluster* ConstructedCluster(Int_t iclus);   //!< If iclus exists in #mClusters then return it otherwise create one, i.e. calls mHits->ConstructedAt(ihit);
  StFcsPicoPoint* ConstructedPoint(Int_t ipoint);      //!< If ipoint exists in #mPoints then return it otherwise create one, i.e. calls mHits->ConstructedAt(ihit);

  void ClearAll();        //!< Clear #mHits, #mClusters, and #mPoints
  void ClearHits();       //!< Clear #mHits
  void ClearClusters();   //!< Clear #mClusters
  void ClearPoints();     //!< Clear #mPoints

  void DeleteAll();       //!< Check and delete all TClonesArrays
  void DeleteHits();      //!< Check and delete #mHits
  void DeleteClusters();  //!< Check and delete #mClusters
  void DeletePoints();    //!< Check and delete #mPoints
  
    
 protected:
  TClonesArray* mHits = 0;         //!< Array of #StFcsPicoHit
  TClonesArray* mClusters = 0;     //!< Array of #StFcsPicoCluster
  TClonesArray* mPoints = 0;       //!< Array of #StFcsPicoPoint

  void InitBranches();             //!< Set up the branch objects for *this* tree
  
  ClassDef(StFcsPicoTree,1);
};

#endif
