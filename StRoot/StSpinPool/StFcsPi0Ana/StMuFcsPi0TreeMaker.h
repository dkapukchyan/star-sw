/*
  AUTHOR
  David Kapukchyan

  SYNOPSIS 
  The purpose of this class is to generate a TTree of PI0s to use for transverse single spin asymmetry analysis.

  LOG
  @[March 5, 2024] > Copied from *StMuFcsPi0Maker* and modified to write pi0s to a tree

  @[May 24, 2024] > Added #FcsPi0Info class to hold information about the reconstructed Pi0s. Got rid of the histograms.

*/


#ifndef STMUFCSPI0TREEMAKER_HH
#define STMUFCSPI0TREEMAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

//STAR Headers
#include "StEnumerations.h"
#include "StMaker.h"
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
//#include "StEventTypes.h" (not in header)
//#include "StEvent/StEvent.h" (not in header)
//#include "StEvent/StFmsCollection.h" (not in header)
//#include "StEvent/StFmsPoint.h" (not in header)
//#include "StEvent/StFmsPointPair.h" (not in header)
#include "StEvent/StTriggerData.h"
#include "StEvent/StTriggerId.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "Stypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsHit.h"
#include "StMuDSTMaker/COMMON/StMuFcsCluster.h"
#include "StMuDSTMaker/COMMON/StMuFcsPoint.h"

//Class to hold basic info for reconstructed pi0 candidates
class FcsPi0Info : public TObject
{
  FcsPi0Info() {}
  ~FcsPi0Info() {}

  ULong_t mDetId = 0;   //!< Unique Pi0 identifier
  ULong_t mRun = 0;     //!< Run number where pi0 was found
  ULong_t mEvent = 0;   //!< STAR Event Id where pi0 was found
  UShort_t mSpin = 0;   //!< Spin state of proton when pi0 was reconstructed
  Double_t mPx = 0;     //!< X-Momentum
  Double_t mPy = 0;     //!< Y-Momentum
  Double_t mPz = 0;     //!< Z-Momentum
  Double_t mE  = 0;     //!< Energy
  Double_t mEta = -1;   //!< Pseudorapidity
  Double_t eta()  { if( mEta<0 ){ return asinh(mPz/pt()); }else{ return mEta; } }
  Double_t phi()  { return atan2(mPy,mPx); }
  Double_t pt()   { return sqrt( mPx*mPx + mPy*mPy ); }
  Double_t ptot() { return sqrt( mPx*mPx + mPy*mPy + mPz*mPz ); }
  Double_t theta(){ return 2.0*atan(exp(-1.0*eta())); }
  Double_t mass() { return sqrt(mE*mE - ptot()*ptot()); }
  //Need to project using momentum
  //Double_t mStarX = 0;  //!< Global STAR x postion
  //Double_t mStarY = 0;  //!< Global STAR y postion
  //Double_t mStarZ = 0;  //!< Global STAR z postion
  Double_t mAlpha = 0;  //!< Opening angle of pi0

  Double_t mE1 = 0;    //!< Energy of particle 1
  Double_t mX1 = 0;    //!< Lab frame x-position of particle 1, taking into account the beamline and z-vertex
  Double_t mY1 = 0;    //!< Lab frame y-position of particle 1, taking into account the beamline and z-vertex
  Double_t mZ1 = 0;    //!< Lab frame z-position of particle 1, taking into account the beamline and z-vertex
  Double_t mE2 = 0;    //!< Energy of particle 2
  Double_t mX2 = 0;    //!< Lab frame x-position of particle 2, taking into account the beamline and z-vertex
  Double_t mY2 = 0;    //!< Lab frame y-position of particle 2, taking into account the beamline and z-vertex
  Double_t mZ2 = 0;    //!< Lab frame z-position of particle 2, taking into account the beamline and z-vertex
  //TLorentzVector LV_P1(){ TLorentzVector v; v.SetPxPyPzE(mPX1,mPY1,mPZ1,mE1); return v; }
  //TLorentzVector LV_P2(){ TLorentzVector v; v.SetPxPyPzE(mPX2,mPY2,mPZ2,mE2); return v; }
  Double_t zgg {return fabs(mE1-mE2)/(mE1+mE2);}    //!< Energy asymmetry of pi0
  Double_t mDgg = 0;        //!< distance between the two particles (cm)
  Double_t mTpcVz = -999;   //!< TPC z vertex
  Double_t mVpdVz = -999;   //!< VPD z vertex
  Double_t mBbcVz = -999;   //!< BBC z Vertex

  ClassDef( FcsPi0Info, 1 );
}

class StMuFcsPi0TreeMaker : public StMaker {
public:

  StMuFcsPi0TreeMaker(const Char_t* name = "MuFcsPi0Maker");
  ~StMuFcsPi0TreeMaker();
  virtual Int_t Init();
  virtual Int_t InitRun(int runnumber);
  virtual Int_t Make();
  virtual Int_t Finish();
  void setOutFileName(const char* name) { mFilename = name; }
  #ifndef __CINT__
  template<typename T, typename... Args>
  void SetTrigs(int trigid, Args... restargs){ mTargetTrig.push_back(trigid); return SetTrigs(restargs...); } //function to set trigger ids to use. Does not check for repetition so need to be a good user
  #endif
  void IgnoreTrig(bool value=true){ mIgnoreTrig = value; }
  static void LoadDataFromFile(TFile* file, TTree* tree, TClonesArray* arr, TH1* hist=0);
  
protected:
  //UInt_t mEvent;  //Keeps track of number of events
  //int mTrig = -1;   //Found trigger index in vector 'mTargetTrig'
  //int mXing = 0;    //Bunch Crossing Id
  StMuDstMaker* mMuDstMkr = 0;
  StMuDst* mMuDst = 0;
  StMuEvent* mMuEvent = 0;
  const StTriggerData* mTrigData = 0;
  StRunInfo* mRunInfo = 0;

  StFcsDb* mFcsDb = 0;
  StMuFcsCollection* mMuFcsColl = 0;
  //StSpinDbMaker* mSpinDbMkr;
  StEpdGeom* mEpdGeo=0;
  
  std::vector<int> mTargetTrig;  //For Target Trigger ID
  bool mIgnoreTrig;  //!< flag to check if ignoring triggers or not
  //bool mReadMuDst;   //!< flag to check if reading from mudst or not (This can be used to turn off populating event info)
  //bool mReadSim;     //!< flag to check if reading mudst from simulations

  //Data to save
  TString mFilename = "";
  TFILE* mFile_Output = 0; //!< TFile to save all the data
  TTree* mPi0Tree = 0; //An internal tree that you can use to add any desired branches.  It will get written to the output file if it exists.
  TClonesArray* mPi0Arr = 0; //!< Array of #FcsPi0Info to store the valid pi0 events for analysis  
  TH1* mH1F_Entries = 0;           //!< Number of events processed no cuts (i.e. "Make" calls)
  
  // int mFilter = 0;
  // int mNEvents = -1;
  // int mNAccepted = 0;
  // int mMaxEvents = 30000;
  // int bins = 150;
  // float m_low = 0;
  // float m_up = 0.4;
  // float E_up = 10;
  // float E_low = 8;
  // float E_min = 1;

  ClassDef(StMuFcsPi0TreeMaker, 1);
};

#endif
