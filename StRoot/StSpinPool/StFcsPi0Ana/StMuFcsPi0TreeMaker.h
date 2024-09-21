/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpose of this class is to generate a #TTree of photon and pi0 candidates to use in a transverse single spin assymmetry analysis. It will also store basic event level information as well

  DESCRIPTION
  This Maker will loop over all the clusters and points in the MuDst tree and fill a new #TTree with clusters and points as photon candidates (#FcsPhotonCandidate). It will loop over those candidates and pair them together to generate a list of pi0 candidates (#FcsPi0Candidate). It only combines clusters with clusters and points with points. It will only apply an energy cut and trigger selection to filter out what events get saved. Trigger Filter and energy cut can be turned off with appropriate function call. Also stores event level information (#FcsEventInfo) and trigger information (#mTriggers).

  LOG
  @[March 5, 2024] > Copied from *StMuFcsPi0Maker* and modified to write pi0s to a tree

  @[May 24, 2024] > Added #FcsPi0Info class to hold information about the reconstructed Pi0s. Got rid of the histograms.

  @[September 19, 2024] > Major revisions to write the code in the DESCRIPTION. Added two new classes; #FcsEventInfo which will hold event level information like run numbers, fill numbers, etc; #FcsPhotonCandidate which holds basic cluster and point information like the lorentz vector and the x,y,z information. Changed the class name #FcsPi0Info to #FcsPi0Candidate. #FcsPi0Candidate holds the index of the #FcsPhotonCandidate that was used to create the Pi0. It stores basic information like the lorentz vector, mass, eta, alpha, and dgg. Added ability to generate a random spin bit. Added energy cut option for what clusters and points the tree will be filled with. Added members to utilize EPD information from #StMuEpdRun22QaMaker and other EPD related classes. Added #mEvtInfo and #mPhArr as #TClonesArray to add as branches to the analysis #TTree #mPi0Tree.

  @[September 21, 2024] > Moved #FcsEventInfo, #FcsPhotonCandidate, and #FcsPi0Candidate classes to their own file in *StSpinPool/StFcsTreeManager* so that I can load that library without having to load other STAR libraries; I would have to load the STAR libraries If I had to load this one instead. Mainly, this simplifies the rootmap file I would need write. Fixed how points are being used to make Pi0s. Also changed comments to ROOT friendly style.

*/


#ifndef STMUFCSPI0TREEMAKER_HH
#define STMUFCSPI0TREEMAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1F.h"

//STAR Headers
#include "StEnumerations.h"
#include "StMaker.h"
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
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

#include "StSpinPool/StFcsTreeManager/StMuFcsPi0Data.h"
#include "StMuEpdRun22QaMaker.h"
#include "StFcsRun22TriggerMap.h"

class StEpdGeom;

class StMuFcsPi0TreeMaker : public StMaker {
public:

  StMuFcsPi0TreeMaker(const Char_t* name = "MuFcsPi0Maker");
  ~StMuFcsPi0TreeMaker();
  virtual Int_t Init();
  virtual Int_t InitRun(int runnumber);
  virtual Int_t Make();
  virtual Int_t Finish();
  void setOutFileName(const char* name) { mFilename = name; }
  void setEnergyCut(Double_t val){ mEnCut = val; }
  void setRandomSeed(ULong_t seed){ mSpinRndm.SetSeed(seed); }
  UInt_t getRandomSeed(){ return mSpinRndm.GetSeed(); }
  #ifndef __CINT__
  template<typename T, typename... Args>
  void SetTrigs(const char* trigname, Args... restargs){ mTargetTrig.emplace_back(trigname); return SetTrigs(restargs...); } //function to set trigger ids to use. Does not check for repetition so need to be a good user
  #endif
  void IgnoreTrig(bool value=true){ mIgnoreTrig = value; }
  static void LoadDataFromFile(TFile* file, TTree* tree, TClonesArray* arr, TH1* hist=0);
  
protected:  
  StMuDstMaker* mMuDstMkr        = 0;
  StMuDst* mMuDst                = 0;
  StMuEvent* mMuEvent            = 0;
  const StTriggerData* mTrigData = 0;
  StRunInfo* mRunInfo            = 0;
  StSpinDbMaker* mSpinDbMkr      = 0;

  StFcsDb* mFcsDb                   = 0;
  StMuFcsCollection* mMuFcsColl     = 0;
  StFcsRun22TriggerMap* mFcsTrigMap = 0;

  StEpdGeom* mEpdGeo             = 0;
  TClonesArray* mMuEpdHits       = 0;
  StEpdHitMaker* mEpdHitMkr      = 0;
  StEpdCollection* mEpdColl      = 0;
  StMuEpdRun22QaMaker* mEpdQaMkr = 0;
  
  std::vector<std::string> mTargetTrig;  //For Target Trigger ID
  bool mIgnoreTrig = false;  ///< flag to check if ignoring triggers or not
  //bool mReadMuDst;   ///< flag to check if reading from mudst or not (This can be used to turn off populating event info)
  //bool mReadSim;     ///< flag to check if reading mudst from simulations

  //Data to save
  TString mFilename = "";
  TFile* mFile_Output = 0;              ///< #TFile to save all the data
  TTree* mPi0Tree = 0;                  ///< #TTree with desired branches
  TClonesArray* mEvtInfo = 0;           ///< #TClonesArray of #FcsEventInfo
  //For trigger branch of tree 
  static const UShort_t mMaxTrigs = 65; ///< 64 FCS triggers + 1 for any other not found
  Int_t mNTrig  = 0;                    ///< Total triggers in the event
  Int_t mTriggers[mMaxTrigs];           ///< Array of Triggers in the event 
  //void ResetTrigs();                    ///< Reset #mTriggers to default values
  TClonesArray* mPhArr = 0;             ///< #TClonesArray of #FcsPhotonCandidate
  TClonesArray* mPi0Arr = 0;            ///< Array of #FcsPi0Candidate to store the valid pi0 events for analysis  
  TH1* mH1F_Entries = 0;                ///< Number of events processed no cuts (i.e. "Make" calls)

  Double_t mEnCut = 1;                  ///< Energy Cut for #FcsPhotonCandidates
  
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

private:
  TRandom3 mSpinRndm;

  ClassDef(StMuFcsPi0TreeMaker, 1)
};

#endif
