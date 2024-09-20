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

#include "StMuEpdRun22QaMaker.h"
#include "StFcsRun22TriggerMap.h"

class StEpdGeom;

class FcsEventInfo : public TObject
{
public:
  FcsEventInfo();
  ~FcsEventInfo();

  Int_t mRunTime = -1;       //!< Time of the run
  Int_t mRunNum = -1;        //!< Run number for event
  UInt_t mFill = 0;          //!< Fill number for event
  UInt_t mEvent = -1;        //!< STAR Event Id
  Int_t mBx48Id = -1;        //!< 48 bit bunch Id for event
  Int_t mBx7Id = -1;         //!< 7 bit bunch Id for event
  UShort_t mSpin = 0;        //!< [Spin bit, this is source polarization](https://drupal.star.bnl.gov/STAR/blog/oleg/spin-patterns-and-polarization-direction)
  Short_t spinFrom4BitSpin(); //!< Correctly accounts for the spin flip when working with STAR data

  Int_t mTofMultiplicity = -1; //!< TOF Multiplicity
  
  Double_t mVpdVz = -999;      //!< VPD z Vertex
  Double_t mBbcVz = -999;      //!< BBC z Vertex
  Double_t mBbcTacDiff = 0; //!< BBC TAC difference
  Double_t mEpdTacEarlyW = 0;  //!< Earliest EPD TAC for West with cuts 1<adcnmip<15 && TAC>50
  Double_t mEpdTacEarlyE = 0;  //!< Earliest EPD TAC for East with cuts 1<adcnmip<15 && TAC>50
  Double_t mEpdAvgW = 0;    //!< Average EPD TAC for West with cuts 1<adcnmip<15 && TAC>50
  Double_t mEpdAvgE = 0;    //!< Average EPD TAC for East with cuts 1<adcnmip<15 && TAC>50
  //Double_t EpdTacDiffEarly();
  //Double_t EpdTacDiffAvg();
  Double_t mEpdVz = -999;      //!< EPD z Vertex
  Double_t mZdcVz = -999;      //!< ZDC z Vertex
  Short_t mFoundVertex = 0;    //!< Bit vector encoding for which vertex was best; 0 means no vertex, 1=Vpd,2=Epd,4=Bbc

  //This will be used to indicate how many clusters are in the #TClonesArray of #FcsPhotonCandidate. Everything from this number to the size of the array will be points for a given detector Id. I did it this way so I don't have to create a separate branch holding these two numbers and there should only be one #FcsEventInfo object. Also, didn't want a seperate class for clusters and points since they will store the same information. It is kind of a hack since I know that I am only looping up to detector id 2.
  Int_t mClusterSize = 0;       //!< Size of clusters in #mPhArr in #StMuFcsPi0TreeMaker. This means 0 to <#mClusterSize is cluster photon candidates
  //Int_t mPointSizeDet0   = 0;   //!< Size of points for detectorid 0.This means from (including) #mClusterSizeDet0 to <#mPointSizeDet0 is points for detector id 0
  //Int_t mClusterSizeDet1 = 0;   //!< Size of clusters for detector id 1 (Ecal South). This means from (including) #mPointSizeDet0 to <#ClusterSizeDet1 is clusters for detector id 1.

  ClassDef( FcsEventInfo, 1 )
};

//Class to hold basic particle info from which pi0s can be reconstructed
class FcsPhotonCandidate : public TObject
{
public:
  FcsPhotonCandidate();
  ~FcsPhotonCandidate();

  bool mFromCluster = false;  //!< True if from an FCS cluster
  Short_t mDetId = -1;        //!< Detector Id where candidate was found

  Double_t mX = 0;           //!< STAR global x coordinate
  Double_t mY = 0;           //!< STAR global y coordinate
  Double_t mZ = 0;           //!< STAR global z coordinate

  Double_t mEn = 0;          //!< Energy
  Double_t mPxRaw = 0;       //!< X momentum assuming 0,0,0 vertex
  Double_t mPyRaw = 0;       //!< Y momentum assuming 0,0,0 vertex
  Double_t mPzRaw = 0;       //!< Z momentum assuming 0,0,0 vertex

  Double_t mPxVert = 0;       //!< X momentum using best found vertex
  Double_t mPyVert = 0;       //!< Y momentum using best found vertex
  Double_t mPzVert = 0;       //!< Z momentum using best found vertex

  Double_t mEpdHitNmip = 0;   //!< NMIP value from EPD hit

  TLorentzVector lvRaw();        //!< TLorentz vector for this condidate with 0,0,0 vertex momentum
  TLorentzVector lvVert();       //!< TLorentz vector for this candidate with vertex momentum
  Double_t magPosition();        //!< Magnitude of postiion vector i.e. sqrt(#mX^2+#mY^2+#mZ^2)

  ClassDef( FcsPhotonCandidate, 1 )
};

//Class to hold basic info for reconstructed pi0 candidates
class FcsPi0Candidate : public TObject
{
public:
  FcsPi0Candidate();
  ~FcsPi0Candidate();

  bool mFromCluster = false;  //!< Pi0 reconstructed from clusters or points
  Int_t mPhoton1Idx = -1;     //!< Index in #TClonesArray of #FcsPhotonCandidate 1 that was used to reconstruct this Pi0
  Int_t mPhoton2Idx = -1;     //!< Index in #TClonesArray of #FcsPhotonCandidate 2 that was used to reconstruct this Pi0

  //#StMuFcsPi0TreeMaker will only store information for the Lorentz vector, and other data from the best vertex. To switch to 0,0,0 vertex; use the photon index
  Double_t mPx = 0;            //!< X-Momentum from Lorentz vector of two reconstructed candidates
  Double_t mPy = 0;            //!< Y-Momentum from Lorentz vector of two reconstructed candidates
  Double_t mPz = 0;            //!< Z-Momentum from Lorentz vector of two reconstructed candidates
  Double_t mEn = 0;            //!< Energy from Lorentz vector of two reconstructed candidates
  
  Double_t mEta = -1;          //!< Pseudorapidity from the lorentz vector
  Double_t eta();              //!< Needed since in simulations the stored eta was not there when reconstructing from simulated photons
  Double_t phi();              //!< Angle in x,y plane
  Double_t pt();               //!< Transverse momentum
  Double_t ptot();             //!< Total momentum
  Double_t theta();            //!< azimuthal angle (angle from z-axis to x-y plane)
  Double_t mass();             //!< Invariant mass of the Pi0
  //Need to project using momentum
  //Double_t mStarX = 0;     //!< Global STAR x postion from best vertex
  //Double_t mStarY = 0;     //!< Global STAR y postion from best vertex
  //Double_t mStarZ = 0;     //!< Global STAR z postion from best vertex
  Double_t mDgg = 0;        //!< Euclidean Distance between the two particles (cm) 
  Double_t mZgg = 0;       //!< Energy Asymmetry |E1-E2|/(E1+E2) of pi0
  Double_t mAlpha = 0;     //!< Opening angle of pi0
  Double_t mInvMass = -1;  //!< invariant mass using best vertex as a variable to make it easier for analysis
  
  static Double_t zgg(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2);        //!< Energy asymmetry of pi0
  static Double_t dgg(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2);        //!< Distance between the two particles (cm)
  static Double_t alpha(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2);      //!< Opening angle of the two photons

  ClassDef( FcsPi0Candidate, 1 )
};

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
  bool mIgnoreTrig = false;  //!< flag to check if ignoring triggers or not
  //bool mReadMuDst;   //!< flag to check if reading from mudst or not (This can be used to turn off populating event info)
  //bool mReadSim;     //!< flag to check if reading mudst from simulations

  //Data to save
  TString mFilename = "";
  TFile* mFile_Output = 0;              //!< #TFile to save all the data
  TTree* mPi0Tree = 0;                  //!< #TTree with desired branches
  TClonesArray* mEvtInfo = 0;           //!< #TClonesArray of #FcsEventInfo
  //For trigger branch of tree 
  static const UShort_t mMaxTrigs = 65; //!< 64 FCS triggers + 1 for any other not found
  Int_t mNTrig  = 0;                    //!< Total triggers in the event
  Int_t mTriggers[mMaxTrigs];           //!< Array of Triggers in the event 
  //void ResetTrigs();                    //!< Reset #mTriggers to default values
  TClonesArray* mPhArr = 0;             //!< #TClonesArray of #FcsPhotonCandidate
  TClonesArray* mPi0Arr = 0;            //!< Array of #FcsPi0Candidate to store the valid pi0 events for analysis  
  TH1* mH1F_Entries = 0;                //!< Number of events processed no cuts (i.e. "Make" calls)

  Double_t mEnCut = 1;                  //!< Energy Cut for #FcsPhotonCandidates
  
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
