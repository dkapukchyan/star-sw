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

  @[September 26, 2024] > Couldn't get the variable argument function to work so just made a singular #AddTrig() function. Added accessor functions for the internal #TTree and #TClonesArrays. Implemented #LoadDataFromFile() in order to open and fill this object's data members with the data from the file. Added a #Print() function to print the information in #mPi0Tree. Got rid of the #TClonesArray for #FcsEventInfo and instead created a single pointer #mEvtInfo; because all I need is a single instance so I made just made a single branch for #FcsEventInfo pointing to #mEvtInfo in #mPi0Tree and related changes. Properly reseting #mNTrig in event if trigger was not found in #TargetTrig.

  @[September 27, 2024] > Added a function #ProjectToEpd() that will project an x,y,z position on the fcs and a zvertex onto the EPD plane. Used this function to fill FcsPhotonCandidate::mEpdHitNmip Changed found vertex from the hex to integer representation.

  @[October 4, 2024] > Added #mEpdNmipCut variable to make it easier to vary the EPD nmip cut. Added and implemented #mHists, a #HistManager, to manage all the histograms from 'AnaPi0Tree.cc' where the histograms were being filled by reading the tree. It was moved here to speed up processing the data. Implemented #LoadHists() to load all the histograms from a file and the tree. Added various "Paint" functions for the histograms. Because of this change the version was upgraded to version 2.
  
  @[October 8, 2024] > Found spikes in the point energy distribution and implemented #mH2F_PhotonHeatMapG and #mH2F_PhotonHeatMapB to see if the energy spikes are happening in a particular region. Also implemented #PaintEnergyZoom() to zoom in on the energy region and to plot the "G" and "B" histograms. Also increased the energy range to 200 GeV since the tower maximum was designed to go to 180 GeV.
  + @[October 10, 2024] > The data does show that there is a hot spot in the energy spike region when compared to a region without a spike. I need to analyze these histograms and the hit distribution histograms on a run by run basis.

  @[October 10, 2024] > Added Pi0 Pt histograms. Changed number of bins on the invariant mass histograms for "EpdCh" to be half that of "Best" and "EpdPh". Added #mH2F_EpdProjHitMap_Vcut that is the same as #mH2F_EpdProjHitMap except that it only gets filled when |vertex|<150cm. This was done to check if the bad distributions were coming from bad vertex projections.
  
  @[October 11, 2024] > It was true that bad vertex projections are causing the bad EPD projections and it was also caused by analyzing certain run numbers. Small fixes and some changes to plotting. Also, now #FcsPhotonCandidate array is properly sorted in descending order of energy but accompanying fix in "StMuFcsPi0Data". It was previously being sorted in ascending order of energy.
  
  @[October 21, 2024] > Added an additional loop over EPD tiles to check if a photon candidate actually intersects with an EPD tile before checking EPD hits. Modified the check condition to make sure an intersection with an EPD tile occurred before checking nmip. Also changed histogram binning and x range of #mH1F_EpdChInvMass and #mH1F_EpdChAllPoints to match #mH1F_BestPi0InvMass and #mH1F_EpdPhInvMass; also related changes for plotting. Also, added the unnormalized highest energy pair invariant mass histogram to the overlap plots

  @[November 1, 2024] > Added #mTreeOnBitMap to allow easily setting on and off certain branches of the Pi0 Tree or to even completely turn it off. Implemented the corresponding setter and accessor functions for the this bit map. Also, changed how the file is loaded so that it checks for and excludes missing branches when loading from a file.

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
#include "TLegend.h"

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
  void setEpdNmipCut(Double_t val ){ mEpdNmipCut = val; }
  void setRandomSeed(ULong_t seed){ mSpinRndm.SetSeed(seed); }
  UInt_t getRandomSeed(){ return mSpinRndm.GetSeed(); }
  void setHistManager( HistManager* hm );

  void setTreeOnBit(UShort_t bitmap){ mTreeOnBitMap = bitmap; }
  void setEventBit(bool val=1);
  void setPhotonOn(bool val=1);
  void setPi0On(bool val=1);

  UShort_t checkTreeOnBit() const { return mTreeOnBitMap; }
  bool isEventOn() const;
  bool isPhotonOn() const;
  bool isPi0On() const;
  
  TTree* getPi0Tree()const{ return mPi0Tree; }

  FcsEventInfo* getEvtInfo()const{ return mEvtInfo; }
  Int_t getTreeEntries()const{ return mPi0Tree->GetEntriesFast(); }
  Int_t getEntry(Int_t ientry){ return mPi0Tree->GetEntry(ientry); }
  Int_t getNTrig()const{ return mNTrig; }
  const Int_t* getTrig()const{ return mTriggers; }
  Int_t getTrig(Int_t itrig)const{ return mTriggers[itrig]; }

  TClonesArray* getPhArr()const{ return mPhArr; }
  Int_t getNPhoton()const{ return mPhArr->GetEntriesFast(); }
  FcsPhotonCandidate* getPhoton(Int_t iph)const{ return dynamic_cast<FcsPhotonCandidate*>(mPhArr->UncheckedAt(iph)); }

  TClonesArray* getPi0Arr()const{ return mPi0Arr; }
  Int_t getNPi0()const{ return mPi0Arr->GetEntriesFast(); }
  FcsPi0Candidate* getPi0(Int_t ipi0)const{ return dynamic_cast<FcsPi0Candidate*>(mPi0Arr->UncheckedAt(ipi0)); }
  //#ifndef __CINT__
  //void SetTrigs(const char* trigname,...);//{ mTargetTrig.emplace_back(trigname); }
  //template<typename... Args>
  //void SetTrigs(Args... restargs){ SetTrigs(restargs...); } //function to set trigger ids to use. Does not check for repetition so need to be a good user
  //#endif
  void AddTrig(const char* trigname ){ mTargetTrig.emplace_back(trigname); }  //function to set trigger ids to use. Does not check for repetition so need to be a good user
  void IgnoreTrig(bool value=true){ mIgnoreTrig = value; }
  void LoadDataFromFile(TFile* file);//, TTree&* tree, FcsEventInfo&* evt,Int_t& ntrig, Int_t&* triggers,  TClonesArray&* pharr, TClonesArray&* pi0arr, TH1&* hist=0):
  virtual UInt_t LoadHists(TFile* file);

  virtual void Print(Option_t* opt="") const; //"e" for event, "t" for trigger, "g" for photon, "p" for pi0, "a" for all

  static std::vector<Double_t> ProjectToEpd(Double_t xfcs, Double_t yfcs, Double_t zfcs, Double_t zvertex);

  void PaintEventQa(TCanvas* canv,  const char* savename="testevent.png")    const;
  void PaintPhotonQa(TCanvas* canv, const char* savename="testphoton.png")   const;
  void PaintBestPi0(TCanvas* canv,  const char* savename="testbestpi0.png")  const;
  void PaintEpdPhPi0(TCanvas* canv, const char* savename="testepdphpi0.png") const;
  void PaintEpdChPi0(TCanvas* canv, const char* savename="testepdchpi0.png") const;
  void PaintPi0Overlap(TCanvas* canv, const char* savename = "testpi0overlap.png") const;
  void PaintEnergyZoom(TCanvas* canv, const char* savename = "testenergyzoom.png") const;

  static void AddHistStatsOneline( TLegend* HistLeg, const TH1* h1, const std::string &title="" );
  
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
  FcsEventInfo* mEvtInfo = 0;           ///< #FcsEventInfo object for TTree
  //For trigger branch of tree 
  static const UShort_t mMaxTrigs = 65; ///< 64 FCS triggers + 1 for any other not found
  Int_t mNTrig  = 0;                    ///< Total triggers in the event
  Int_t mTriggers[mMaxTrigs];           ///< Array of Triggers in the event 
  //void ResetTrigs();                    ///< Reset #mTriggers to default values
  TClonesArray* mPhArr = 0;             ///< #TClonesArray of all #FcsPhotonCandidate
  //TClonesArray* mBestPharr = 0;         ///< #TClonesArray of #FcsPhotonCandidates which are highest energy pairs
  TClonesArray* mPi0Arr = 0;            ///< Array of #FcsPi0Candidate to store the best pi0 candidates for analysis
  
  TH1* mH1F_Entries = 0;                ///< Number of events processed no cuts (i.e. "Make" calls)
  TH1* mH1F_Triggers = 0;               ///< Triggers used in analysis
  TH1* mH2F_foundVvertex = 0;           ///< found vertex bit vs. Vertex

  TH1* mH2F_PhotonHeatMap = 0;          ///< Distribution of photons in STAR x,y space
  TH1* mH2F_PhotonHeatMapG = 0;         ///< Distribution of photons in STAR x,y space when energy has a specific value near an energy spike
  TH1* mH2F_PhotonHeatMapB = 0;         ///< Distribution of photons in STAR x,y space when energy has a specific value on an energy spike
  TH1* mH2F_EpdProjHitMap = 0;          ///< Distribution of x,y projections of photon candidates onto STAR EPD plane in x,y space
  TH1* mH2F_EpdProjHitMap_Vcut = 0;     ///< Distribution of x,y projections of photon candidates onto STAR EPD plane in x,y space with cut |vertex|<150cm
  TH1* mH2F_EpdNmip = 0;                ///< Nmip distributions for EPD matched projected clusters (x-axis bin 1) and points (y-axis bin 2)

  TH1* mH1F_ClusterEnergy = 0;          ///< All Cluster energy
  TH1* mH1F_PointEnergy = 0;            ///< All Point energy
  TH1* mH2F_Energy_ph1Vph2 = 0;         ///< Histogram of two photons energy used in reconstruction
  
  TH1* mH1F_BestPi0Mass = 0;            ///< Invariant mass with just highest energy points
  //TH1* mH2F_Pi0HeatMap = 0;             ///< x,y locations of BestPi0
  TH1* mH1F_PointMult = 0;              ///< Raw point multiplicity in event
  TH1* mH1F_BestPi0Zgg = 0;             ///< Zgg of for BestPi0
  TH1* mH1F_BestPi0Phi = 0;             ///< Azimuthal angle for BestPi0
  TH1* mH1F_BestPi0Eta = 0;             ///< Psuedorapidity for BestPi0
  TH1* mH1F_BestPi0En = 0;              ///< Energy for BestPi0
  TH1* mH1F_BestPi0Pt = 0;              ///< Pt for BestPi0
  TH1* mH1F_AllPointPairMass = 0;       ///< Invariant mass of all pairs of points

  TH1* mH1F_EpdPhInvMass = 0;           ///< Invariant mass with EPD nmip cut to isolate uncharged particles (EpdPh stands for Epd Photon, as in uncharged particle)
  //TH1* mH2F_EpdPhHeatMap = 0;           ///< x,y locations for EpdPh
  TH1* mH1F_EpdPhPointMult = 0;         ///< Point Multiplicity for EpdPh
  TH1* mH1F_EpdPhZgg = 0;               ///< Zgg of for EpdPh
  TH1* mH1F_EpdPhPhi = 0;               ///< Azimuthal angle for EpdPh
  TH1* mH1F_EpdPhEta = 0;               ///< Psuedorapidity for EpdPh
  TH1* mH1F_EpdPhEn = 0;                ///< Energy for EpdPh
  TH1* mH1F_EpdPhPt = 0;                ///< Pt for EpdPh
  TH1* mH1F_EpdPhAllPoints = 0;         ///< Invariant mass of all point pairs with EPD nmip cut to isolate uncharged particles

  TH1* mH1F_EpdChInvMass = 0;           ///< Invariant mass with EPD nmip cut to isolate charged particles (EpdCh stands for Epd Charged, as in charged particle)
  //TH1* mH2F_EpdChHeatMap = 0;           ///< x,y locations for EpdCh
  TH1* mH1F_EpdChPointMult = 0;         ///< Point multiplicty for EpdCh
  TH1* mH1F_EpdChZgg = 0;               ///< Zgg for EpdCh
  TH1* mH1F_EpdChPhi = 0;               ///< Azimuthal angle for EpdCh
  TH1* mH1F_EpdChEta = 0;               ///< Psuedorapidity for EpdCh
  TH1* mH1F_EpdChEn = 0;                ///< Energy for EpdCh
  TH1* mH1F_EpdChPt = 0;                ///< Pt for EpdCh
  TH1* mH1F_EpdChAllPoints = 0;         ///< Invariant mass of all point pairs with EPD nmip cut to isolate charged particles  

  //Also separate low point multiplicity events
  //TH1* mH1F_2ndBestPi0Mass = 0;
  //TH1* mH1F_3rdBestPi0Mass = 0;
  //TH1* mH1F_LowPointMult = 0;
  
  Double_t mEnCut = 1;                  ///< Energy Cut for #FcsPhotonCandidates
  Double_t mEpdNmipCut = 0.7;           ///< Cut on EPD nmip to classify cluster or point as charged or uncharged
  UShort_t mTreeOnBitMap = 0x7;         /// Turn on or off branches in the pi0 tree. first bit is events, second bit is photon branch, third bit is pi0 branch. Turn on all branches by default
  
private:
  HistManager* mHists = 0;            ///< Manage loading and saving histograms
  bool mInternalHists = false;        ///< Boolean to keep track if mHists was added externally or an internal one was created
  TRandom3 mSpinRndm;                 ///< Spin state randomizer

  ClassDef(StMuFcsPi0TreeMaker, 2)
};

#endif
