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

  @[November 26, 2024] > Added an array of histograms #mH1F_InvMassEpdCuts and related code to check how different EPD nmip values are changing the pi0 signal region. also added #NEPDCUTS to make filling and looping over array easier. Implemented several cuts for the pi0 A_N analysis and added an array of histograms #mH1F_NPi0ByEnByPhi to keep track of the number of pi0s found in each energy bin and phi bin for a given spin configuration; it is not "filled" but populated after counting the pi0s from the different spin states. In this way I no longer need to loop over a tree of pi0s after processing MuDsts but can just merge the histograms to get the total number of pi0s for computing A_N. To this effect added static integers #NENERGYBIN and #NPHIBIN that keeps track of the energy and phi bins I will use. Added #mH1F_InvMassAllCuts to see the invariant mass distribution of the found pi0s after all the cuts and also the #mH1F_Pi0MultAllCuts to see how many potential pi0s are left in each event after all the cuts. I thought this would be 1 pi0 but turns out it is quite significant. It's mostly coming from one high energy "photon" that gets paired with other "photons". However, this may be a bit biased since I cut out other combinations for the sake of making the trees take up less space. Also, in #Make() added #emtrigfound which is a boolean to indicate one of the EM triggers had fired. This is used to separate the EM trigger events from others and is one of the cuts applied. The cuts are in order of checking in the cod: emtrig,|vertex|<100,fcspoint,Zgg<0.7,both points of pi0 passed epd nmip cut, pi0 Pt larger than trigger threshold. Changed code to grab 4 bit spin from bx7. Also wrote paint functions for the new histograms.

  @[December 20, 2025] > Added several histograms related to filling invariant mass with cuts to isolate good pi0s and making histograms of number of pi0s, along with a few QA other histograms. Implemented drawing functions for them as well.

  @[January 8, 2025] > Fixed how histograms are loaded so it can work with multiple files. Added two graphs and methods for processing, filling, and plotting them to look at data over many runs; one for the peak of the invariant mass, the other for pi0 energy.

  @[January 15, 2025] > Made #NENERGYBIN and #NPHIBIN public so I can use it outside the class. Added #MergeForTssa() that will also merge the NPi0 histograms so that I don't need to hadd after running the QA. In order for the code to work properly I had to fix the second argument in the array. I didn't seem to work with NPHIBIN so I need to be careful moving forward.

  @[February 3, 2025] > Changed the energy and phi binning. Changed the "NPi0" histogram arrays of phi and energy to 2D histograms where the histogram now bins the phi and energy bins rather than by an array. It is still an array but now the array represents different beam and spin states. Also, changed the invariant mass histogram with all the cuts from an array of different energies and phi values to a 3d histogram where the x and y bins reprsent the phi and energy bins respectively. Changed #MergeTssa() to merge the new 2d and 3d histograms of #mH2F_NPi0_enVphi and #mH3F_AllCutsInvMass_enVphi. Added static function #DoTssaAna() which will take as input a 2x2 array of 2d histograms (matching the form of #mH2F_NPi0_enVphi) and compute the transverse single spin asymmetry for each energy bin and beam species generating a new array of 1D histograms that are the phi and raw asymmetry. Modified the draw functions to compensate for the changed histograms.

  @[March 17, 2025] > Finalizedish code for dissertation with all neccessary functionality added to get a  result including calculating luminosities, polarizations, and signal and background fits and subtraction. Modified #LoadDataFromFile() to be able to handle more histograms. Added more painting functions for the different analysis functionality and other histograms that was added. MergeTssa(), DoTssaAna(), DoTssaFit(), DoPi0Fits(), DoBgCorrectedAn(), ReadPolFile() are all functions that were added and modified to allow the A_N calculation to happen in both pre and post data processing. Added histograms for the polarization. Added more invariant mass histograms at different cut criteria. Added "NPi0" histograms for the signal region and two background regions. Added struct #PolData which is filled with the polarization data from the RHIC polarimeter website (and is saved elsewhere specially formatted to be read by this class)

  @[June 29, 2025] > Made a hack so that rather than energy it is binned in xf. The names and symbols still say energy but it is really xf. Added Zgg and Dgg histograms for the pi0 after all cuts.

  @[July 31, 2025] > Properly changed all the energy naming to x_{Feynmann} (x_F,xf) and changed slightly the xf binning ranges

  @[August 4, 2025] > Added histograms to help understand the data

  @[August 11, 2025] > Changed some histograms to be split by the 4 different EM triggers. Added code and histograms to look at pi0s where one of the points satisfy the photon or charged criteria from the EPD. Fixed bug where code wasn't properly reading in the last line of the polarization data file.

  @[August 12, 2025] > Added histograms and code to check for cases when single photon candidate passes the EPD nmip cut. Needs fixing since I check mFromPh==0 but only store the best pi0 so need to fix.

  @[August 23, 2025 > Fixed the cases for single photon candidate passing the cut by getting rid of the loops where I fill in a vector of "good" photons and electrons and now keep all pi0s in the array and loop over all those pi0s and only select the two photon cut when filling for the A_N analysis. This means that I got rid of the histograms related to the "Best" pi0 and the photon multiplicity with the EPD nMIP cut. Also, cleaned up extraneous code that is no longer needed with these changes. Added more paint functions and modified existing ones to paint the more relevant histograms. Added histograms to check how well the projected point from the FCS gives the correct matching EPD tile to check by computing various differences between the projected point and all other EPD hit tiles. Found that nMIP 0.4 has best valley so changed code to use nMIP=0.4 for the cut because it may be working better. Added #EpdTilePoly() as a precursor to drawing EPD tiles.

  @[September 26, 2025] > Implemented mixed events to compare how the EPD projections are working in real vs. mixed events. Implemented #Clear() function to use to clean up variables rather than at the end of #Make(). Implemented plotting functions for the mixed event histograms. Wrote #DrawEpdProjection() to work as an event display of projected FCS points and clusters on top of EPD hits. Looked at low point multiplicity events based on the event display. Implemented several functions to create polylines corresponding to different EPD adjacencies and also to draw the adjacencies.

*/


#ifndef STMUFCSPI0TREEMAKER_HH
#define STMUFCSPI0TREEMAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TString.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGeoPolygon.h"

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
  virtual void Clear(Option_t* option="");           ///< Gets called before #Make() in StChain::EventLoop()
  virtual Int_t Make();
  virtual Int_t Finish();
  void setOutFileName(const char* name) { mFilename = name; }
  void setEnergyCut(Double_t val){ mEnCut = val; }
  void setEpdNmipCut(Double_t val ){ mEpdNmipCut = val; }
  void setRandomSeed(ULong_t seed){ mSpinRndm.SetSeed(seed); }
  UInt_t getRandomSeed(){ return mSpinRndm.GetSeed(); }
  HistManager* getHists()const{ return mHists; }
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
  UInt_t LoadDataFromFile(TFile* file);//, TTree&* tree, FcsEventInfo&* evt,Int_t& ntrig, Int_t&* triggers,  TClonesArray&* pharr, TClonesArray&* pi0arr, TH1&* hist=0):
  virtual UInt_t LoadHists(TFile* file);

  //void AnalyzePi0s();

  virtual void Print(Option_t* opt="") const; //"e" for event, "t" for trigger, "g" for photon, "p" for pi0, "a" for all

  static std::vector<Double_t> ProjectToEpd(Double_t xfcs, Double_t yfcs, Double_t zfcs, Double_t zvertex);

  void PaintEventQa(TCanvas* canv,  const char* savename = "testevent.png")    const;
  void PaintPolarization(TCanvas* canv, const char* savename = "testpolarization.png") const;
  void PaintPhotonQa(TCanvas* canv, const char* savename = "testphoton.png")   const;
  void PaintAllPi0(TCanvas* canv,  const char* savename = "testallpi0.png")  const;
  void PaintNoEpdCut(TCanvas* canv,  const char* savename="testnoepdcutpi0.png")  const;
  void PaintEpdPhPi0(TCanvas* canv, const char* savename = "testepdphpi0.png") const;
  void PaintEpdChPi0(TCanvas* canv, const char* savename = "testepdchpi0.png") const;
  void PaintEpdSinglePh(TCanvas* canv, const char* savename = "testepdsingleph.png") const;
  void PaintEpdSingleCh(TCanvas* canv, const char* savename = "testepdsinglech.png") const;
  void PaintPi0Overlap(TCanvas* canv, const char* savename = "testpi0overlap.png") const;
  void PaintInvMassEpdQa(TCanvas* canv, const char* savename = "testinvmasscutqa.png") const;
  void PaintEpdAllDistQa(TCanvas* canv, const char* savename = "testepdalldistqa.png") const;
  void PaintEpdAllDistQaLowMult(TCanvas* canv, const char* savename = "testepdalldistqalowmult.png" ) const;
  void PaintEpdTileDistQa(TCanvas* canv, const char* savename = "testepdtiledistqa.png") const;
  void PaintEpdDistAnaQa(TCanvas* canv, const char* savename = "testepddistanaqa.png") const;
  void PaintEpdQa(TCanvas* canv, const char* savename = "testepdsingleqa.png") const;
  void PaintEnergyZoom(TCanvas* canv, const char* savename = "testenergyzoom.png") const;
  void PaintEpdNmipCuts(TCanvas* canv, const char* savename = "testepdnmipcut.png") const;
  void PaintPi0Cuts(TCanvas* canv, const char* savename = "testpi0cuts.png") const;
  void PaintInvMassCuts(TCanvas* canv, const char* savename = "testinvmasscuts.pdf") const;
  void PaintNpi0Inc(TCanvas* canv, const char* savename = "testnpi0inc.png") const;
  void PaintNpi0Bg1(TCanvas* canv, const char* savename = "testnpi0bg1.png") const;
  void PaintNpi0Bg2(TCanvas* canv, const char* savename = "testnpi0bg2.png") const;
  
  void PaintAllHistOneTrigger(TCanvas* canv, int trigidx, const char* savename) const;
  
  void PaintOneHistAllTrigger(TCanvas* canv, TObjArray* histarr, const char* drawoption, const char* savename) const;
  void PaintAllTrigInvMass(TCanvas* canv, const char* savename="testAllTrigInvMass.png") const{ PaintOneHistAllTrigger(canv,mH1F_InvMassAllCuts,"hist e",savename); }
  void PaintAllTrigPi0Mult(TCanvas* canv, const char* savename="testAllTrigPi0Mult.png") const{ PaintOneHistAllTrigger(canv,mH1F_Pi0MultAllCuts,"hist e",savename); }
  void PaintAllTrigxF(TCanvas* canv, const char* savename="testAllTrigxF.png") const{ PaintOneHistAllTrigger(canv,mH1F_AllCuts_xF,"hist e",savename); }
  void PaintAllTrigxFZoom(TCanvas* canv, const char* savename="testAllTrigxFZoom.png") const{ PaintOneHistAllTrigger(canv,mH1F_AllCuts_xFZoom,"hist e",savename); }
  void PaintAllTrigZgg(TCanvas* canv, const char* savename="testAllTrigZgg.png") const{ PaintOneHistAllTrigger(canv,mH1F_AllCuts_Zgg,"hist e",savename); }
  void PaintAllTrigDgg(TCanvas* canv, const char* savename="testAllTrigDgg.png") const{ PaintOneHistAllTrigger(canv,mH1F_AllCuts_Dgg,"hist e",savename); }
  void PaintAllTrigPi0En(TCanvas* canv, const char* savename="testAllTrigPi0En.png") const{ PaintOneHistAllTrigger(canv,mH1F_AllCuts_Pi0En,"hist e",savename); }
  void PaintAllTrigPi0massVen(TCanvas* canv, const char* savename="testAllTrigPi0massVen.png") const{ PaintOneHistAllTrigger(canv,mH2F_AllCuts_Pi0_massVen,"colz",savename); }
  void PaintAllTrigPi0xfVen(TCanvas* canv, const char* savename="testAllTrigPi0xfVen.png") const{ PaintOneHistAllTrigger(canv,mH2F_AllCuts_Pi0_xfVen,"colz",savename); }
  void PaintAllTrigPi0ptVeta(TCanvas* canv, const char* savename="testAllTrigPi0ptVeta.png") const{ PaintOneHistAllTrigger(canv,mH2F_AllCuts_Pi0_ptVeta,"colz",savename); }
  void PaintAllTrigPi0etaVphi(TCanvas* canv, const char* savename="testAllTrigPi0etaVphi.png") const{ PaintOneHistAllTrigger(canv,mH2F_AllCuts_Pi0_etaVphi,"colz",savename); }
  void PaintAllTrigPi0yVx(TCanvas* canv, const char* savename="testAllTrigPi0yVx.png") const{ PaintOneHistAllTrigger(canv,mH2F_AllCuts_Pi0_yVx,"colz",savename); }

  void PaintPhotonQaForDefense(TCanvas* canv, const char* savename)   const;
  void PaintPi0QaForDefense(TCanvas* canv, const char* savename) const;

  static void AddHistStatsOneline( TLegend* HistLeg, const TH1* h1, const std::string &title="" );

  Int_t LoadGraphsFromFile(TFile* file, TObjArray* graphs );
  void FillGraphs(Int_t irun);

  void DrawQaGraphs(TCanvas* canv, const char* savename="testGraphPi0.png");

  Int_t DrawEpdProjection(TCanvas* canvas, const char* savename);  //!< If pass all cuts and drawn then returns 1, failure returns 0

  void DrawSelectEpdAdjTiles(TCanvas* canv) const;
  void DrawAllEpdAdjTiles(TCanvas* canv, Int_t tt) const;

  //static const short NENERGYBIN = 8;     ///< Number of energy bins
  static const short NXFBIN = 8;         ///< Number of x_F bins (should match energy binning)
  static const short NPHIBIN = 24;       ///< Number of phi bins

  //Double_t ebins[NENERGYBIN+1] = {0, 15, 20, 25, 30, 40, 55, 70, 100};
  //static const Double_t xfbins[NXFBIN+1] = {0, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.5}; //@[June 28, 2025] > xF binning such that statistics are even across bins
  //static constexpr Double_t xfbins[NXFBIN+1] = {0, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.5}; //@[June 28, 2025] > Trying as static data member
  //const static double xfbins[NXFBIN+1];  //@[July 31, 2025] > Can't use constexpr with ROOT 5 but this is the right c++ 11 way to do it
  //Double_t xfbins[NXFBIN+1] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.3, 0.5}; //@[June 28, 2025] > Really fine xF binning to check which is the best xF binning
  static const Double_t xfbins[NXFBIN+1];// = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.3, 0.5}; //@[August 4, 2025] > New xF binning

  void MergeForTssa( TH1* totalhistinc[][2], TH1* totalhistbg1[][2], TH1* totalhistbg2[][2], TH3* mergedinvmass, TH1* mergedpolblue, TH1* mergedpolyell, TH1* mergedpolblueerr, TH1* mergedpolyellerr );
  static void DoTssaAna( TH1* npi0[][2], TH1* h1_rawasymVphi[][StMuFcsPi0TreeMaker::NXFBIN] );    ///< Compute asymmetry from npi0 array of [blue,yellow][spin up,down] and return an array of raw assymtries by [blue,yellow][energybin]
  static void DoTssaFit(TH1* h1_rawasymVphi[][StMuFcsPi0TreeMaker::NXFBIN], TH1* h1_bluepoldata, TH1* h1_yellowpoldata, TH1* h1_AnResult[]);
  static void DoPi0Fits(TH3* mH1F_invmass, TH1* hist_proj[] );
  static void DoBgCorrectedAn(TH1* h1_invmass_en[], TH1* h1_an_inc, TH1* h1_an_bg, TH1* h1_anresult );
  Int_t ReadPolFile(const char* filename="Run22PolForJobs.txt");        ///< Function to read the polarization data file that is custom made for this class

  static Double_t pol4bg(Double_t* x, Double_t* par);      ///< disjoint 4th order polynomial to exclude signal range
  static Double_t skewgaus(Double_t*x, Double_t* par);     ///< Skewed Gaussian
  static Double_t pol4skewgaus(Double_t*x, Double_t* par); ///< fourth order polynomial + skewed Gaussian
  
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
  TClonesArray* mMixedPhArr = 0;        ///< #TClonesArray of #FcsPhotonCandidate from last event used for event mixing
  //TClonesArray* mBestPharr = 0;         ///< #TClonesArray of #FcsPhotonCandidates which are highest energy pairs
  TClonesArray* mPi0Arr = 0;            ///< Array of #FcsPi0Candidate to store the best pi0 candidates for analysis
  
  TH1* mH1D_Entries = 0;                ///< bin1=Number of events processed no cuts (i.e. "Make" calls), bin2=blue polarization sum no event cut, bin3=yellow polarization sum no event cut, bin4=Number of events trigger cut, bin5=blue polarization sum event cut,  , bin6=yellow polarization sum event cut
  TH1* mH1D_BluePol = 0;                ///< Distribution of Blue beam polarization in %
  TH1* mH1D_YellowPol = 0;              ///< Distribution of Yellow beam Polarization in %
  TH1* mH1D_BluePolErr = 0;             ///< Distribution of Blue beam polarization error in %
  TH1* mH1D_YellowPolErr = 0;           ///< Distribution of Yellow beam Polarization error in %
  TH1* mH1F_Triggers = 0;               ///< Triggers used in analysis
  TH1* mH2F_foundVvertex = 0;           ///< found vertex bit vs. Vertex
  TH1* mH1F_RndmSpin = 0;               ///< Single bin histogram to know if random spins are being used or grabbing from database

  TH1* mH2F_PhotonHeatMap = 0;          ///< Distribution of photons in STAR x,y space
  TH1* mH2F_PhotonHeatMapG = 0;         ///< Distribution of photons in STAR x,y space when energy has a specific value near an energy spike
  TH1* mH2F_PhotonHeatMapB = 0;         ///< Distribution of photons in STAR x,y space when energy has a specific value on an energy spike
  TH1* mH2F_EpdProjHitMap = 0;          ///< Distribution of x,y projections of photon candidates onto STAR EPD plane in x,y space
  TH1* mH2F_EpdProjHitMap_Vcut = 0;     ///< Distribution of x,y projections of photon candidates onto STAR EPD plane in x,y space with cut |vertex|<150cm
  TH1* mH2F_EpdNmip = 0;                ///< Nmip distributions for EPD matched projected clusters (x-axis bin 1) and points (y-axis bin 2)

  TH1* mH1F_ClusterEnergy = 0;          ///< All Cluster energy
  TH1* mH1F_PointEnergy = 0;            ///< All Point energy
  TH1* mH2F_Energy_ph1Vph2 = 0;         ///< Histogram of two photons energy used in reconstruction
  TH1* mH1F_NBadEpdProj = 0;            ///< Number of points that did not have a valid projection to an EPD tile in a given event
  TH1* mH1F_NBadEpdProjVcut = 0;        ///< Number of points that did not have a valid projection to an EPD tile in a given event with cut |vertex|<150

  TH1* mH1F_PointMult = 0;              ///< Raw point multiplicity in event

  TH1* mH2F_PointProj_nmipValldx=0;        ///< nMIP vs. FCS projected point to EPD x-position minus EPD x-position of all hits
  TH1* mH2F_PointProj_nmipValldy=0;        ///< nMIP vs. FCS projected point to EPD y-position minus EPD y-position of all hits
  TH1* mH2F_PointProj_nmipValldr=0;        ///< nMIP vs. FCS projected point to EPD, r difference to all other hits
  TH1* mH2F_PointProj_nmipValldphi=0;      ///< nMIP vs. FCS projected point to EPD, angle difference between all other hits
  TH1* mH2F_MixedPointProj_nmipValldr=0;        ///< Mixed event nMIP vs. FCS projected point to EPD, r difference to all other hits
  TH1* mH2F_MixedPointProj_nmipValldphi=0;      ///< Mixed event nMIP vs. FCS projected point to EPD, angle difference between all other hits
  TH1* mH2F_PointProj_nmipVtiledx=0;       ///< nMIP vs. FCS projected point to EPD x-position minus EPD x-position of center tile
  TH1* mH2F_PointProj_nmipVtiledy=0;       ///< nMIP vs. FCS projected point to EPD y-position minus EPD y-position of center tile
  TH1* mH2F_PointProj_nmipVtiledr=0;       ///< nMIP vs. FCS projected point to EPD, r difference to tile hit
  TH1* mH2F_PointProj_nmipVtiledphi=0;     ///< nMIP vs. FCS projected point to EPD, angle difference to tile hit
  TH1* mH2F_MixedPointProj_nmipVtiledr=0;       ///< Mixed event nMIP vs. FCS projected point to EPD, r difference to tile hit
  TH1* mH2F_MixedPointProj_nmipVtiledphi=0;     ///< Mixed event nMIP vs. FCS projected point to EPD, angle difference to tile hit

  TH1* mH2F_PointProj_LowMult_nmipValldr=0;          ///< nMIP vs. FCS projected point to EPD, r difference to tile hit
  TH1* mH2F_PointProj_LowMult_nmipValldphi=0;        ///< nMIP vs. FCS projected point to EPD, angle difference to tile hit
  TH1* mH2F_MixedPointProj_LowMult_nmipValldr=0;     ///< Mixed event nMIP vs. FCS projected point to EPD, r difference to tile hit
  TH1* mH2F_MixedPointProj_LowMult_nmipValldphi=0;   ///< Mixed event nMIP vs. FCS projected point to EPD, angle difference to tile hit
  
  //TH1* mH2F_Pi0HeatMap = 0;             ///< x,y locations of BestPi0
  TH1* mH1F_AllPi0Mult = 0;             ///< pi0 multiplicity for all pairs of points
  TH1* mH1F_AllPi0Zgg = 0;              ///< Zgg of pi0 for all pairs of points
  TH1* mH2F_AllPi0_etaVphi = 0;         ///< eta vs. phi for all pairs of points
  TH1* mH1F_AllPi0En = 0;               ///< Energy for all pairs of points
  TH1* mH1F_AllPi0Pt = 0;               ///< Pt for all pairs of points
  TH1* mH1F_AllPi0Mass = 0;             ///< Invariant mass for all pairs of points

  TH1* mH1F_NoEpdCutPi0Mult = 0;        ///< Pi0 multiplicity after all cuts except Epd nMIP one
  TH1* mH1F_NoEpdCutZgg = 0;         ///< Pi0 Zgg after all cuts except Epd nMIP one
  TH1* mH2F_NoEpdCut_etaVphi = 0;       ///< Pi0 eta vs. phi after all cuts except Epd nMIP one
  TH1* mH1F_NoEpdCutEn = 0;             ///< Pi0 energy after all cuts except Epd nMIP one
  TH1* mH1F_NoEpdCutPt = 0;             ///< Pi0 pt after all cuts except Epd nMIP one
  TH1* mH1F_NoEpdCutAllMass = 0;        ///< Invariant Mass of all point pairs after all cuts except the EPD nmip cuts

  //TH1* mH1F_EpdPhInvMass = 0;           ///< Invariant mass with EPD nmip cut to isolate uncharged particles (EpdPh stands for Epd Photon, as in uncharged particle)
  //TH1* mH2F_EpdPhHeatMap = 0;           ///< x,y locations for EpdPh
  TH1* mH1F_EpdPhPi0Mult = 0;           ///< Pi0 Multiplicity for EpdPh
  TH1* mH1F_EpdPhZgg = 0;               ///< Zgg of for EpdPh
  TH1* mH2F_EpdPh_etaVphi = 0;          ///< Azimuthal angle for EpdPh
  //TH1* mH1F_EpdPhEta = 0;             ///< Psuedorapidity for EpdPh
  TH1* mH1F_EpdPhEn = 0;                ///< Energy for EpdPh
  TH1* mH1F_EpdPhPt = 0;                ///< Pt for EpdPh
  TH1* mH1F_EpdPhAllMass = 0;           ///< Invariant mass of all point pairs with EPD nmip cut to isolate uncharged particles

  //TH1* mH1F_EpdChInvMass = 0;           ///< Invariant mass with EPD nmip cut to isolate charged particles (EpdCh stands for Epd Charged, as in charged particle)
  //TH1* mH2F_EpdChHeatMap = 0;           ///< x,y locations for EpdCh
  TH1* mH1F_EpdChPi0Mult = 0;           ///< Pi0 multiplicty for EpdCh
  TH1* mH1F_EpdChZgg = 0;               ///< Zgg for EpdCh
  TH1* mH2F_EpdCh_etaVphi = 0;          ///< Psuedorapidity (eta) vs. Azimuthal (phi) angle for EpdCh
  TH1* mH1F_EpdChEn = 0;                ///< Energy for EpdCh
  TH1* mH1F_EpdChPt = 0;                ///< Pt for EpdCh
  TH1* mH1F_EpdChAllMass = 0;           ///< Invariant mass of all point pairs with EPD nmip cut to isolate charged particles

  TH1* mH1F_EpdSinglePhPi0Mult = 0;     ///< Pi0 multiplicty for all point pairs that pass a single photon requirement on Epd
  TH1* mH1F_EpdSinglePhZgg = 0;         ///< Zgg for all point pairs that pass a single photon requirment on Epd
  TH1* mH2F_EpdSinglePh_etaVphi = 0;    ///< eta V phi for all point pairs that pass a single photon requirement on Epd
  TH1* mH1F_EpdSinglePhEn = 0;          ///< Energy for all point pairs that pass a single photon requirement on Epd
  TH1* mH1F_EpdSinglePhPt = 0;          ///< Pt for all point pairs that pass a single photon requirement on Epd
  TH1* mH1F_EpdSinglePhAllMass = 0;     ///< Invariant mass of all point pairs with EPD nmip cut on a single photon that passed the photon level cut

  TH1* mH1F_EpdSingleChPi0Mult = 0;     ///< Pi0 multiplicty for all point pairs that pass a single electron requirement on Epd
  TH1* mH1F_EpdSingleChZgg = 0;         ///< Zgg for all point pairs that pass a single electron requirment on Epd
  TH1* mH2F_EpdSingleCh_etaVphi = 0;    ///< eta V phi for all point pairs that pass a single electron requirement on Epd
  TH1* mH1F_EpdSingleChEn = 0;          ///< Energy for all point pairs that pass a single electron requirement on Epd
  TH1* mH1F_EpdSingleChPt = 0;          ///< Pt for all point pairs that pass a single electron requirement on Epd
  TH1* mH1F_EpdSingleChAllMass = 0;     ///< Invariant mass of all point pairs with EPD nmip cut on a single electron that passed the electron level cut

  //Also separate low point multiplicity events
  //TH1* mH1F_2ndBestPi0Mass = 0;
  //TH1* mH1F_3rdBestPi0Mass = 0;
  //TH1* mH1F_LowPointMult = 0;
  //Get pt thresholds; Plot x,y positions of pi0s using the projection from px,py,pz;make plot of epd bad projections on same pad also plot vs. npointmult?; do run by run qa

  static const short NEPDCUTS = 8;
  TObjArray* mH1F_InvMassEpdCuts[2];         ///< Invariant Mass using different epd nmip cuts and all triggers or only EM triggers
  //TH1* mH1F_InvMassZggCuts[2][8];          ///< Invariant Mass using different zgg cuts and all triggers or only EM triggers
  //TH1* mH1F_InvMassEnBins[2][4][3];        ///< Invariant Mass in different energy bins after each cut, fidicual volume cut, Zgg cut, Epd photon cut, and all triggers or only EM triggers
  
  TObjArray* mH1F_InvMassAllCuts = 0;              ///< Invariant Mass of all potential pi0s after all cuts applied and the "neutral" particle criteria for EPD nmip, one for every EM trigger
  TObjArray* mH1F_Pi0MultAllCuts = 0;              ///< Number of "good pi0s" i.e. number of potential pi0s after all cuts applied
  //(TH1* mH1F_AllPi0_xF = 0;             ///< Feynman-x (xF) of all pi0s
  //TH1* mH1F_NFoundPhiBin = 0;           ///< Number of valid phi bins found. This is a cross check to make sure that I am not double counting pi0s and finding more than one valid bin when I loop over the phi bins
  TObjArray* mH1F_AllCuts_xF = 0;               ///< Feynman-x (xF) of the pi0s that pass all the cuts
  TObjArray* mH1F_AllCuts_xFZoom = 0;           ///< Feynman-x (xF) of the pi0s that pass all the cuts on a smaller scale and more bining
  TObjArray* mH1F_AllCuts_Zgg = 0;              ///< Zgg after all cuts
  TObjArray* mH1F_AllCuts_Dgg = 0;              ///< Dgg after all cuts
  TObjArray* mH1F_AllCuts_Pi0En = 0;            ///< Pi0 energy distribution with all cuts
  TObjArray* mH2F_AllCuts_Pi0_massVen = 0;      ///< Pi0 Invariant mass vs. energy with all cuts
  TObjArray* mH2F_AllCuts_Pi0_xfVen = 0;        ///< Pi0 xF vs. Energy
  TObjArray* mH2F_AllCuts_Pi0_ptVeta = 0;       ///< Pi0 p_T vs. psuedorapidity (eta)
  //TH1* mH1F_AllCuts_Pi0Phi = 0;         ///< Pi0 phi distribution with all cuts
  TObjArray* mH2F_AllCuts_Pi0_etaVphi = 0;      ///< Pi0 eta vs. phi distributions, phi binning matches #NPHIBIN after all cuts
  //TH1* mH2F_AllCuts_Poi_yVx = 0;        ///< Point y vs. x distributions after all cuts this gets complicated because you have two point
  TObjArray* mH2F_AllCuts_Pi0_yVx = 0;          ///< Pi0 FCS projected y vs. x distributions after all cuts
  //TH1* mH1F_InvMassAllCutsByEnByPhi[NENERGYBIN][NPHIBIN]; ///< Invariant mass of reconstructed pions using all cuts by energy and phi bin
  //TH1* mH1F_NPi0ByEnByPhi[NENERGYBIN][NPHIBIN];  ///< Number of pions in a given energy and phi bin. The hisotram represents the split by spin state (0-10, 20-40, 40-60, 80-100, 100+)
  TH1* mH3F_AllCutsInvMass_xfVphi = 0;      ///< Pi0 invariant mass after all cuts vs. xf vs. phi used for the fits to estimate background
  TH1* mH2F_NPi0Inc_xfVphi[2][2];           ///< x_F and phi where pi0 candidate was found [blue,yellow][up,down] for invariant mass range of 0.1-0.2
  TH1* mH2F_NPi0Bg1_xfVphi[2][2];           ///< x_F and phi where pi0 candidate was found [blue,yellow][up,down] for invariant mass range of 0.3-0.4
  TH1* mH2F_NPi0Bg2_xfVphi[2][2];           ///< x_F and phi where pi0 candidate was found [blue,yellow][up,down] for invariant mass range of 0.7-0.9
  //Add mass vs. energy hisotgram??

  TGraphErrors* mGE_AllCuts_InvMass = 0;
  TGraphErrors* mGE_AllCuts_Pi0En = 0;
  //TGraph* mG_TSSA = 0;

  bool mValidTrigFound = false;
  bool mEmTrigFound = false;
  // Variables to keep track if a given trigger was found in an event and can be used to fill the trigger sepcific histograms. The logic is that less than 0 means trigger was not fired. greater than or equal to 0 means trigger was fired. If trigger is fired then the value gets set to the value that corresponds to the index of a histogram array to be used later for filling. The 0 index is for all triggers
  short mTrigEm0 = -1;        ///< Gets set to 1 if EM0 or EM0_tpc trigger was fired
  short mTrigEm1 = -1;        ///< Gets set to 2 if EM1 or EM1_tpc trigger was fired
  short mTrigEm2 = -1;        ///< Gets set to 3 if EM2 or EM2_tpc trigger was fired
  short mTrigEm3 = -1;        ///< Gets set to 4 if EM3 or EM3_tpc trigger was fired

  short mFoundVertex = 0;              ///< Bit vector for knowing which vertex was used 0 means no vertex, 1=Vpd,2=Epd,4=Bbc
  Double_t mUseVertex = -999.0;        ///< Vertex to use in analysis
  Double_t mVertexCutLow = -150.0;     ///< Variables for z vertex cuts at low end in cm
  Double_t mVertexCutHigh = 150.0;     ///< Variables for z vertex cuts at high end in cm
  
  Double_t mEnCut = 1;                  ///< Energy Cut for #FcsPhotonCandidates
  Double_t mEpdNmipCut = 0.4;           ///< Cut on EPD nmip to classify cluster or point as charged or uncharged
  UShort_t mTreeOnBitMap = 0x7;         ///< Turn on or off branches in the pi0 tree. first bit is events, second bit is photon branch, third bit is pi0 branch. Turn on all branches by default
  

  /*Simple data struct to hold polarization information from file only important values are kept*/
  struct PolData
  {
    Int_t mFillNum = 0;              //! Fill Number
    Int_t mBeamEn = 0;               //! Beam energy
    Int_t mStartTime = 0;            //! Start time of fill
    Double_t mBlueP0 = 0;            //! Blue beam intial polarization in %
    Double_t mBlueErrP0 = 0;         //! Blue beam intial polarization error in %
    Double_t mBluedPdT = 0;          //! Blue beam polarization decay in %/hour
    Double_t mBlueErrdPdT = 0;       //! Blue beam polarization decay error in %/hour
    Double_t mYellowP0 = 0;          //! Yellow beam intial polarization in %
    Double_t mYellowErrP0 = 0;       //! Yellow beam intial polarization error in %
    Double_t mYellowdPdT = 0;        //! Yellow beam polarization decay in %/hour
    Double_t mYellowErrdPdT = 0;     //! Yellow beam polarization decay error in %/hour

    void Print() const;
  };
  std::map<Int_t,PolData*> mPolarizationData;   ///< Map of polarization data from file with fill number as key to quickly look up from data structure

  //void GetAdjacentEpdTile(int pp, int tt, int& pp_adj, int& tt_adj) const;
  TPolyLine* EpdTilePoly(short pp, short tt) const;       ///< Returns a new polyline using corners from StEpdGeom
  TPolyLine* EpdCCWOuterCorner(short pp, short tt) const;  ///< Return a new polyline  using adjacency of outer CCW (see #StMuEpdRun22QaMaker)
  TPolyLine* EpdCWOuterCorner(short pp, short tt) const;
  TPolyLine* EpdCWInnerCorner(short pp, short tt) const;
  TPolyLine* EpdCCWInnerCorner(short pp, short tt) const;
  std::map<Int_t,TPolyLine*> mEpdTileMap;    ///< EPD "tile key" to polyline for drawing
  
private:
  HistManager* mHists = 0;            ///< Manage loading and saving histograms
  bool mInternalHists = false;        ///< Boolean to keep track if mHists was added externally or an internal one was created
  TRandom3 mSpinRndm;                 ///< Spin state randomizer

  Double_t mOldVertex = -999.0;                ///< Vertex from last event needed for event mixing
  Int_t mNOldPoints = 0;                       ///< Number of points in previous event

  static Int_t GetColor(Double_t Value, Double_t MinVal, Double_t MaxVal);

  ClassDef(StMuFcsPi0TreeMaker, 4)
};

#endif
