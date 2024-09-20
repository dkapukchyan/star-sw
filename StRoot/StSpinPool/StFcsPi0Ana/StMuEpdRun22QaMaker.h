/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To generate a ROOT file of histograms related to the FCS and event data in the produced MuDsts from RHIC Run 22 that can be used to do quality assurance (QA) on the data.

  DESCRIPTION
  This class inherits from StMaker and contains many histograms to be used for Quality Assurance (QA) of MuDst files that contain Forward Calorimeter System (FCS) data. It uses #LoadHists() together with the functions #AddH1F(), #AddH2F(), #AddH1FArr(), #AddH2FArr() to ease histogram creation and management. These functions should be used exclusively in this and inherited classes to automate histogram management. #FillEventInfo() is used (and can be re-implemented) to fill QA histograms related to event information. #FillFcsInfo() is used (and can be re-implemented) to fill QA histograms related to FCS information. The idea is that inherited classes don't need to implement the #Init() function only the #LoadHists() function. Also can be used to do some EPD Qa with the #FillEpdInfo() function.

  LOG
  @[August 14, 2024] > Copied from StMuFcsRun22QaMaker.h into its own class
  @[August 29, 2024] > Added some histograms to check the correlation between various vertex methods. Added a variable #mEpdVertex to store the found vertex so it can be retrieved outside the maker if needed. Added a draw function for the vertex histograms. Added multiplicity histograms with cuts on nMIP. Some cleanup.
  @[August 30, 2024] > Changed default vertex for all detectors to -999. Modified some drawing options. Small fixes
  @[September 9, 2024] > Clean up extraneous code and added #getFileName() for mHists
  @[September 18, 2024] > Made the variables related to computing the vertex members of the class and added get functions for them
  
  Do DEP calib of EPD chs, bunch xing analysis for spin. Change some plots so they use logz and move/remove the stats box for some of hte 2d histograms when plotting. Show on the fly EPD MIP peak locations and valleys
 */

#ifndef STMUEPDRUN22QAMAKER_HH
#define STMUEPDRUN22QAMAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TRandom3.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

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
#include "StEpdDbMaker/StEpdDbMaker.h"
#include "StEpdHitMaker/StEpdHitMaker.h"

//Custom headers in this folder
#include "HistManager.h"


class StMuEpdRun22QaMaker : public StMaker
{
 public:
  StMuEpdRun22QaMaker(const char* name = "FcsEpdRun22Qa");
  ~StMuEpdRun22QaMaker();

  virtual Int_t Init();
  virtual Int_t InitRun(int runnumber);
  virtual Int_t Make();
  virtual Int_t Finish();

  const char* getFileName(){ if( mHists!=0 ){ return mHists->getFileName(); }else{return 0;} }
  void setHistManager( HistManager* hm );
  virtual UInt_t LoadHists(TFile* file);
  
  void setEpdTacAdcOn(bool value=true) { mEpdTacAdcOn = value; }

  Int_t epdTacEarlyE(){ return mCutEarliestTacE; }
  Int_t epdTacEarlyW(){ return mCutEarliestTacW; }
  Double_t epdTacAvgE()  { return mCutAvgTacE; }
  Double_t epdTacAvgW()  { return mCutAvgTacW; }
  Double_t epdVertex()   { return mEpdVertex; }
  
  void DrawEpdAllQa(TCanvas* canv, const char* savename);      //!< Call all the draw functions in same order as below, should be some kind of pdf as it will encompass many pages
  void DrawVertex(TCanvas* canv, const char* savename);        //!< Draw 2D vertex correlation histograms, single page
  void DrawEpdHitQa(TCanvas* canv, const char* savename);      //!< Draw Multiplicty and nMIP histograms, single page
  void DrawEpdTacQa(TCanvas* canv, const char* savename);      //!< Draw the histograms related to the TAC difference, single page
  void DrawEpdTacCutQa(TCanvas* canv, const char* savename);   //!< Draw the histograms related to the TAC difference with cuts, single page
  void DrawEpdTacAdcQa(TCanvas* canv, const char* savename);   //!< Draw the TAC vs. ADC histograms for all channels, savename should be some kind of pdf as it will encompass many pages

protected:
  StMuDstMaker* mMuDstMkr = 0;
  StMuDst* mMuDst = 0;
  StMuEvent* mMuEvent = 0;
  const StTriggerData* mTrigData = 0;
  StRunInfo* mRunInfo = 0;
  
  //StEpdGeom* mEpdGeo=0;
  TClonesArray* mMuEpdHits = 0;
  StEpdHitMaker* mEpdHitMkr = 0;
  StEpdCollection* mEpdColl = 0;
  
  //Data to save
  //TString mFileName;
  //TFile* mFile_Output = 0; //!< TFile to save all the data

  virtual Int_t FillEpdInfo();

  /*
  TH1* mH1F_Entries = 0;              //!< Number of events processed no cuts (i.e. "Make" calls)
  TH1* mH1F_Triggers = 0;             //!< Triggers in the events
  TH1* mH1F_VertexPrimZ = 0;          //!< Vertex histograms from Primary Vertex
  TH1* mH1F_VertexVpd = 0;            //!< Vertex histograms from VPD
  TH1* mH1F_VertexBbc = 0;            //!< Vertex histograms from BBC
  TH1* mH1F_BbcTimeDiff = 0;          //!< BBC Time difference used to compute the vertex
  TH1* mH1F_VertexZdc = 0;            //!< Vertex from ZDC
  */
  TH1* mH1F_VertexEpd = 0;            //!< Vertex from EPD

  TH1* mH2F_VertexZ_vpdVepd = 0;      //!< Correlation histogram between VPD z vertex vs. computed EPD z vertex
  TH1* mH2F_VertexZ_zdcVepd = 0;      //!< Correlation histogram between ZDC z vertex vs. computed EPD z vertex
  TH1* mH2F_VertexZ_bbcVepd = 0;      //!< Correlation histogram between BBC z vertex vs. computed EPD z vertex
  TH1* mH2F_VertexZ_vpdVbbc = 0;      //!< Correlation histogram between VPD z vertex vs. the BBC z vertex
  TH1* mH2F_VertexZ_vpdVzdc = 0;      //!< Correlation histogram between VPD z vertex vs. the ZDC z vertex
  TH1* mH2F_VertexZ_zdcVbbc = 0;      //!< Correlation histogram between ZDC z vertex vs. the BBC z vertex

  TH1* mH1F_Epd_NHits = 0;              //!< Number of hits from EPD collection
  TH1* mH1F_Epd_NHits_Cut = 0;          //!< Number of hits in EPD with nMIP>0.7
  TH1* mH1F_Epd_NHitsWest = 0;          //!< Number of hits from EPD collection only west side
  TH1* mH1F_Epd_NHitsWest_Cut = 0;      //!< Number of hits from EPD collection only west side with nMIP>0.7

  TH1* mH2F_HitEpd_nmipVchkey[2];        //!< Special for checking nmip of a given channel in the EPD. The "chkey" is (supersector-1)*31+(tileid-1)
  TObjArray* mH2F_HitEpd_tacVadcmip[2]; //!< Special for checking EPD TAC vs. ADC/ADC_1mip histograms which may help with slew corrections in the EPD

  TH1* mH2F_Epd_earlywVearlye = 0;      //!< EPD hit with Earliest West TAC vs. Earliest East TAC (Since using common stop early means largest TAC value)
  TH1* mH2F_Epd_avgwVavge = 0;          //!< EPD Averaged West TAC vs. Averaged East TAC
  TH1* mH2F_EpdTacDiff_avgVearly = 0;       //EPD tac difference between West and East tiles (in that order) computed from average tac vs. differences computed using earliest TAC
  TH1* mH2F_EpdCut_earlywVearlye = 0;   //!< EPD hit with Earliest West TAC vs. Earliest East TAC with channel 1<nMIP<15 (Since using common stop early means largest TAC value)
  TH1* mH2F_EpdCut_avgwVavge = 0;       //!< EPD Averaged West TAC vs. Averaged East TAC with channel 1<nMIP<15
  TH1* mH2F_EpdCutTacDiff_avgVearly = 0;  //!< EPD tac difference between West and East tiles (in that order)  computed from average tac vs. differences computed using earliest TAC; computed with cut 1<adcnmip<15 && TAC>50

  bool mEpdTacAdcOn = true;            //!< For turning on/off TAC vs. ADC histograms for all EPD channels

private:
  HistManager* mHists = 0;            //!< Manage loading and saving histograms
  bool mInternalHists = false;        //!< Boolean to keep track if mHists was added externally or an internal one was created
  //TFile* mFileOutput = 0;           //!< For saving histograms not loading
  Double_t mEpdVertex = -999;         //!< Saved vertex from event
  Int_t mCutEarliestTacE  = 0;        //!< Stores the largest found TAC value in EPD East with cuts
  Int_t mCutEarliestTacW  = 0;        //!< Stores the largest found TAC value in EPD West with cuts
  Double_t mCutAvgTacE = 0;           //!< Stores the average TAC value in EPD East with cuts
  Double_t mCutAvgTacW = 0;           //!< Stores the average TAC value in EPD West with cuts

  double mEpdScale = 15.6;             //!< picoSecond/TAC for EPD

  ClassDef(StMuEpdRun22QaMaker,1)

};

#endif
