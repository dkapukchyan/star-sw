/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpose of this class is to do the the mixed event testing to compare hit EPD tiles to FCS points from the previous event

  DESCRIPTION
  This Maker inherits from #StMuFcsPi0TreeMaker and adds histogram,s and code to do the mixed event analysis. In particular the make function will only do mixed events and exit

  LOG
  @[December 23, 2025] > Copied from #StMuFcsPi0TreeMaker and modified to use new class name
  @[January 19, 2026] > Obsolete class as part of the experimentation with using #StMuFcsTreeMaker for doing the EPD mixed event analysis. Still haven't implemented into the new framework and want to keep this for keeping track of the history of how the analysis worked and as a learning tool for how I was experimenting to find a new framework

*/


#ifndef STMUFCSPOINTEPDMIXEDEVENTMAKER_HH
#define STMUFCSPOINTEPDMIXEDEVENTMAKER_HH

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
#include "StSpinPool/StFcsPi0Ana/StMuFcsPi0TreeMaker.h"

class StEpdGeom;

class StMuFcsPointEpdMixedEventMaker : public StMuFcsPi0TreeMaker {
public:
  
  StMuFcsPointEpdMixedEventMaker(const Char_t* name = "MuFcsPointEpdMixedEvent");
  ~StMuFcsPointEpdMixedEventMaker();
  virtual Int_t Init();
  virtual Int_t Make();

  //virtual Int_t Make_VertexInfo();         ///< Get and set members related to the collision vertex
  //virtual Int_t Make_FillFcsClusPoint();   ///< Get and set the clusters and points to #mPhArr
  virtual Int_t Make_CheckAndSetEpdHit();  ///< Loop over #mPhArr and check which photon candidates pair with which EPD tiles
  //virtual Int_t Make_PointPairs(TClonesArray* pointpairs);     ///< Loop over all the photon (point) candidates and make photon pairs with and save them to the TClonesArray that will hold the particle candidates
  //virtual Int_t Make_TssaAna(TClonesArray* pointpairs);        ///< Loop over the point pairs and apply some cut criteria to then make histograms for the Transverse Single Spin Asymmetry (TSSA) analysis
  
  //virtual Int_t Finish();

  //UInt_t LoadDataFromFile(TFile* file);//, TTree&* tree, FcsEventInfo&* evt,Int_t& ntrig, Int_t&* triggers,  TClonesArray&* pharr, TClonesArray&* pi0arr, TH1&* hist=0):
  virtual UInt_t LoadHists(TFile* file);

  //virtual void Print(Option_t* opt="") const; //"e" for event, "t" for trigger, "g" for photon, "p" for pi0, "a" for all

  void PaintEpdAllDistQa(TCanvas* canv, const char* savename = "testepdalldistqa.png") const;
  void PaintEpdAllDistQaLowMult(TCanvas* canv, const char* savename = "testepdalldistqalowmult.png" ) const;
  void PaintEpdTileDistQa(TCanvas* canv, const char* savename = "testepdtiledistqa.png") const;
  void PaintEpdDistAnaQa(TCanvas* canv, const char* savename = "testepddistanaqa.png") const;
  
protected:
  TClonesArray* mMixedPhArr = 0;        ///< #TClonesArray of #FcsPhotonCandidate from last event used for event mixing

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

private:
  Double_t mOldVertex = -999.0;                ///< Vertex from last event needed for event mixing
  Int_t mNOldPoints = 0;                       ///< Number of points in previous event
  
  ClassDef(StMuFcsPointEpdMixedEventMaker, 1)
};

#endif

