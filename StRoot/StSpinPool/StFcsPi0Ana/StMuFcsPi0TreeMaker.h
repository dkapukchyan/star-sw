/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpose of this class is to generate a #TTree of photon and pi0 candidates to use in a transverse single spin assymmetry (TSSA) analysis. It will also store basic event level information as well

  DESCRIPTION
  This Maker will loop over all the clusters and points in the MuDst tree and fill a new #TTree with clusters and points as photon candidates (#FcsPhotonCandidate). It will loop over those candidates and pair them together to generate a list of pi0 candidates (#FcsPi0Candidate). It only combines clusters with clusters and points with points. It will only apply an energy cut and trigger selection to filter out what events get saved. Trigger Filter and energy cut can be turned off with appropriate function call. Also stores event level information (#FcsEventInfo) and trigger information (#mTriggers).

  LOG
  @[December 30, 2025] > Copied from "new" StMuFcsTreeMaker and modified to only keep the stuff related to actual pi0 analysis
  @[February 11, 2026] > Obsolete class as part of the experimentation with using #StMuFcsTreeMaker for doing just the pi0 TSSA analysis. Want to keep this for keeping track of the history of how the analysis worked and as a learning tool for how I was experimenting to find a new framework

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
#include "StMuFcsTreeMaker.h"

class StEpdGeom;

class StMuFcsPi0TreeMaker : public StMuFcsTreeMaker {
public:
  
  StMuFcsPi0TreeMaker(const Char_t* name = "MuFcsPi0Maker");
  ~StMuFcsPi0TreeMaker();
  //virtual Int_t Init();
  //virtual Int_t InitRun(int runnumber);
  //virtual void Clear(Option_t* option="");           ///< Gets called before #Make() in StChain::EventLoop()
  virtual Int_t Make();
  
  //virtual Int_t Make_CheckAndSetEpdHit(){return StMuFcsTreeMaker::Make_CheckAndSetEpdHit();}  ///< Loop over #mPhArr and check which photon candidates pair with which EPD tiles
  //virtual Int_t Make_PointPairs(TClonesArray* pointpairs);     ///< Loop over all the photon (point) candidates and make photon pairs with and save them to the TClonesArray that will hold the particle candidates
  virtual Int_t Make_TssaAna(TClonesArray* pointpairs);        ///< Loop over the point pairs and apply some cut criteria to then make histograms for the Transverse Single Spin Asymmetry (TSSA) analysis
  
  //virtual Int_t Finish();

  virtual UInt_t LoadHists(TFile* file);

  void PaintAllPi0(TCanvas* canv,  const char* savename = "testallpi0.png")  const;
  void PaintNoEpdCut(TCanvas* canv,  const char* savename="testnoepdcutpi0.png")  const;
  void PaintEpdPhPi0(TCanvas* canv, const char* savename = "testepdphpi0.png") const;
  void PaintEpdChPi0(TCanvas* canv, const char* savename = "testepdchpi0.png") const;
  void PaintEpdSinglePh(TCanvas* canv, const char* savename = "testepdsingleph.png") const;
  void PaintEpdSingleCh(TCanvas* canv, const char* savename = "testepdsinglech.png") const;
  void PaintPi0Overlap(TCanvas* canv, const char* savename = "testpi0overlap.png") const;
  void PaintInvMassEpdQa(TCanvas* canv, const char* savename = "testinvmasscutqa.png") const;
  void PaintEpdQa(TCanvas* canv, const char* savename = "testepdsingleqa.png") const;
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
  
  void PaintPi0QaForDefense(TCanvas* canv, const char* savename) const;

  //virtual void MergeForTssa( TH1* totalhistinc[][2], TH1* totalhistbg1[][2], TH1* totalhistbg2[][2], TH3* mergedinvmass, TH1* mergedpolblue, TH1* mergedpolyell, TH1* mergedpolblueerr, TH1* mergedpolyellerr );
  
protected:
  
  //TH1* mH1F_PointProjPh_dx=0;        ///< FCS projected point to EPD x-position minus EPD x-position of all hits
  //TH1* mH1F_PointProjPh_dy=0;        ///< projected point to EPD y-position minus EPD y-position of all hits
  //TH1* mH1F_PointProjPh_dr=0;        ///< nMIP vs. FCS projected point to EPD, r difference to all other hits
  //TH1* mH1F_PointProjPh_dphi=0;      ///< nMIP vs. FCS projected point to EPD, angle difference between all other hits

  //TH1* mH1F_PointProjCh_dx=0;        ///< FCS projected point to EPD x-position minus EPD x-position of all hits
  //TH1* mH1F_PointProjCh_dy=0;        ///< projected point to EPD y-position minus EPD y-position of all hits
  //TH1* mH1F_PointProjCh_dr=0;        ///< nMIP vs. FCS projected point to EPD, r difference to all other hits
  //TH1* mH1F_PointProjCh_dphi=0;      ///< nMIP vs. FCS projected point to EPD, angle difference between all other hits
  
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

  //TGraphErrors* mGE_AllCuts_Pi0En = 0;
  //TGraph* mG_TSSA = 0;

  //virtual void CheckInsideEpdTile(FcsPhotonCandidate* photon, Double_t projx, Double_t projy) const;  ///< This algorithm will first check if a point lies inside any of the EPD tiles. If yes, then the first entry in the #FcsCandidate::mEpdHitNmip array gets set to the nMIP value of that tile. If the point does not fall into any of the EPD tiles then it checks all the CCW regions for which the point can be in. These will get saved to the second and beyond entries. It will then check all the second to beyond entries for the tile that is closest to the point and set that nMIP value to the first entry. This ensures that I only need to check the first entry of the nMIP array when distinguishing photon candidates.

  //Since "old" StMuFcsPi0TreeMaker is "new" StMuFcsTreeMaker, increment this to 6 to prevent clashing with older versions of "StMuFcsPi0TreeMaker
  ClassDef(StMuFcsPi0TreeMaker, 6)
};

#endif
