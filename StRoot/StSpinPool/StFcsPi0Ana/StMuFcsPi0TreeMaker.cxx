#include "StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StEventTypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StSpinPool/StFcsQaMaker/StFcsQaMaker.h"
#include "StSpinPool/StFcsRawDaqReader/StFcsRawDaqReader.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StThreeVectorF.hh"
#include "Stypes.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TMarker.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"

#include "StMuFcsPi0TreeMaker.h"

ClassImp(StMuFcsPi0TreeMaker)

//const double StMuFcsPi0TreeMaker::xfbins[] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52};
const Double_t StMuFcsPi0TreeMaker::xfbins[] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.3, 0.5}; //@[August 4, 2025] > New xF binning

StMuFcsPi0TreeMaker::StMuFcsPi0TreeMaker(const Char_t* name) : StMaker(name)
{
  mSpinRndm.SetSeed(0);
  memset(mTriggers,0,sizeof(mTriggers));
  memset(mH1F_InvMassEpdCuts,0,sizeof(mH1F_InvMassEpdCuts));
  //memset(mH1F_InvMassZggCuts,0,sizeof(mH1F_InvMassZggCuts));
  //memset(mH1F_InvMassPtCuts, 0,sizeof(mH1F_InvMassPtCuts));
  //memset(mH1F_InvMassAllCutsByEnByPhi, 0,sizeof(mH1F_mH1F_InvMassAllCutsByEnByPhi));
  //memset(mH1F_NPi0ByEnByPhi, 0,sizeof(mH1F_NPi0ByEnByPhi));
  memset(mH2F_NPi0Inc_xfVphi, 0,sizeof(mH2F_NPi0Inc_xfVphi));
  memset(mH2F_NPi0Bg1_xfVphi, 0,sizeof(mH2F_NPi0Bg1_xfVphi));
  memset(mH2F_NPi0Bg2_xfVphi, 0,sizeof(mH2F_NPi0Bg2_xfVphi));
}

StMuFcsPi0TreeMaker::~StMuFcsPi0TreeMaker()
{
  delete mEvtInfo;
  delete mPhArr;
  delete mPi0Arr;
  delete mPi0Tree;
  if( mInternalHists ){ delete mHists; } //This deletes the file too
  for( auto itr=mPolarizationData.begin(); itr!=mPolarizationData.end(); ++itr){ delete itr->second; }
  mPolarizationData.clear();
  //if( mFile_Output!=0 ){ mFile_Output->Close(); }
  //delete mFile_Output;
}

/*#ifndef __CINT__
void StMuFcsPi0TreeMaker::SetTrigs(const char* trigname,...)
{
  va_list args;
  va_start(args,trigname);
  
  char* name = va_arg(args,char*);
  mTargetTrig.emplace_back(name);

  va_end(args);
}
#endif*/

void StMuFcsPi0TreeMaker::setEventBit(bool val)
{
  if( val ){ mTreeOnBitMap |= 0x01; }
  else{ mTreeOnBitMap &= ~(0x01); }
}

void StMuFcsPi0TreeMaker::setPhotonOn(bool val)
{
  if( val ){ mTreeOnBitMap |= 0x02; }
  else{ mTreeOnBitMap &= ~(0x02); }
}

void StMuFcsPi0TreeMaker::setPi0On(bool val)
{
  if( val ){ mTreeOnBitMap |= 0x04; }
  else{ mTreeOnBitMap &= ~(0x04); }
}

/*bool StMuFcsPi0TreeMaker::checkTreeOnBit(UShort_t bit)
{
  if( mTreeOnBitMap & (1<<bit) ){ return true; }
  else{ return false; }
  }*/

bool StMuFcsPi0TreeMaker::isEventOn() const
{
  if( mTreeOnBitMap & 0x01 ){ return true; }
  else{ return false; }
}

bool StMuFcsPi0TreeMaker::isPhotonOn() const
{
  if( mTreeOnBitMap & 0x02 ){ return true; }
  else{ return false; }
}

bool StMuFcsPi0TreeMaker::isPi0On() const
{
  if( mTreeOnBitMap & 0x04 ){ return true; }
  else{ return false; }
}

void StMuFcsPi0TreeMaker::setHistManager( HistManager* hm )
{
  if( mInternalHists ){ delete mHists; mHists = 0; }
  mInternalHists = false;
  if( mHists!=0 ){ LOG_WARN << "StMuFcsPi0TreeMaker::setHistManager() - HistManager exists and is external - no changes made" << endm; return; }
  else{ mHists = hm; }
}

UInt_t StMuFcsPi0TreeMaker::LoadHists( TFile* file )
{
  if( mHists==0 ){ return 0; }
  UInt_t loaded = 0;
  //loaded += mHists->AddH1F(file,mH1F_Entries,"H1_Entries", "Number of entries;", 3,-0.5, 2.5);
  loaded += mHists->AddH1D(file,mH1D_Entries,     "H1D_Entries", "Number of entries", 1,-0.5,0.5 );
  loaded += mHists->AddH1D(file,mH1D_BluePol,     "H1D_BluePol",     "Blue Beam Polarization;%", 1000,0,100 );
  loaded += mHists->AddH1D(file,mH1D_YellowPol,   "H1D_YellowPol",   "Yellow Beam Polarization;%", 1000,0,100 );
  loaded += mHists->AddH1D(file,mH1D_BluePolErr,  "H1D_BluePolErr",  "Blue Beam Polarization Error;%", 500,0,10 );
  loaded += mHists->AddH1D(file,mH1D_YellowPolErr,"H1D_YellowPolErr","Yellow Beam Polarization Error;%", 500,0,10 );
  //Below trigger histogram copied from StMuFcsRun22QaMaker
  loaded += mHists->AddH1F(file,mH1F_Triggers,"H1F_Triggers","Triggers;;",65,0,65);
  if( mFcsTrigMap!=0 ){
    int trigsize = mFcsTrigMap->sizeOfTriggers();
    //std::cout << "|trigsize:"<<trigsize << std::endl;
    if( trigsize==64 ){ //Check to make sure all triggers are in the map and it matches the bin size
      for( int i=0; i<trigsize; ++i ){
	//std::cout << "|itrig:"<<i << "|trigname:"<< mFcsTrigMap->triggerName(i) << std::endl;
	mH1F_Triggers->GetXaxis()->SetBinLabel(i+1,mFcsTrigMap->triggerName(i)); //Bin numbers are offset by 1
      }
    }
    mH1F_Triggers->GetXaxis()->SetBinLabel(65,"NF");  //Last bin is named "NF" for not found which is what nameFromId() returns if trigger is not in the map. This way if no map was loaded but an StFcsRun22TriggerMap was found, searching for triggers will return "NF"
  }
  loaded += mHists->AddH2F(file,mH2F_foundVvertex,"H2F_foundVvertex", "Used vertex bit vs. Vertex that was used in Pi0 Reconstruction;Vertex (cm);found (bit) 0=NA,1=VPD,2=EPD,4=BBC", 50,-200,200, 5,0,5);
  loaded += mHists->AddH1F(file,mH1F_RndmSpin,"H1F_RndmSpin","Histogram to know if using random spin or not",2,0,2);

  loaded += mHists->AddH2F(file,mH2F_PhotonHeatMap,"H2F_PhotonHeatMap","Distribution of photons in STAR x,y space;x (cm);y (cm)", 400,-200,200, 300,-150,150);
  loaded += mHists->AddH2F(file,mH2F_PhotonHeatMapG,"H2F_PhotonHeatMapG","Distribution of photons in STAR x,y space (No Spike);x (cm);y (cm)", 400,-200,200, 300,-150,150);
  loaded += mHists->AddH2F(file,mH2F_PhotonHeatMapB,"H2F_PhotonHeatMapB","Distribution of photons in STAR x,y space (Spike);x (cm);y (cm)", 400,-200,200, 300,-150,150);
  loaded += mHists->AddH2F(file,mH2F_EpdProjHitMap,"H2F_EpdProjHitMap","Distribution of x,y projections of photon candidates onto STAR EPD plane;x (cm);y (cm)", 300,-150,150, 200,-100,100);
  loaded += mHists->AddH2F(file,mH2F_EpdProjHitMap_Vcut,"H2F_EpdProjHitMap_Vcut","Distribution of x,y projections of photon candidates onto STAR EPD plane with |vertex|<150cm;x (cm);y (cm)", 300,-150,150, 200,-100,100);
  loaded += mHists->AddH2F(file,mH2F_EpdNmip,"H2F_EpdNmip","Distribution of nmip values from matching projected clusters and points;;Nmip",2,0,2, 70,0,7);
  mH2F_EpdNmip->GetXaxis()->SetBinLabel(1,"Clusters");
  mH2F_EpdNmip->GetXaxis()->SetBinLabel(2,"Points");

  loaded += mHists->AddH1F(file,mH1F_ClusterEnergy,"H1F_ClusterEnergy","Energy of FCS Clusters;Energy (GeV);", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_PointEnergy,"H1F_PointEnergy","Energy of FCS Points;Energy (GeV);", 1000,0,200);
  loaded += mHists->AddH2F(file,mH2F_Energy_ph1Vph2,"H2F_Energy_ph1Vph2","Energy of photon 1 vs. photon 2 no EPD cuts;Energy Highest Photon (GeV);Energy Next Highest (GeV)", 1000,0,200, 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_NBadEpdProj,"H1F_NBadEpdProj","Number of points in an event that did not project back to a valid EPD tile;;",15,0,15);
  loaded += mHists->AddH1F(file,mH1F_NBadEpdProjVcut,"H1F_NBadEpdProjVcut","Number of points in an event that did not project back to a valid EPD tile with |vertex|<150;;",15,0,15);

  loaded += mHists->AddH1F(file,mH1F_PointMult,"H1F_PointMult","Point Multiplicity with only an energy cut;Point Multiplicity", 30,0,30);

  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldx,"H2F_PointProj_nmipValldx","nMIP vs. FCS projected point to EPD x minus EPD x of all hits;dX (cm);nmip", 200,-100,100, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldy,"H2F_PointProj_nmipValldy","nMIP vs. FCS projected point to EPD y minus EPD y of all hits;dY (cm);nmip", 200,-100,100, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldr,"H2F_PointProj_nmipValldr","nMIP vs. FCS projected point to EPD, polar distance to all other hits;dR (cm);nmip", 200,-100,100, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldphi,"H2F_PointProj_nmipValldphi","nMIP vs. FCS projected point to EPD, angle difference to all other hits;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipValldr,"H2F_MixedPointProj_nmipValldr","Mixed event nMIP vs. FCS projected point to EPD, polar distance to all other hits;dR (cm);nmip", 200,-100,100, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipValldphi,"H2F_MixedPointProj_nmipValldphi","Mixed event nMIP vs. FCS projected point to EPD, angle difference to all other hits;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);

  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledx,"H2F_PointProj_nmipVtiledx","nMIP vs. FCS projected point to EPD x minus EPD x of proj tile;dX (cm);nmip", 80,-20,20, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledy,"H2F_PointProj_nmipVtiledy","nMIP vs. FCS projected point to EPD y minus EPD y of proj tile;dY (cm);nmip", 80,-20,20, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledr,"H2F_PointProj_nmipVtiledr","nMIP vs. FCS projected point to EPD, polar distance to projected hit;dR (cm);nmip", 80,-20,20, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledphi,"H2F_PointProj_nmipVtiledphi","nMIP vs. FCS projected point to EPD, angle difference to proj tile;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipVtiledr,"H2F_MixedPointProj_nmipVtiledr","Mixed event nMIP vs. FCS projected point to EPD, polar distance to projected hit;dR (cm);nmip", 80,-20,20, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipVtiledphi,"H2F_MixedPointProj_nmipVtiledphi","Mixed event nMIP vs. FCS projected point to EPD, angle difference to proj tile;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);

  //loaded += mHists->AddH2F(file,mH2F_PointProj_LowMult_nmipValldr,"H2F_PointProj_LowMult_nmipValldr","Low Mult nMIP vs. FCS projected point to EPD, polar distance to all other hits;dR (cm);nmip", 200,-100,100, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_PointProj_LowMult_nmipValldphi,"H2F_PointProj_LowMult_nmipValldphi","Low Mult nMIP vs. FCS projected point to EPD, angle difference to all other hits;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_MixedPointProj_LowMult_nmipValldr,"H2F_MixedPointProj_LowMult_nmipValldr","Low Mult Mixed event nMIP vs. FCS projected point to EPD, polar distance to projected hit;dR (cm);nmip", 200,-100,100, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_MixedPointProj_LowMult_nmipValldphi,"H2F_MixedPointProj_LowMult_nmipValldphi","Low Mult Mixed event nMIP vs. FCS projected point to EPD, angle difference to proj all;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);

  loaded += mHists->AddH1F(file,mH1F_AllPi0Mult,"H1F_AllPi0Mult","Pi0 Multiplicity with only an energy cut;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_AllPi0Zgg,"H1F_AllPi0Zgg","Zgg of all Pi0s with only an energy;Zgg;", 100,0,1);
  loaded += mHists->AddH2F(file,mH2F_AllPi0_etaVphi,"H1F_AllPi0_etaVphi","#eta vs. #phi of all Pi0s with only energy cut;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  loaded += mHists->AddH1F(file,mH1F_AllPi0En,"H1F_AllPi0En","Energy of all Pi0s with only an energy cut;Energy (GeV)", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_AllPi0Pt,"H1F_AllPi0Pt","Pt of Pi0s with only an energy cut;p_{T} (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_AllPi0Mass,"H1F_AllAllMass","Invariant mass of all point pair combinations;Invariant Mass (GeV);", 500,0,1);
  
  loaded += mHists->AddH1F(file,mH1F_NoEpdCutPi0Mult,"H1F_NoEpdCutPi0Mult","Pi0 Multiplicity all cuts except EPD nmip;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_NoEpdCutZgg,"H1F_NoEpdCutZgg","Zgg of all Pi0s, all cuts except EPD nmip;Zgg;", 100,0,1);
  loaded += mHists->AddH2F(file,mH2F_NoEpdCut_etaVphi,"H1F_NoEpdCut_etaVphi","#eta vs. #phi of all Pi0s, all cuts except EPD nmip;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  loaded += mHists->AddH1F(file,mH1F_NoEpdCutEn,"H1F_NoEpdCutEn","Energy of all Pi0s, all cuts except EPD nmip;Energy (GeV)", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_NoEpdCutPt,"H1F_NoEpdCutPt","Pt of all Pi0s, all cuts except EPD nmip;p_{T} (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_NoEpdCutAllMass,"H1F_NoEpdCutAllMass","Invariant mass all Pi0s, all cuts except EPD nmip;Invariant Mass (GeV);", 500,0,1);
  
  loaded += mHists->AddH1F(file,mH1F_EpdPhPi0Mult,"H1F_EpdPhPi0Mult","Point Multiplicity with an energy cut and EPD cut on photon;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_EpdPhZgg,"H1F_EpdPhZgg","Zgg of Pi0s using highest energy pairs and Epd Cut Photons;Z_{#gamma#gamma};", 100,0,1);
  loaded += mHists->AddH2F(file,mH2F_EpdPh_etaVphi,"H1F_EpdPh_etaVphi","#eta vs. #phi of Pi0s with all cuts and EPD photon cut on both points;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  loaded += mHists->AddH1F(file,mH1F_EpdPhEn,"H1F_EpdPhEn","Energy of Pi0s using highest energy pairs and Epd Cut Photons;Energy (GeV)", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_EpdPhPt,"H1F_EpdPhPt","Pt of Pi0s using highest energy pairs and Epd Cut Photons;Pt (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_EpdPhAllMass,"H1F_EpdPhAllMass","Invariant mass of all point pair combinations with Epd Cut Photons;Invariant Mass (GeV);", 500,0,1);

  loaded += mHists->AddH1F(file,mH1F_EpdChPi0Mult,"H1F_EpdChPi0Mult","Point Multiplicity with an energy cut and EPD cut on photon;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_EpdChZgg,"H1F_EpdChZgg","Zgg of Pi0s using highest energy pairs and Epd Cut Charged;Zgg;", 100,0,1);
  loaded += mHists->AddH2F(file,mH2F_EpdCh_etaVphi,"H1F_EpdCh_etaVphi","#eta vs. #phi of Pi0s with all cuts and EPD electron cut on both points;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  loaded += mHists->AddH1F(file,mH1F_EpdChEn,"H1F_EpdChEn","Energy of Pi0s using highest energy pairs and Epd Cut Charged;Energy (GeV)", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_EpdChPt,"H1F_EpdChPt","Pt of Pi0s using highest energy pairs and Epd Cut Charged;Pt (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_EpdChAllMass,"H1F_EpdChAllMass","Invariant mass of all point pair combinations with Epd Cut Charged;Invariant Mass (GeV);", 500,0,1); //This makes it such that this bin size is twice that of the 0,1 range with 500 bins

  loaded += mHists->AddH1F(file,mH1F_EpdSinglePhPi0Mult,"H1F_EpdSinglePhPi0Mult","Pi0 Multiplicity with all cuts and EPD cut on only one photon;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_EpdSinglePhZgg,"H1F_EpdSinglePhZgg","Z_{gg} of Pi0s with all cuts and EPD cut on only one photon;Z_{gg};", 100,0,1);
  loaded += mHists->AddH2F(file,mH2F_EpdSinglePh_etaVphi,"H1F_EpdSinglePh_etaVphi","#eta vs. #phi of Pi0s with all cuts and EPD cut on only one photon;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  loaded += mHists->AddH1F(file,mH1F_EpdSinglePhEn,"H1F_EpdSinglePhEn","Energy of Pi0s with all cuts and EPD cut on only one photon;Energy (GeV)", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_EpdSinglePhPt,"H1F_EpdSinglePhPt","p_{T} of Pi0s with all cuts and EPD cut on only one photon;#p_{T} (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_EpdSinglePhAllMass,"H1F_EpdSinglePhAllMass","Invariant mass of all point pair combinations after all cuts and Epd cut on single photon;Invariant Mass (GeV);", 500,0,1); //This makes it such that this bin size is twice that of the 0,1 range with 500 bins

  loaded += mHists->AddH1F(file,mH1F_EpdSingleChPi0Mult,"H1F_EpdSingleChPi0Mult","Pi0 Multiplicity with all cuts and EPD cut on only one photon;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_EpdSingleChZgg,"H1F_EpdSingleChZgg","Z_{gg} of Pi0s with all cuts and EPD cut on only one electron;Z_{gg};", 100,0,1);
  loaded += mHists->AddH2F(file,mH2F_EpdSingleCh_etaVphi,"H1F_EpdSingleCh_etaVphi","#eta vs. #phi of Pi0s with all cuts and EPD cut on only one electron;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  loaded += mHists->AddH1F(file,mH1F_EpdSingleChEn,"H1F_EpdSingleChEn","Energy of Pi0s with all cuts and EPD cut on only one electron;Energy (GeV)", 1000,0,200);
  loaded += mHists->AddH1F(file,mH1F_EpdSingleChPt,"H1F_EpdSingleChPt","p_{T} of Pi0s with all cuts and EPD cut on only one electron;#p_{T} (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_EpdSingleChAllMass,"H1F_EpdSingleChAllMass","Invariant mass of all point pair combinations after all cuts and Epd cut on single electron;Invariant Mass (GeV);", 500,0,1); //This makes it such that this bin size is twice that of the 0,1 range with 500 bins

  //loaded += mHists->AddH1F(file,SanityCheck,"SanityCheck","SanityCheck",500,0,1);

  if( mH1F_InvMassEpdCuts[0]==0 ){ mH1F_InvMassEpdCuts[0] = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_InvMassEpdCuts[0],NEPDCUTS,"H1F_InvMassEpdCuts_AllTrig","Different EPD NMIP cuts all triggers",500,0,1);
  if( mH1F_InvMassEpdCuts[1]==0 ){ mH1F_InvMassEpdCuts[1] = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_InvMassEpdCuts[1],NEPDCUTS,"H1F_InvMassEpdCuts_EmTrig","Different EPD NMIP cuts EM triggers",500,0,1);
  //loaded += mHists->AddH1FArr(file,mH1F_InvMassZggCuts,8,"H1F_InvMassZggCuts",500,0,1);
  //loaded += mHists->AddH1FArr(file,mH1F_InvMassPtCuts,8,"H1F_InvMassPtCuts",500,0,1);
  //loaded += mHists->AddH1F(file,mH1F_InvMassAllButEpdCut, "H1F_InvMassAllButEpdCut", "Invariant mass of two highest photon pairs with all cuts except EPD nmip",500,0,1);
  //loaded += mHists->AddH1F(file,mH1F_InvMassAllCutsEpdCh, "H1F_InvMassAllCutsEpdCh", "Invariant mass of two highest photon pairs with all cuts but charged epd nmip", 500,0,1);

  if( mH1F_InvMassAllCuts==0 ){ mH1F_InvMassAllCuts = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_InvMassAllCuts,5,"H1F_InvMassAllCuts","Invariant Mass of two photons after all cuts applied;M_{inv} (GeV/c^{2})", 500,0,1);
  if( mH1F_Pi0MultAllCuts==0 ){ mH1F_Pi0MultAllCuts = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_Pi0MultAllCuts,5,"H1F_Pi0MultAllCuts","Number of potential pi0s per event after all cuts;NGoodPi0", 20,0,20);
  //loaded += mHists->AddH2F(file,mH2F_AllCuts_Poi_yVx,"H2F_AllCuts_Poi_yVx","Point distribution y vs. x after all cuts", 400,-200,200, 300,-150,150);
  if( mH2F_AllCuts_Pi0_yVx==0 ){ mH2F_AllCuts_Pi0_yVx = new TObjArray(); }
  loaded += mHists->AddH2FArr(file,mH2F_AllCuts_Pi0_yVx,5,"H2F_AllCuts_Pi0_yVx","#pi^{0} distribution projected to FCS y vs. x after all cuts;x (cm);y (cm)", 400,-200,200, 300,-150,150);
  //loaded += mHists->AddH1F(file,mH1F_NFoundPhiBin,"H1F_NFoundPhiBin","Number of found phi bins",3,0,3);
  if( mH1F_AllCuts_xF==0 ){ mH1F_AllCuts_xF = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_AllCuts_xF,5,"H1F_AllCuts_xF","xF of pi0s after all cuts applied;x_{F}", 100,0,1);
  if( mH1F_AllCuts_xFZoom==0 ){ mH1F_AllCuts_xFZoom = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_AllCuts_xFZoom,5,"H1F_AllCuts_xFZoom","xF of pi0s after all cuts applied;x_{F}", 200,0,0.5);
  if( mH1F_AllCuts_Zgg==0 ){ mH1F_AllCuts_Zgg = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_AllCuts_Zgg,5,"H1F_AllCuts_Zgg","Z_{gg} of pi0s after all cuts applied;Z_{gg}",100,0,1);
  if( mH1F_AllCuts_Dgg==0 ){ mH1F_AllCuts_Dgg = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_AllCuts_Dgg,5,"H1F_AllCuts_Dgg","D_{gg} of pi0s after all cuts applied;D_{gg} (cm)",100,0,100);
  if( mH1F_AllCuts_Pi0En==0 ){ mH1F_AllCuts_Pi0En = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_AllCuts_Pi0En,5,"H1F_AllCuts_Pi0En","Energy of pi0s after all cuts applied;(GeV)", 200,0,200);
  if( mH2F_AllCuts_Pi0_massVen==0 ){ mH2F_AllCuts_Pi0_massVen = new TObjArray(); }
  loaded += mHists->AddH2FArr(file,mH2F_AllCuts_Pi0_massVen,5,"H2F_AllCuts_Pi0_massVen","Invariant mass of pi0 vs. pi0 Energy after all cuts applied;Energy (GeV);M_{inv} (GeV/c^{2})", 200,0,200, 500,0,1);
  if( mH2F_AllCuts_Pi0_xfVen==0 ){ mH2F_AllCuts_Pi0_xfVen = new TObjArray(); }
  loaded += mHists->AddH2FArr(file,mH2F_AllCuts_Pi0_xfVen,5,"H2F_AllCuts_Pi0_xfVen","Pi0 x_{F} vs. energy after all cuts applied;Energy (GeV);x_{F}", 200,0,200, 200,0,0.5);
  if( mH2F_AllCuts_Pi0_ptVeta==0 ){  mH2F_AllCuts_Pi0_ptVeta = new TObjArray(); }
  loaded += mHists->AddH2FArr(file,mH2F_AllCuts_Pi0_ptVeta,5,"H2F_AllCuts_Pi0_ptVeta","Pi0 p_{T} vs. eta after all cuts applied;#eta;p_{T}", 100,0,10, 100,0,50);
  
  //loaded += mHists->AddH1F(file,mH1F_AllCuts_Pi0Phi,"H1F_AllCuts_Pi0Phi","Phi of pi0s after all cuts applied", NPHIBIN,0,TMath::TwoPi());
  if( mH2F_AllCuts_Pi0_etaVphi==0 ){ mH2F_AllCuts_Pi0_etaVphi = new TObjArray(); }
  loaded += mHists->AddH2FArr(file,mH2F_AllCuts_Pi0_etaVphi,5,"H2F_AllCuts_Pi0_etaVphi","Eta vs. Phi distriubtion of pi0s;#phi;#eta", NPHIBIN,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0, 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_EpdNmip,"H2F_EpdNmip","EpdNmip;cluster;nmip", 2,0,2, 50,0,5);
			   
  //TString entitletext[NENERGYBIN] = { "En<=10", "10<En<=30", "30<En<=50", "50<En<=70", "70<En<=100", "En>100" };
  //TString phititletext[NPHIBIN] = { "0<=#phi<#frac{#pi}{4}", "#frac{#pi}{4}<=#phi<#frac{#pi}{2}", "#frac{#pi}{2}<=#phi<#frac{3#pi}{4}", "#frac{3#pi}{4}<=#phi<#pi" };
  double pi  = TMath::Pi();
  //Double_t xfbins[NXFBIN+1] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52};
  //Double_t xfbins[NXFBIN+1] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.3, 0.5};

  loaded += mHists->AddH2F(file,mH2F_NPi0Inc_xfVphi[0][0],"H2F_NPi0Inc_xfVphi_blue_up","Number of Pi0s in a given energy and phi bin for blue beam spin up 0.1<=M<=0.2", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Inc_xfVphi[0][1],"H2F_NPi0Inc_xfVphi_blue_down","Number of Pi0s in a given energy and phi bin for blue beam spin down 0.1<=M<=0.2", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Inc_xfVphi[1][0],"H2F_NPi0Inc_xfVphi_yellow_up","Number of Pi0s in a given energy and phi bin for yellow beam spin up 0.1<=M<=0.2", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Inc_xfVphi[1][1],"H2F_NPi0Inc_xfVphi_yellow_down","Number of Pi0s in a given energy and phi bin for yellow beam spin down 0.1<=M<=0.2", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );

  loaded += mHists->AddH2F(file,mH2F_NPi0Bg1_xfVphi[0][0],"H2F_NPi0Bg1_xfVphi_blue_up","Number of Pi0s in a given energy and phi bin for blue beam spin up 0.3<=M<=0.4", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Bg1_xfVphi[0][1],"H2F_NPi0Bg1_xfVphi_blue_down","Number of Pi0s in a given energy and phi bin for blue beam spin down 0.3<=M<=0.4", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Bg1_xfVphi[1][0],"H2F_NPi0Bg1_xfVphi_yellow_up","Number of Pi0s in a given energy and phi bin for yellow beam spin up 0.3<=M<=0.4", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Bg1_xfVphi[1][1],"H2F_NPi0Bg1_xfVphi_yellow_down","Number of Pi0s in a given energy and phi bin for yellow beam spin down 0.3<=M<=0.4", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );

  loaded += mHists->AddH2F(file,mH2F_NPi0Bg2_xfVphi[0][0],"H2F_NPi0Bg2_xfVphi_blue_up","Number of Pi0s in a given energy and phi bin for blue beam spin up 0.7<=M<=0.9", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Bg2_xfVphi[0][1],"H2F_NPi0Bg2_xfVphi_blue_down","Number of Pi0s in a given energy and phi bin for blue beam spin down 0.7<=M<=0.9", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Bg2_xfVphi[1][0],"H2F_NPi0Bg2_xfVphi_yellow_up","Number of Pi0s in a given energy and phi bin for yellow beam spin up 0.7<=M<=0.9", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );
  loaded += mHists->AddH2F(file,mH2F_NPi0Bg2_xfVphi[1][1],"H2F_NPi0Bg2_xfVphi_yellow_down","Number of Pi0s in a given energy and phi bin for yellow beam spin down 0.7<=M<=0.9", NPHIBIN,-pi/2.0,3.0*pi/2.0, NXFBIN,xfbins );

  //[January 29, 2025] > [How to create a TH3F with fixed and variable bin sizes](https://root-forum.cern.ch/t/how-to-make-th3-histograms-with-variable-bin-edges/38789)
  TAxis tmpphi(NPHIBIN,-pi/2.0,3.0*pi/2.0);
  Double_t phiedges[NPHIBIN+1] = {0};
  for( int i=0; i<25; ++i ){ phiedges[i] = tmpphi.GetBinUpEdge(i); }
  TAxis tmpmass(500,0,1);
  Double_t massedges[501] = {0};
  for( int i=0; i<501; ++i ){ massedges[i] = tmpmass.GetBinUpEdge(i); }
  loaded += mHists->AddH3F(file,mH3F_AllCutsInvMass_xfVphi,"H3F_InvMass_envVphi","Invariant mass after all cuts by phi and energy binning", NPHIBIN,phiedges, NXFBIN,xfbins, 500,massedges );
  /*
  for( int ebin=0; ebin<NENERGYBIN; ++ebin ){
    for( int phibin=0; phibin<NPHIBIN; ++phibin ){
      TString histname = "H1F_InvMassEn" + TString::Itoa(ebin,10) + "Phi" + TString::Itoa(phibin,10);
      TString histtitle = "Invariant Mass of Pi0s with " + entitletext[ebin] + " and " + phititletext[phibin] + ";InvMass (GeV);";
      loaded += mHists->AddH1F(file,mH1F_InvMassAllCutsByEnByPhi[ebin][phibin],histname.Data(),histtitle.Data(),500,0,1);
      histname = "H1F_NPi0En" + TString::Itoa(ebin,10) + "Phi" + TString::Itoa(phibin,10);
      histtitle = "Number of Pi0s with " + entitletext[ebin] + " and " + phititletext[phibin];
      loaded += mHists->AddH1F(file,mH1F_NPi0ByEnByPhi[ebin][phibin],histname.Data(),histtitle.Data(),4,0,4);
      mH1F_NPi0ByEnByPhi[ebin][phibin]->GetXaxis()->SetBinLabel(1,"NPi0UpAtPhi");
      mH1F_NPi0ByEnByPhi[ebin][phibin]->GetXaxis()->SetBinLabel(2,"NPi0UpAtPhiPlusPi");
      mH1F_NPi0ByEnByPhi[ebin][phibin]->GetXaxis()->SetBinLabel(3,"NPi0DownAtPhi");
      mH1F_NPi0ByEnByPhi[ebin][phibin]->GetXaxis()->SetBinLabel(4,"NPi0DownAtPhiPlusPi");
    }
  }
  */
  return loaded;
}

Int_t StMuFcsPi0TreeMaker::Init()
{
  if( mFilename.Length() == 0){ mFilename="StMuFcsPi0Ana.root"; } //Ensure a TFile is always created
  if( mHists==0 ){
    LOG_INFO << "StMuFcsPi0TreeMaker::Init() - No HistManager specified. Creating a new one with file name " << mFilename << ". Potential conflicts exist if a HistManager exists with same file name" << endm;
    mHists = new HistManager();
    mInternalHists = true;
    mFile_Output = mHists->InitFile(mFilename.Data(),"RECREATE");//new TFile(mFilename.Data(), "RECREATE");
  }
  else{
    mFile_Output = mHists->InitFile(); //No arguments just returns the internal file pointer
    mInternalHists = false;
  }
  mFile_Output->cd(); //File expected to be nonzero here if everything initialized correctly
  if( mTreeOnBitMap!=0 ){
    mPi0Tree     = new TTree("Pi0Tree","Tree with FcsPi0Candidate");
  }
  //These are still needed in Make so even if you are not writing the tree still need these objects
  mEvtInfo     = new FcsEventInfo();
  mPhArr       = new TClonesArray("FcsPhotonCandidate");
  mPi0Arr      = new TClonesArray("FcsPi0Candidate");
  //mMixedPhArr  = new TClonesArray("FcsPhotonCandidate");

  if( mTreeOnBitMap!=0 ){
    if( isEventOn() ){
      mPi0Tree->Branch("EventInfo","FcsEventInfo",&mEvtInfo);
      mPi0Tree->Branch("TriggerInfo",0,"NTrig/I:Trig[NTrig]/I");
      ((TLeaf*)mPi0Tree->GetBranch("TriggerInfo")->GetListOfLeaves()->At(0))->SetAddress(&mNTrig);
      ((TLeaf*)mPi0Tree->GetBranch("TriggerInfo")->GetListOfLeaves()->At(1))->SetAddress(&mTriggers);
    }
    if( isPhotonOn() ){ mPi0Tree->Branch("Photon",&mPhArr); }
    if( isPi0On() ){ mPi0Tree->Branch("Pi0",&mPi0Arr); }
  }
  
  mFcsTrigMap = (StFcsRun22TriggerMap*)GetMaker("fcsRun22TrigMap");
  
  if( mFcsTrigMap==0 ){ LOG_WARN << "StMuFcsPi0TreeMaker::Init() - No Trigger Map found" << endm; }
  else{ if( mFcsTrigMap->sizeOfTriggers()<=0 ){ LOG_WARN << "StMuFcsPi0TreeMaker::Init() - Trigger Map is empty" << endm; } }

  Int_t npoldata = ReadPolFile();
  //The 259 fills is special for my Run 22 text file after filtering out data from nonexistent fills
  if( npoldata!=259 ){ LOG_WARN << "Incorrect number of polarizations found in file 'Run22PolForJobs.txt' either because file is missing or file is improperly formatted" << endm; }

  UInt_t totalhists = this->LoadHists(0); //This is total of histograms loaded from a file not created. Don't use mFileOutput as you are not trying to load from #mFileOutput
  mHists->SetOwner(kTRUE);
  LOG_INFO << "StMuFcsPi0TreeMaker::Init() - Loaded " << totalhists << " histograms" << endm;

  return kStOk;
}

UInt_t StMuFcsPi0TreeMaker::LoadDataFromFile(TFile* file) //, TTree&* tree, FcsEventInfo&* evt,Int_t& ntrig, Int_t&* triggers,  TClonesArray&* pharr, TClonesArray&* pi0arr, TH1&* hist)
{
  if( file==0 || file->IsZombie() ){ std::cout << "LoadDataFromFile - ERROR:Unable to load from null file or zombie file" << std::endl; return 0; }
  mFile_Output = file;
  //if( tree!=0 ){ std::cout << "LoadDataFromFile - WARNING:Overwriting TTree pointer" << std::endl; }
  if( mPi0Tree!=0 ){ std::cout << "LoadDataFromFile - WARNING:Internal TTree pointer not zero must have been intialized elsewhere\n -> Deleting old data for new" << std::endl; delete mPi0Tree; mPi0Tree=0; }
  mPi0Tree = (TTree*)file->Get("Pi0Tree");
  //tree = mPi0Tree;
  if( mPi0Tree!=0 ){
    //Set event branches
    if( mEvtInfo!=0 ){ std::cout << "LoadDataFromFile - WARNING:Internal #FcsEventInfo not zero must have been intialized elsewhere\n -> Deleting old data for new" << std::endl; delete mEvtInfo; mEvtInfo=0; }
    if( isEventOn() ){
      if( mPi0Tree->Branch("EventInfo")!=0 && mPi0Tree->Branch("TriggerInfo")!=0 ){
	mEvtInfo     = new FcsEventInfo();
	mPi0Tree->SetBranchAddress("EventInfo",&mEvtInfo);
	//Set trigger branches
	((TLeaf*)mPi0Tree->GetBranch("TriggerInfo")->GetListOfLeaves()->At(0))->SetAddress(&mNTrig);
	((TLeaf*)mPi0Tree->GetBranch("TriggerInfo")->GetListOfLeaves()->At(1))->SetAddress(&mTriggers);
      }
      else{ std::cout << "LoadDataFromFile - WARNING:No \"EventInfo\" and no \"TriggerInfo\" branch found in mPi0Tree it could be that the tree was generated without this option." << std::endl;
      }
    }
    
    if( mPhArr!=0 ){ std::cout << "LoadDataFromFile - WARNING:Internal #FcsPhotonCandidate array not zero must have been intialized elsewhere\n -> Deleting old data for new" << std::endl; delete mPhArr; mPhArr=0; }
    if( isPhotonOn() ){
      if( mPi0Tree->Branch("Photon")!=0 ){
	mPhArr       = new TClonesArray("FcsPhotonCandidate");
	mPi0Tree->SetBranchAddress("Photon",&mPhArr);
      }
      else{
	std::cout << "LoadDataFromFile - WARNING:No \"Photon\" branch found in mPi0Tree it could be that the tree was generated without this option." << std::endl;
      }
    }
    
    if( mPi0Arr!=0 ){ std::cout << "LoadDataFromFile - WARNING:Internal #FcsPi0Candidate array not zero must have been intialized elsewhere\n -> Deleting old data for new" << std::endl; delete mPi0Arr; mPi0Arr=0; }
    if( isPi0On() ){
      if( mPi0Tree->Branch("Pi0")!=0 ){
	mPi0Arr      = new TClonesArray("FcsPi0Candidate");  
	mPi0Tree->SetBranchAddress("Pi0",&mPi0Arr);
      }
      else{ std::cout << "LoadDataFromFile - WARNING:No \"Pi0\" branch found in mPi0Tree it could be that the tree was generated without this option." << std::endl; }
    }
  }
  //else{ std::cout << "LoadDataFromFile - WARNING:Pi0Tree not found in file" << std::endl; }
  
  mFcsTrigMap = (StFcsRun22TriggerMap*)GetMaker("fcsRun22TrigMap"); //Don't use GetMaker since it doesn't use the right name when searching
  if( mFcsTrigMap==0 ){ std::cout << "StMuFcsPi0TreeMaker::LoadDataFromFile() - No Trigger Map found" << std::endl; }
  else{ if( mFcsTrigMap->sizeOfTriggers()<=0 ){ std::cout << "StMuFcsPi0TreeMaker::LoadDataFromFile() - Trigger Map is empty" << std::endl; } }
  
  if( mHists==0 ){ mHists = new HistManager(); }
  UInt_t totalhists = this->LoadHists(file);
  //std::cout << "TotalHistsLoaded: "<<totalhists << std::endl;
  //std::cout << "DUMB CHECK "<<mFcsTrigMap << std::endl;
  //std::string name(mFcsTrigMap->nameFromId(45,22349011));
  //std::cout << name << std::endl;
  
  return totalhists;
}

Int_t StMuFcsPi0TreeMaker::InitRun(int runnumber)
{
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  //mFcsDb->setDbAccess(0);
  if(!mFcsDb){
    LOG_ERROR << "StMuFcsPi0TreeMaker::InitRun - Failed to get StFcsDbMaker" << endm;
    return kStFatal;
  }

  mSpinDbMkr = static_cast<StSpinDbMaker*>(GetMaker("spinDb"));
  if( !mSpinDbMkr ){
    LOG_WARN << "StMuFcsPi0TreeMaker::InitRun - Could not find StSpinDbMaker named 'spinDb'" << endm;
  }
  else{
    if( !mSpinDbMkr->isValid() ){
      LOG_WARN << "StMuFcsPi0TreeMaker::InitRun - Found StSpinDbMaker but contains invalid data so will not use it" << endm;
      mSpinDbMkr=0;
      mH1F_RndmSpin->SetBinContent(1,1);
    }
    else{
      mH1F_RndmSpin->SetBinContent(2,1);
    }
  }
  if( !mEpdGeo ){ mEpdGeo = new StEpdGeom(); }
  else{
    LOG_ERROR << "StMuFcsPi0TreeMaker::InitRun - StEpdGeom Exists!" << endm;
    return kStFatal;
  }
  mEpdQaMkr = (StMuEpdRun22QaMaker*) GetMaker("FcsEpdRun22Qa");
  if( mEpdQaMkr==0 ){
    LOG_WARN << "StMuFcsPi0TreeMaker::InitRun - No StMuEpdRun22QaMaker found, EPD vertex information may be missing" << endm;
  }
  return kStOK;
}

//-----------------------
Int_t StMuFcsPi0TreeMaker::Finish()
{
  if( mFile_Output!=0 ){
    mFile_Output->cd();
    if( mPi0Tree!=0 ){ mPi0Tree->Write(); }
    mHists->Write();
  }
  return kStOK;
}

//----------------------
Int_t StMuFcsPi0TreeMaker::Make()
{
  Int_t result = Make_LoadEvent();
  if( result==kStErr ){ return kStErr; }
  result = Make_Polarization();
  if( result == kStErr ){ return kStErr; }
  result = Make_CheckAndSetTrig();
  if( !mValidTrigFound ){ mNTrig=0; return kStSkip; } //Reset trigger array size before going to next event
  if( result!=kStOk ){ return result; }
  result = Make_SpinInfo();
  if( result!=kStOk ){ return result; }

  result = this->Make_GetEpdColl();
  if( result==kStWarn ){ return kStWarn; }

  result = this->Make_VertexInfo();
  if( result!=kStOk ){ return kStErr; }

  result = this->Make_FillFcsClusPoint();
  result = this->Make_CheckAndSetEpdHit();
  result = this->Make_PointPairs(mPi0Arr);
  result = this->Make_TssaAna(mPi0Arr);

  if( mPi0Tree!=0 ){ mPi0Tree->Fill(); }
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_LoadEvent()
{
  mMuDstMkr = (StMuDstMaker*)GetInputDS("MuDst");
  if( mMuDstMkr==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make_LoadEvent - !MuDstMkr" <<endm; return kStErr; }
  mMuDst = mMuDstMkr->muDst();
  if( mMuDst==0 ){ LOG_ERROR << "StMuFcsPi0TreeMaker::Make_LoadEvent - !MuDst" << endm; return kStErr; }
  mMuEvent = mMuDst->event();
  if( mMuEvent==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make_LoadEvent - !MuEvent" <<endm; return kStErr; }
  mTrigData = mMuEvent->triggerData();
  if( mTrigData==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make_LoadEvent - !TrigData" <<endm; return kStErr; }
  mRunInfo = &(mMuEvent->runInfo());
  if( mRunInfo==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make_LoadEvent - !RunInfo" <<endm; return kStErr; }
  
  mH1D_Entries->Fill(0); //This is just counting valid make calls (i.e. increment bin 1 by 1)

  //FcsEventInfo* evtinfo = (FcsEventInfo*) mEvtInfo->ConstructedAt(0);  
  mEvtInfo->mRunTime         = mRunInfo->beamFillNumber(StBeamDirection::east);    //using yellow beam
  mEvtInfo->mRunNum          = mMuEvent->runNumber();
  mEvtInfo->mFill            = mMuEvent->eventInfo().time();
  mEvtInfo->mEvent           = mMuEvent->eventId();
  mEvtInfo->mBx48Id          = mTrigData->bunchId48Bit();
  mEvtInfo->mBx7Id           = mTrigData->bunchId7Bit();
  mEvtInfo->mTofMultiplicity = mTrigData->tofMultiplicity();
  //std::cout << "|runtime:"<<mEvtInfo->mRunTime << "|runnum:"<<mEvtInfo->mRunNum << "|event:"<<mEvtInfo->mEvent << std::endl;

  
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_Polarization()
{
  Int_t fillnum = mRunInfo->beamFillNumber(StBeamDirection::east);
  Int_t evttime = mMuEvent->eventInfo().time();
  PolData* poldat = 0;
  if( mPolarizationData.find(fillnum) == mPolarizationData.end() ){
    //If no polarization data for this fill then set it equal to zero do it doesn't contribute to the sum
    poldat = new PolData();
    poldat->mFillNum = fillnum;
    poldat->mStartTime = evttime;
    poldat->mBlueP0 = 0;
    poldat->mBluedPdT = 0;
    poldat->mYellowdPdT = 0;
    poldat->mYellowP0 = 0;
    mPolarizationData[fillnum] = poldat;
  }
  else{
    poldat = mPolarizationData[fillnum];
  }

  //Divide by 3600 to convert seconds to hours since dP/dT is in %/hour
  Double_t timeelapsed = Double_t(evttime-poldat->mStartTime)/3600.0;
  //Sum because dPdT is already negative
  Double_t polblue      = timeelapsed * poldat->mBluedPdT + poldat->mBlueP0;
  Double_t polyellow    = timeelapsed * poldat->mYellowdPdT + poldat->mYellowP0;
  Double_t polblueerr   = sqrt(timeelapsed*poldat->mBlueErrdPdT*timeelapsed*poldat->mBlueErrdPdT + poldat->mBlueErrP0*poldat->mBlueErrP0);
  Double_t polyellowerr = sqrt(timeelapsed*poldat->mYellowErrdPdT*timeelapsed*poldat->mYellowErrdPdT + poldat->mYellowErrP0*poldat->mYellowErrP0);
  mH1D_BluePol->Fill(polblue);
  mH1D_YellowPol->Fill(polyellow);
  mH1D_BluePolErr->Fill(polblueerr);
  mH1D_YellowPolErr->Fill(polyellowerr);
  //std::cout  << " + "<<"|eventnum:"<< mH1D_Entries->GetBinContent(1)  <<"|fillnum:"<<fillnum << "|evttime:"<<evttime << "|polblue:"<<polblue << "|polyellow:"<<polyellow << "|totpolblue:"<<mH1D_Entries->GetBinContent(2) << "|totpolyellow:"<<mH1D_Entries->GetBinContent(3) << std::endl;
  
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_CheckAndSetTrig()
{
  //Filter with Trigger Information first
  if( mIgnoreTrig ){ mValidTrigFound = true; return kStOk; } //If ignoring triggers then set mValidTrigFound to true so event is not skipped
  StMuTriggerIdCollection* TrigMuColl = &(mMuEvent->triggerIdCollection());
  if( !TrigMuColl ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::FillEventInfo - !TrigMuColl" <<endl; return kStErr; }
  const StTriggerId& trgIDs = TrigMuColl->nominal();
  Int_t ntrig = trgIDs.triggerIds().size();
  for( Int_t i=0; i<ntrig; ++i ){
    unsigned int trig = trgIDs.triggerId(i);
    mTriggers[i] = trig;
    if( mFcsTrigMap!=0 ){
      std::string thistrig = mFcsTrigMap->nameFromId(trig,mMuEvent->runNumber());      
      if( !mIgnoreTrig ){
	for( unsigned int j=0; j<mTargetTrig.size(); ++j ){
	  if( mTargetTrig.at(j)==thistrig ){
	    mNTrig++;
	    mValidTrigFound=true;
	    mH1F_Triggers->Fill( thistrig.c_str(), 1);
	  }
	}
      }
      if( thistrig=="fcsEM0" || thistrig=="fcsEM1" || thistrig=="fcsEM2" || thistrig=="fcsEM3"
	  || thistrig=="fcsEM0_tpc" || thistrig=="fcsEM1_tpc" || thistrig=="fcsEM2_tpc" || thistrig=="fcsEM3_tpc" )
	{
	  if( thistrig=="fcsEM0" || thistrig=="fcsEM0_tpc" ){ mTrigEm0 = 1; }
	  if( thistrig=="fcsEM1" || thistrig=="fcsEM1_tpc" ){ mTrigEm1 = 2; }
	  if( thistrig=="fcsEM2" || thistrig=="fcsEM2_tpc" ){ mTrigEm2 = 3; }
	  if( thistrig=="fcsEM3" || thistrig=="fcsEM3_tpc" ){ mTrigEm3 = 4; }
	  mEmTrigFound = true;
	}
    }
    else{ mNTrig = ntrig; }
  }
  /* Debugging why some events were skipped. The events that were skipped were fcsHad triggered events because I didn't turn them on in the runMuDst.C macro
  if( !mValidTrigFound ){
    mNTrig=0; std::cout << "No ValidTrigFound" << std::endl;
    for( Int_t i=0; i<ntrig; ++i ){
      unsigned int trig = trgIDs.triggerId(i);
      std::cout << "|trig:"<<trig << std::endl;
    }
    return kStSkip;
    }*/
  return kStOk;
}
  
Int_t StMuFcsPi0TreeMaker::Make_SpinInfo()
{
  //Spin information
  if( mSpinDbMkr==0 ){
    Double_t rndm = mSpinRndm.Rndm(); //random number between 0 and 1
    //The numbers below are chosen for their bit representation as described and correspond to the source polarization. Need to flip to convert it to STAR polarization direction because of the Siberian Snakes; i.e. '+' -> '-' and '-' -> '+'
    if( rndm<=0.25 ){ mEvtInfo->mSpin = 5; }                    //Bits 0101 is B+ and Y+
    else if( 0.25<rndm && rndm<=0.5 ){ mEvtInfo->mSpin = 6; }   //Bits 0110 is B+ and Y-
    else if( 0.5<rndm && rndm<=0.75 ){ mEvtInfo->mSpin = 9; }   //Bits 1001 is B- and Y+
    else{ mEvtInfo->mSpin = 10; }                               //Bits 1010 is B- and Y-
  }
  else{
    mEvtInfo->mSpin = mSpinDbMkr->spin4usingBX7( mTrigData->bunchId7Bit() ); //This is also source polarization
  }
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_VertexInfo()
{
  //Vertex Information
  mEvtInfo->mVpdVz = -999;
  if( mMuDst->btofHeader() ){ mEvtInfo->mVpdVz = mMuDst->btofHeader()->vpdVz(); }
  //@[April 7, 2021] > No Slewing correction for BBC yet, see StFmsJetMaker2015 in BrightSTAR??
  mEvtInfo->mBbcTacDiff = mTrigData->bbcTimeDifference() - 4096; //subtract 4096 since 0 means bad event and distribution is Gaussian around 4096
  mEvtInfo->mBbcVz = -999;
  if( fabs(mEvtInfo->mBbcTacDiff)>1.e-6 ){ mEvtInfo->mBbcVz = mEvtInfo->mBbcTacDiff * -0.2475; } //0.2475 = 0.0165*30/2x
  //StZdcTriggerDetector& zdc = mMuEvent->zdcTriggerDetector();
  //std::cout <<"|ZdcV:"<<zdc.vertexZ() << std::endl;
  mEvtInfo->mZdcVz =  mTrigData->zdcVertexZ();

  if( mEpdQaMkr!=0 ){
    mEvtInfo->mEpdTacEarlyE = mEpdQaMkr->epdTacEarlyE();
    mEvtInfo->mEpdTacEarlyW = mEpdQaMkr->epdTacEarlyW();
    mEvtInfo->mEpdAvgE = mEpdQaMkr->epdTacAvgE();
    mEvtInfo->mEpdAvgW = mEpdQaMkr->epdTacAvgW();
    mEvtInfo->mEpdVz = mEpdQaMkr->epdVertex();
  }
  else{
    mEvtInfo->mEpdVz = -999;
  }

  if( mEvtInfo->mVpdVz > -998 )     { mFoundVertex = 1; mUseVertex = mEvtInfo->mVpdVz; }
  else if( mEvtInfo->mEpdVz > -998 ){ mFoundVertex = 2; mUseVertex = mEvtInfo->mEpdVz; }
  else if( mEvtInfo->mBbcVz > -998 ){ mFoundVertex = 4; mUseVertex = mEvtInfo->mBbcVz; }
  else{ mFoundVertex = 0; mUseVertex = 0; } //If no vertex found use 0
  mEvtInfo->mFoundVertex = mFoundVertex;
  mH2F_foundVvertex->Fill(mUseVertex,mFoundVertex);

  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_GetEpdColl()
{
  //Get EPD collection and/or hits
  mMuEpdHits = 0;
  mEpdColl = 0;
  mMuEpdHits = mMuDst->epdHits();
  if( mMuEpdHits!=0 ){ if( mMuEpdHits->GetEntriesFast()==0 ){mMuEpdHits=0;} }//If mMuEpdHits is not zero but has no hits set it to zero so rest of code processes from StEpdHitMaker
  if( mMuEpdHits==0 ){ LOG_INFO << "StMuFcsPi0TreeMaker::Make - No MuEPD hits" << endm;
    mEpdHitMkr = (StEpdHitMaker*)GetMaker("epdHit");
    if( mEpdHitMkr==0 ){ LOG_WARN << "StMuFcsPi0TreeMaker::Make - No StEpdHitMaker(\"epdHit\")" << endm; return kStWarn; }
    else{ mEpdColl = mEpdHitMkr->GetEpdCollection(); }
    if( mEpdColl==0 ){ LOG_WARN << "StMuEpdRun22QaMaker::FillFcsInfo - No Epd hit information found" << endm; mEpdHitMkr=0; return kStWarn; }//Set the hit maker back to zero so it can be used as a check that the epd collection doesn't exist
  }
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_FillFcsClusPoint()
{
  //Fcs Collection
  mMuFcsColl = mMuDst->muFcsCollection();
  if (!mMuFcsColl) { LOG_ERROR << "StMuFcsPi0TreeMaker::Make did not find MuFcsCollection" << endm; return kStErr; }

  //TClonesArray* hits = mMuFcsColl->getHitArray();
  //if( hits==0 ){ LOG_INFO << "StMuFcsPi0TreeMaker::FillFcsInfo - No FCS hits" << endm; }
  TClonesArray* clusters = mMuFcsColl->getClusterArray();
  //if( clusters==0 ){ LOG_INFO << "StMuFcsPi0TreeMaker::FillFcsInfo - No FCS clusters" << endm; }
  if( clusters==0 ){ std::cout << "StMuFcsPi0TreeMaker::FillFcsInfo - No FCS clusters" << std::endl; }
  TClonesArray* points = mMuFcsColl->getPointArray();
  //if( points==0 ){ LOG_INFO << "StMuFcsPi0TreeMaker::FillFcsInfo - No FCS points" << endm; }
  if( points==0 ){ std::cout << "StMuFcsPi0TreeMaker::FillFcsInfo - No FCS points" << std::endl; }
  
  //std::cout << "|hits:"<<hits << "|clusters:"<<clusters << "|points:"<<points << std::endl;
  Int_t ncandidates = 0;
  if( clusters!=0 ){
    for( UInt_t idet=0; idet<=kFcsEcalSouthDetId; ++idet ){
      //This is to do a fiducial volume cut on the clusters and points, since points and clusters store local coordinates use row, column space designation
      float mincolumncut = 1.0; //Column 1 will give x value between 0 and 1 so want greater than 1
      float minrowcut = 1.0;    //Row 1 will give xy value between 0 and 1 so want greater than 1
      float maxcolumncut = mFcsDb->nColumn(idet) - 1; //subtract 1 from max column to get rid of the edge towers
      float maxrowcut = mFcsDb->nRow(idet) - 1;       //subtract 1 from max rows to get rid of the edge towers
      //std::cout << "+ |idet:"<<idet << "|maxdet:"<<kFcsNDet;
      unsigned int nc = mMuFcsColl->numberOfClusters(idet);
      unsigned int iclus=mMuFcsColl->indexOfFirstCluster(idet);
      nc += iclus;
      for( ; iclus<nc; ++iclus){
	StMuFcsCluster* clu = (StMuFcsCluster*)clusters->At(iclus);
	float iclu_x = clu->x();
	float iclu_y = clu->y();
	float iclu_energy = clu->energy();
	if( iclu_energy<mEnCut ){ continue; }
	if( !(mincolumncut<=iclu_x && iclu_x<=maxcolumncut) ){ continue; }
	if( !(minrowcut<=iclu_y && iclu_y<=maxrowcut) ){ continue; }

	StThreeVectorD iclu_pos = mFcsDb->getStarXYZfromColumnRow( idet, iclu_x, iclu_y );
	StLorentzVectorD iclu_p = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, 0 );

	FcsPhotonCandidate* ph = (FcsPhotonCandidate*) mPhArr->ConstructedAt(ncandidates++);
	ph->mFromCluster = true;
	ph->mDetId = idet;
	ph->mX = iclu_pos[0];
	ph->mY = iclu_pos[1];
	ph->mZ = iclu_pos[2];

	ph->mEn = iclu_energy;
	ph->mPxRaw = iclu_p.px();
	ph->mPyRaw = iclu_p.py();
	ph->mPzRaw = iclu_p.pz();
	mH1F_ClusterEnergy->Fill(iclu_energy);

	//std::cout << "Cluster|detid:"<<ph->mDetId << "|mX:"<<ph->mX << "|mY:"<<ph->mY << "|mZ:"<<ph->mZ << std::endl;

	if( mFoundVertex > 0 ){
	  StLorentzVectorD iclu_p_withz = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, mUseVertex );
	  ph->mPxVert = iclu_p_withz.px();
	  ph->mPyVert = iclu_p_withz.py();
	  ph->mPzVert = iclu_p_withz.pz();
	}
	else{
	  ph->mPxVert = 0;
	  ph->mPyVert = 0;
	  ph->mPzVert = 0;
	}
      }
    }
  }

  //std::cout << "===== EventId:"<< mEvtInfo->mEvent <<" =====" << std::endl;
  mEvtInfo->mClusterSize = ncandidates;
  Int_t clustersize = ncandidates; //local copy of mEvtInfo->mClusterSize
  
  if( points!=0 ){
    for( UInt_t idet=0; idet<=kFcsEcalSouthDetId; ++idet ){
      //This is to do a fidicul volume cut on the clusters and points, since points and clusters store local coordinates use row, column space designation
      float mincolumncut = 1.0; //Column 1 will give x value between 0 and 1 so want greater than 1
      float minrowcut = 1.0;    //Row 1 will give xy value between 0 and 1 so want greater than 1
      float maxcolumncut = mFcsDb->nColumn(idet) - 1; //subtract 1 from max column to get rid of the edge towers
      float maxrowcut = mFcsDb->nRow(idet) - 1;       //subtract 1 from max rows to get rid of the edge towers
      unsigned int np = mMuFcsColl->numberOfPoints(idet);
      unsigned int ipoint=mMuFcsColl->indexOfFirstPoint(idet);
      np += ipoint;
      for( ; ipoint<np; ++ipoint ){
	StMuFcsPoint* point = (StMuFcsPoint*)points->At(ipoint);
	float ipoi_x = point->x();
	float ipoi_y = point->y();
	float ipoi_energy = point->energy();
	if( ipoi_energy<mEnCut ){ continue; }
	if( !(mincolumncut<=ipoi_x && ipoi_x<=maxcolumncut) ){ continue; }
	if( !(minrowcut<=ipoi_y && ipoi_y<=maxrowcut) ){ continue; }

	StThreeVectorD ipoi_pos = mFcsDb->getStarXYZfromColumnRow( idet, ipoi_x, ipoi_y );
	StLorentzVectorD ipoi_p = mFcsDb->getLorentzVector(ipoi_pos, ipoi_energy, 0);

	FcsPhotonCandidate* ph = (FcsPhotonCandidate*) mPhArr->ConstructedAt(ncandidates++);
	ph->mFromCluster = false;
	ph->mDetId = idet;
	ph->mX = ipoi_pos[0];
	ph->mY = ipoi_pos[1];
	ph->mZ = ipoi_pos[2];
	mH2F_PhotonHeatMap->Fill(ipoi_pos[0],ipoi_pos[1]);
	if( 78.2<=ipoi_energy && ipoi_energy<79.4 ){ mH2F_PhotonHeatMapB->Fill(ipoi_pos[0],ipoi_pos[1]); } //Region with spike
	if( 76.8<=ipoi_energy && ipoi_energy<78.0 ){ mH2F_PhotonHeatMapG->Fill(ipoi_pos[0],ipoi_pos[1]); } //Region near spike

	ph->mEn = ipoi_energy;
	ph->mPxRaw = ipoi_p.px();
	ph->mPyRaw = ipoi_p.py();
	ph->mPzRaw = ipoi_p.pz();
	mH1F_PointEnergy->Fill(ipoi_energy);

	if( mFoundVertex > 0 ){ 	StLorentzVectorD ipoi_p_withz = mFcsDb->getLorentzVector( ipoi_pos, ipoi_energy, mUseVertex );
	  ph->mPxVert = ipoi_p_withz.px();
	  ph->mPyVert = ipoi_p_withz.py();
	  ph->mPzVert = ipoi_p_withz.pz();
	}
	else{
	  ph->mPxVert = 0;
	  ph->mPyVert = 0;
	  ph->mPzVert = 0;
	}
	//std::cout << "Point|detid:"<<ph->mDetId << "|mX:"<<ph->mX << "|mY:"<<ph->mY << "|mZ:"<<ph->mZ << std::endl;	
      }//i point
    }//fcs dets
  }

  //std::cout << "|ncandidates:"<<ncandidates <<"|clustersize:"<<clustersize <<"|Size:"<<mPhArr->GetEntriesFast() << std::endl;
  Int_t npoints = ncandidates - clustersize; //Don't need to add 1 since including clustersize but not ncandidates
  mH1F_PointMult->Fill(npoints);
  mPhArr->Sort(); //Since this is properly sorted with clusters showing up first clustersize is unchanged. Also sorts by energy
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_CheckAndSetEpdHit()
{
  //Check photon candidates if they have any hits in the EPD. Use a separate loop so that this information could be used in the pi0 checking loop if needed. In future may also want to check against FCS preshower (EPD) hits
  Int_t ntotal = mPhArr->GetEntriesFast();
  Int_t npoints = ntotal - mEvtInfo->mClusterSize;
  unsigned int nepdhits = 0;
  StSPtrVecEpdHit* epdhits = 0;
  if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
  else if( mEpdColl!=0 ){
    epdhits = &(mEpdColl->epdHits());
    nepdhits = epdhits->size();
  }
  else{ LOG_ERROR << "StMuEpdRun22QaMaker::FillEpdinfo() - If you see this error then there is a bug that is setting EPD hits improperly" << endm; return kStErr; }
  Int_t nepdwesthits = 0;
  for( Int_t iph = 0; iph<ntotal; ++iph ){
    //std::cout << "|iph:"<<iph << "|iphnew:"<<iph-noldhits << std::endl;
    FcsPhotonCandidate* ph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(iph);
    if( ph==0 ){ std::cout << "==========I=CANNOT=BE=ZERO==========" << std::endl; return kStErr; }
    std::vector<Double_t> epdproj = ProjectToEpd(ph->mX,ph->mY,ph->mZ,mUseVertex);
    
    mH2F_EpdProjHitMap->Fill( epdproj.at(0),epdproj.at(1) );
    if( mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh ){ mH2F_EpdProjHitMap_Vcut->Fill(epdproj.at(0),epdproj.at(1)); }
    //std::cout << " + |phx:"<<ph->mX << "|phy:"<<ph->mY << "|phz:"<<ph->mZ << "|v:"<<mUseVertex << std::endl;
    //std::cout << " + |epdx:"<<epdproj.at(0) << "|epdy:"<<epdproj.at(1) << "|epdz:"<<epdproj.at(2) << std::endl;
    //std::cout << " ** |iph:"<< iph-noldhits << std::endl;
    CheckInsideEpdTile(ph,epdproj.at(0),epdproj.at(1));  //Function that will check which EPD tiles photon candidate overlaps with and sets the appropriate variables for it
    /*
      if( ph->mEpdMatch[0]==0 ){ //If no intersection found it would be -1 so now check all the CCW adjacencies
      std::cout << "     ** |iph:"<<iph-noldhits <<"|projx:"<<epdproj.at(0) << "|projy:"<<epdproj.at(1) << "|nmip:"<< ph->mEpdHitNmip[0] << "|epdkey:"<<ph->mEpdMatch[0] << std::endl;
      }*/
    
    //loop over all hits and if an nmip value exists set for the point
    StMuEpdHit* muepdhit = 0;
    StEpdHit* epdhit = 0;
    for(unsigned int i=0; i<nepdhits; ++i ){
      if( mMuEpdHits!=0 ){ muepdhit = (StMuEpdHit*)mMuEpdHits->UncheckedAt(i); } //To match similar in StMuDstMaker->epdHit(int i)
      else if( epdhits!=0 ){ epdhit = (StEpdHit*)((*epdhits)[i]); }
      else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A VERY SERIOUS BUG IN THE CODE" << endm; return kStErr; } 
      //std::cout << "|i:"<<i << "|muepdhit:"<<muepdhit << "|epdhit:"<<epdhit << std::endl;
      int ew    = muepdhit!=0 ? muepdhit->side()    : epdhit->side();      //east=-1, west=1
      if( ew==-1 ){ continue; }
      if( iph==0){ ++nepdwesthits; }
      int epdpp = muepdhit!=0 ? muepdhit->position(): epdhit->position();  //Supersector runs [1,12]
      int epdtt = muepdhit!=0 ? muepdhit->tile()    : epdhit->tile();      //Tile number [1,31]
      //int adc = muepdhit!=0 ? muepdhit->adc() : epdhit->adc();
      float nmip = muepdhit!=0 ? muepdhit->nMIP(): epdhit->nMIP();         //The ADC value of the hit divided by the MIP peak position; e.g. if nmip==1 then adc value sits at the MIP peak
      /*TVector3 epdhitxyz = mEpdGeo->TileCenter(epdpp,epdtt,ew);
      Double_t dx = epdproj.at(0)-epdhitxyz[0];
      Double_t dy = epdproj.at(1)-epdhitxyz[1];
      double rpoint = sqrt(epdproj.at(0)*epdproj.at(0) + epdproj.at(1)*epdproj.at(1));
      double rhit = sqrt(epdhitxyz[0]*epdhitxyz[0] + epdhitxyz[1]*epdhitxyz[1]);
      Double_t phipoint = TMath::ATan2(epdproj.at(1),epdproj.at(0));
      Double_t phihit = TMath::ATan2(epdhitxyz[1],epdhitxyz[0]);
      Double_t diffphi = phipoint-phihit;
      if( diffphi>TMath::Pi() ){ diffphi = diffphi - TMath::Pi(); }
      if( diffphi<(-1.0*TMath::Pi()) ){ diffphi = diffphi + TMath::Pi(); }
      */
      //std::cout << "|epdpp:"<<epdpp <<"|epdtt:"<<epdtt <<"|nmip:"<<nmip << std::endl;
      //std::cout << "|epdz:"<<epdhitxyz[2] << std::endl;
      if( ! ph->mFromCluster ){
	//if( mTrigEm2==3 && mTrigEm0<0 && mTrigEm1<0 ){
	if( mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh ){
	  //mH2F_PointProj_nmipValldx->Fill(dx,nmip);
	  //mH2F_PointProj_nmipValldy->Fill(dy,nmip);
	  //mH2F_PointProj_nmipValldr->Fill(rpoint-rhit,nmip);
	  //mH2F_PointProj_nmipValldphi->Fill(diffphi,nmip);
	  if( npoints<=5 ){
	    //mH2F_PointProj_LowMult_nmipValldr->Fill(rpoint-rhit,nmip);
	    //mH2F_PointProj_LowMult_nmipValldphi->Fill(diffphi,nmip);
	  }
	  if( ph->mEpdMatch[0] == (100*epdpp+epdtt)  ){
	    //mH2F_PointProj_nmipVtiledx->Fill(dx,nmip);
	    //mH2F_PointProj_nmipVtiledy->Fill(dy,nmip);
	    //mH2F_PointProj_nmipVtiledr->Fill(rpoint-rhit,nmip);
	    //mH2F_PointProj_nmipVtiledphi->Fill(diffphi,nmip);
	    ph->mEpdHitNmip[0] = nmip;
	  }
	}
      }
    }
    mH2F_EpdNmip->Fill(ph->mFromCluster,ph->mEpdHitNmip[0]);
  }
  //std::cout << "|nold:"<<noldhits << "|nnew:"<<nnewhits << "|ntotal:"<<ntotal << "|oldvert:"<<mOldVertex << "|newvert:"<<mUseVertex<< "|nepdhits:"<<nepdwesthits <<"|noldpoints:"<<mNOldPoints << "|npoints:"<<npoints << std::endl;

  //std::cout << "|clustersize:"<<clustersize << "|ncandidates:"<<ncandidates << "|npoints:"<<npoints << std::endl;
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_PointPairs(TClonesArray* pointpairs)
{
  Int_t clustersize = mEvtInfo->mClusterSize;  
  Int_t npi0candidate = 0;
  //Filling cluster pi0s and cluster photon/elecron epd nmip cut. For clusters only store best pair to speed up code
  for( Int_t ic = 0; ic<clustersize; ++ic ){
    FcsPhotonCandidate* iclus = (FcsPhotonCandidate*) mPhArr->UncheckedAt(ic);
    if( !(iclus->mFromCluster) ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
    //std::cout << "|ic:"<<ic << std::endl;
    //std::cout << "  + ";
    //iclus->Print();
    if( ic==(clustersize-1) ){ continue; }

    ///if( iclus->mEpdHitNmip>-0.1){ //Only include candidates who have their nmip value set
    //if( iclus->mEpdHitNmip<mEpdNmipCut ){ goodclusphotonsidx.emplace_back(ic); }
    //else{ goodcluselectronsidx.emplace_back(ic); }
    //}
    for( Int_t jc=ic+1; jc<clustersize; jc++ ){
      FcsPhotonCandidate* jclus = (FcsPhotonCandidate*) mPhArr->UncheckedAt(jc);
      if( !(jclus->mFromCluster) ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
      TLorentzVector pi0Vert_LV = iclus->lvVert() + jclus->lvVert();
      if( ic==0 && jc==ic+1 ){ //Since we have a sorted photon array highest two energies are the first two entries
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) pointpairs->ConstructedAt(npi0candidate++);
	pi0c->mFromCluster = true;
	pi0c->mFromPh = 0;
	pi0c->mPhoton1Idx = ic;
	pi0c->mPhoton2Idx = jc;
	
	pi0c->mPx = pi0Vert_LV.Px();
	pi0c->mPy = pi0Vert_LV.Py();
	pi0c->mPz = pi0Vert_LV.Pz();
	pi0c->mEn = pi0Vert_LV.E();
	
	pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
	pi0c->mDgg     = FcsPi0Candidate::dgg(*iclus,*jclus);
	pi0c->mZgg     = FcsPi0Candidate::zgg(*iclus,*jclus);
	pi0c->mAlpha   = FcsPi0Candidate::alpha(*iclus,*jclus);
	pi0c->mInvMass = pi0Vert_LV.Mag();
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|clustermass:"<<pi0c->mInvMass <<  std::endl;
	break;
      }
      /*
	else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_AllPointPairMass->Fill(pi0Vert_LV.Mag());
	}*/
    }
  }

  //Filling point pi0s and cluster photon/elecron epd nmip cut
  //std::cout << "===== EventId:"<< mEvtInfo->mEvent <<" =====" << std::endl;
  FcsPhotonCandidate* firstphotoncut[NEPDCUTS] = {0};
  FcsPhotonCandidate* secondphotoncut[NEPDCUTS] = {0};
  Int_t n_noepdproj = 0;
  Int_t n_noepdproj_vcut = 0;
  for( Int_t ip = clustersize; ip<mPhArr->GetEntriesFast(); ++ip ){
    FcsPhotonCandidate* ipoi = (FcsPhotonCandidate*) mPhArr->UncheckedAt(ip);
    if( ipoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
    //std::cout << "|ip:"<<ip << std::endl;
    //std::cout << "  + ";
    //ipoi->Print();

    if( ipoi->mEpdHitNmip[0]>-0.1){ //Only include candidates who have their nmip value set
      for( short i=0; i<NEPDCUTS; ++i ){
	if( ipoi->mEpdHitNmip[0] < 0.2+0.1*static_cast<double>(i) ){
	  if( firstphotoncut[i]==0 ){ firstphotoncut[i]=ipoi; }
	  else{ if( secondphotoncut[i]==0 ){ secondphotoncut[i]=ipoi; } }
	}
      }
    }
    else{
      ++n_noepdproj;
      if( mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh ){ ++n_noepdproj_vcut; }
    }
    
    if( ip==(mPhArr->GetEntriesFast()-1) ){ continue; }
    for( Int_t jp=ip+1; jp<mPhArr->GetEntriesFast(); ++jp ){
      FcsPhotonCandidate* jpoi = (FcsPhotonCandidate*) mPhArr->UncheckedAt(jp);
      if( jpoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
      TLorentzVector pi0Vert_LV = ipoi->lvVert() + jpoi->lvVert();
      //if( ip==clustersize && jp==ip+1 ){ //Since we have a sorted photon array highest two energies are the first two entries
      FcsPi0Candidate* pi0c = (FcsPi0Candidate*) pointpairs->ConstructedAt(npi0candidate++);
      pi0c->mFromCluster = false;
      pi0c->mFromPh = 0;
      pi0c->mPhoton1Idx = ip;
      pi0c->mPhoton2Idx = jp;
      
      pi0c->mPx = pi0Vert_LV.Px();
      pi0c->mPy = pi0Vert_LV.Py();
      pi0c->mPz = pi0Vert_LV.Pz();
      pi0c->mEn = pi0Vert_LV.E();
      
      pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
      pi0c->mDgg     = FcsPi0Candidate::dgg(*ipoi,*jpoi);
      pi0c->mZgg     = FcsPi0Candidate::zgg(*ipoi,*jpoi);
      pi0c->mAlpha   = FcsPi0Candidate::alpha(*ipoi,*jpoi);
      pi0c->mInvMass = pi0Vert_LV.Mag();
      //std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
      mH2F_Energy_ph1Vph2->Fill(ipoi->mEn,jpoi->mEn);
    }
  }
  mH1F_NBadEpdProj->Fill(n_noepdproj);
  mH1F_NBadEpdProjVcut->Fill(n_noepdproj_vcut);

  //Handle the epd cuts
  for( short i=0; i<NEPDCUTS; ++i ){
    if( firstphotoncut[i]!=0 && secondphotoncut[i]!=0 ){
      TLorentzVector pi0Vert_LV = firstphotoncut[i]->lvVert() + secondphotoncut[i]->lvVert();
      if( mEmTrigFound ){ ((TH1*)mH1F_InvMassEpdCuts[1]->UncheckedAt(i))->Fill(pi0Vert_LV.Mag()); }
      ((TH1*)mH1F_InvMassEpdCuts[0]->UncheckedAt(i))->Fill(pi0Vert_LV.Mag());
    }
  }
  return kStOk;
}

Int_t StMuFcsPi0TreeMaker::Make_TssaAna(TClonesArray* pointpairs)
{
  std::map<Int_t,PolData*>::iterator  politr =   mPolarizationData.find(mEvtInfo->mFill);
  PolData* poldat = 0;
  if( politr!=mPolarizationData.end() ){ poldat = politr->second; }
  else{ LOG_WARN << "No polariztion data found for fill "<< mEvtInfo->mFill << endm; return kStSkip; }
  
  Int_t nallpi0 = 0;
  Int_t npi0noepdcut = 0 ;
  Int_t ngoodpi0s = 0;
  Int_t ngoodsingleph = 0;
  Int_t ngoodbothph = 0;
  Int_t ngoodsinglech = 0;
  Int_t ngoodbothch = 0;
  //std::cout << "NPi0s:"<<pointpairs->GetEntriesFast() << std::endl;
  short emtrig[5] = {0, mTrigEm0, mTrigEm1, mTrigEm2, mTrigEm3 };
  for( int i=0; i<mPi0Arr->GetEntriesFast(); ++i ){
    FcsPi0Candidate* pi0 = (FcsPi0Candidate*)pointpairs->At(i);
    if( pi0==0 ){ continue; }
    if( pi0->mFromCluster ){ /*std::cout << "StMuFcsPi0TreeMaker::Make() - Not a point - "<<pi0->mFromCluster<< std::endl;*/ continue; }
    ++nallpi0;
    Double_t pi0en = pi0->mEn;
    TLorentzVector pi0_lv = pi0->lv();
    Double_t phi = pi0_lv.Phi(); //Range of this phi is -pi to pi
    Double_t mpi = -TMath::Pi();
    if( mpi<=phi && phi<mpi/2.0 ){ phi += TMath::TwoPi(); } //Since my binning goes from -pi/2 to 3pi/2 need to add 2pi to angles in the region from [-pi,-pi/2)
    //pi0->Print();

    //mH2F_BestPi0HeatMap->Fill();
    mH1F_AllPi0Zgg->Fill(pi0->mZgg);
    mH2F_AllPi0_etaVphi->Fill(phi,pi0_lv.Eta());
    //mH1F_AllPi0Phi->Fill(pi0c->phi());
    //mH1F_AllPi0Eta->Fill(pi0c->mEta);
    mH1F_AllPi0En->Fill(pi0en);
    mH1F_AllPi0Pt->Fill(pi0->pt());
    mH1F_AllPi0Mass->Fill(pi0->mInvMass);

    if( ! (mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh) ){ /*std::cout << " StMuFcsPi0TreeMaker::Make() - Failed vertex cut:"<<mUseVertex << std::endl;*/ continue; }
    if( pi0->mZgg>0.7     ){ /*std::cout << "StMuFcsPi0TreeMaker::Make() - Failed Zgg:"<< pi0->mZgg << std::endl;*/ continue; }
    //Add pt cut based on trigger
    bool exceedtrigpt = false;
    //Float_t trigptthr = -1;
    //std::string trigname = "";
    if( !mIgnoreTrig ){
      if( !mEmTrigFound ){ continue; }
      if( mFcsTrigMap!=0 ){
	for( Int_t itrig=0; itrig<mNTrig; ++itrig ){
	  Float_t ptthres = mFcsTrigMap->GetPtThr(mTriggers[itrig]);
	  std::string thistrig = mFcsTrigMap->nameFromId(mTriggers[itrig],mMuEvent->runNumber());
	  if( pi0->pt()>=ptthres ){ exceedtrigpt=true; /*trigptthr=ptthres; trigname=thistrig;*/ }
	}
      }
    }
    else{ exceedtrigpt = true; }
    if( !exceedtrigpt ){ continue; }

    Double_t pi0xf = static_cast<Double_t>(pi0->mPz) / static_cast<Double_t>(poldat->mBeamEn);
    //if( pi0xf<0.01 && pi0en>5){ std::cout << "  + 4MOM|("<<pi0en<<","<<pi0->mPx<<","<<pi0->mPy<<","<<pi0->mPz << ")"<<"|phi:"<<phi<<"eta:"<<pi0->eta()<<"|beamen:"<<poldat->mBeamEn << "|xf:"<<pi0xf << "|vert:"<<mUseVertex << std::endl; }
    Double_t pi0mass = (Double_t)pi0->mass();
    
    //epd photon cut, mFromPh can only be -1,0,1 for less than epd cut (neutral), no epd cut, greater than epd cut (charged); respectively
    //if( pi0->mFromPh != -1 ){ continue; } //Check to force mFromPh to always be -1 (neutral)
    if( pi0->mFromPh==0 ){
      mH1F_NoEpdCutZgg->Fill(pi0->mZgg);
      mH2F_NoEpdCut_etaVphi->Fill(phi,pi0_lv.Eta());
      mH1F_NoEpdCutEn->Fill(pi0en);
      mH1F_NoEpdCutPt->Fill(pi0->pt());
      mH1F_NoEpdCutAllMass->Fill(pi0mass);   //All but EpdPh cut
      ++npi0noepdcut;
      FcsPhotonCandidate* ph1 = (FcsPhotonCandidate*)mPhArr->UncheckedAt(pi0->mPhoton1Idx);
      FcsPhotonCandidate* ph2 = (FcsPhotonCandidate*)mPhArr->UncheckedAt(pi0->mPhoton2Idx);
      if( ph1->mEpdHitNmip[0]>-0.1 && ph2->mEpdHitNmip[0]>-0.1 ){ //Only include candidates who have their nmip value set
	if( ph1->mEpdHitNmip[0]<mEpdNmipCut || ph2->mEpdHitNmip[0]<mEpdNmipCut ){ //This is negation of mFromPh==1
	  //Fill hisotrams with only one photon satisfying epd cut criteria for nuetral particles (photons)
	  mH1F_EpdSinglePhZgg->Fill(pi0->mZgg);
	  mH2F_EpdSinglePh_etaVphi->Fill(phi,pi0_lv.Eta());
	  mH1F_EpdSinglePhEn->Fill(pi0en);
	  mH1F_EpdSinglePhPt->Fill(pi0_lv.Pt());
	  mH1F_EpdSinglePhAllMass->Fill(pi0mass);
	  ++ngoodsingleph;
	}
	if( ph1->mEpdHitNmip[0]>=mEpdNmipCut || ph2->mEpdHitNmip[0]>=mEpdNmipCut ){ //This is negation of mFromPh==-1
	  //Fill hisotrams with only one photon satisfying epd cut criteria for charged particle
	  mH1F_EpdSingleChZgg->Fill(pi0->mZgg);
	  mH2F_EpdSingleCh_etaVphi->Fill(phi,pi0_lv.Eta());
	  mH1F_EpdSingleChEn->Fill(pi0en);
	  mH1F_EpdSingleChPt->Fill(pi0_lv.Pt());
	  mH1F_EpdSingleChAllMass->Fill(pi0mass);
	  ++ngoodsinglech;
	}
	if( ph1->mEpdHitNmip[0]<mEpdNmipCut && ph2->mEpdHitNmip[0]<mEpdNmipCut ){
	  pi0->mFromPh = -1;  //Label this pi0 as coming from two photons so can store it and analyze
	  mH1F_EpdPhZgg->Fill(pi0->mZgg);
	  mH2F_EpdPh_etaVphi->Fill(phi,pi0_lv.Eta());
	  mH1F_EpdPhEn->Fill(pi0en);
	  mH1F_EpdPhPt->Fill(pi0_lv.Pt());
	  mH1F_EpdPhAllMass->Fill(pi0mass);
	  ++ngoodbothph;
	}
	if( ph1->mEpdHitNmip[0]>=mEpdNmipCut && ph2->mEpdHitNmip[0]>=mEpdNmipCut ){
	  pi0->mFromPh = 1;  //Label this pi0 as coming from two electrons
	  mH1F_EpdChZgg->Fill(pi0->mZgg);
	  mH2F_EpdCh_etaVphi->Fill(phi,pi0_lv.Eta());
	  mH1F_EpdChEn->Fill(pi0en);
	  mH1F_EpdChPt->Fill(pi0_lv.Pt());
	  mH1F_EpdChAllMass->Fill(pi0mass);
	  ++ngoodbothch;
	}
      }
    }
    if( pi0->mFromPh>=0   ){ /*std::cout << "StMuFcsPi0TreeMaker::Make() - Failed photon cut - "<<pi0->mFromPh<< std::endl;*/ continue; } //Got rid of extra photon loops so everything is now mFromPh==0. Do this check so only both nmip requirements is stored for the A_N analysis
    //std::cout << "StMuFcsPi0TreeMaker::Make() - Passed all cuts!" << std::endl;
    ++ngoodpi0s;
    //std::cout << " + |Ntrig:"<<mNTrig << "|trigname:"<<trigname << "|trigpt:"<<trigptthr;
    //pi0->Print();
    /*
    FcsPhotonCandidate* ph1 = mPhArr->UncheckedAt(pi0->mPhoton1Idx);
    FcsPhotonCandidate* ph2 = mPhArr->UncheckedAt(pi0->mPhoton2Idx);
    mH2F_AllCuts_PoiX_1V2->(ph1->mX,ph2->mX);
    mH2F_AllCuts_PoiY_1V2->(ph1->mY,ph2->mY);
    */

    //Get Pi0 projection to Ecal front face @[Jan 27, 20254] > Get rid of don't need to project Just use TLorentzVector and use the TLorentzVector::phi() to get the phi
    //std::cout << "  + |pi0->pt():"<<pi0->pt() << "|pi0_lv.Pt():"<<pi0_lv.Pt() << std::endl;
    //std::cout << "    - |pi0->mPz:"<<pi0->mPz << "|pi0_lv.Pz():"<<pi0_lv.Pz() << std::endl;
    //std::cout << "    - |pi0->eta():"<<pi0->eta() << "|pi0_lv.Eta():"<<pi0_lv.Eta() << std::endl;
    double pi0_momentum[3] = {pi0->mPx,pi0->mPy,pi0->mPz};
    double pi0_vertex[3] = {0,0,mUseVertex};
    int det = 0; //North side if negative px
    //South side for positive px, set px==0 to south side since fcs planes would intersect in that case
    if( pi0_momentum[2]>=0 && pi0_momentum[0]>=0 ){ det=1; }
    if( pi0_momentum[2]<0  && pi0_momentum[0]<0  ){ det=1; }
    StThreeVectorD pi0_xyz = mFcsDb->projectLine(det,pi0_momentum,pi0_vertex,0); //Project to front face of FCS
    
    //@[Jan 27, 2025] > (Keep energy binnng and check energy plot for binning) Make a phi histogram one for blue (up & down) and yellow (up & down) from -Pi/2 to 3/2 Pi with whatever binning. bin 1 which is bottom most on left, with nbins/2, TH1F* mhphi[2][2] = {0}; //[blue,yellow] [up,down], Move to Ana code [mhasym[2] //[blue,yellow beam]. mhphi[0][0]->bin(1)*mhphi[0][1]->bin(1+nbin/2)]
    for( short i=0; i<5; ++i ){
      if( emtrig[i]>=0 ){
	((TH2*) mH2F_AllCuts_Pi0_yVx->UncheckedAt(i))->Fill(pi0_xyz.x(),pi0_xyz.y());
	((TH1*) mH1F_AllCuts_xF->UncheckedAt(i))->Fill( pi0xf );
	((TH1*) mH1F_AllCuts_xFZoom->UncheckedAt(i))->Fill( pi0xf );
	((TH1*) mH1F_AllCuts_Zgg->UncheckedAt(i))->Fill( pi0->mZgg );
	((TH1*) mH1F_AllCuts_Dgg->UncheckedAt(i))->Fill( pi0->mDgg );
	((TH1*) mH1F_AllCuts_Pi0En->UncheckedAt(i))->Fill( pi0en );
	((TH2*) mH2F_AllCuts_Pi0_massVen->UncheckedAt(i))->Fill( pi0en,pi0mass );
	((TH2*) mH2F_AllCuts_Pi0_xfVen->UncheckedAt(i))->Fill( pi0en,pi0xf );
	((TH2*) mH2F_AllCuts_Pi0_ptVeta->UncheckedAt(i))->Fill( pi0_lv.Eta(),pi0_lv.Pt() );
	((TH2*) mH2F_AllCuts_Pi0_etaVphi->UncheckedAt(i))->Fill( phi,pi0_lv.Eta() );
	((TH1*) mH1F_InvMassAllCuts->UncheckedAt(i))->Fill(pi0mass);
      }
    }
    //mH1F_InvMassAllCutsByEnByPhi[enbin][phibin]->Fill(pi0->mass());
    ((TH3*) mH3F_AllCutsInvMass_xfVphi)->Fill(phi,pi0xf,pi0mass);    
    //short nfoundphibin = 0;
    if( 0.1<=pi0mass && pi0mass<=0.2){
      if( mEvtInfo->BlueSpin()==1    ){ mH2F_NPi0Inc_xfVphi[0][0]->Fill(phi,pi0xf); }
      if( mEvtInfo->BlueSpin()==-1   ){ mH2F_NPi0Inc_xfVphi[0][1]->Fill(phi,pi0xf); }
      if( mEvtInfo->YellowSpin()==1  ){ mH2F_NPi0Inc_xfVphi[1][0]->Fill(phi,pi0xf); }
      if( mEvtInfo->YellowSpin()==-1 ){ mH2F_NPi0Inc_xfVphi[1][1]->Fill(phi,pi0xf); }
    }
    if( 0.3<=pi0mass && pi0mass<=0.4){
      if( mEvtInfo->BlueSpin()==1    ){ mH2F_NPi0Bg1_xfVphi[0][0]->Fill(phi,pi0xf); }
      if( mEvtInfo->BlueSpin()==-1   ){ mH2F_NPi0Bg1_xfVphi[0][1]->Fill(phi,pi0xf); }
      if( mEvtInfo->YellowSpin()==1  ){ mH2F_NPi0Bg1_xfVphi[1][0]->Fill(phi,pi0xf); }
      if( mEvtInfo->YellowSpin()==-1 ){ mH2F_NPi0Bg1_xfVphi[1][1]->Fill(phi,pi0xf); }
    }
    if( 0.7<=pi0mass && pi0mass<=0.9){
      if( mEvtInfo->BlueSpin()==1    ){ mH2F_NPi0Bg2_xfVphi[0][0]->Fill(phi,pi0xf); }
      if( mEvtInfo->BlueSpin()==-1   ){ mH2F_NPi0Bg2_xfVphi[0][1]->Fill(phi,pi0xf); }
      if( mEvtInfo->YellowSpin()==1  ){ mH2F_NPi0Bg2_xfVphi[1][0]->Fill(phi,pi0xf); }
      if( mEvtInfo->YellowSpin()==-1 ){ mH2F_NPi0Bg2_xfVphi[1][1]->Fill(phi,pi0xf); }
    }
  }
  mH1F_NoEpdCutPi0Mult->Fill(npi0noepdcut);
  mH1F_AllPi0Mult->Fill(nallpi0);
  mH1F_EpdSinglePhPi0Mult->Fill(ngoodsingleph);
  mH1F_EpdSingleChPi0Mult->Fill(ngoodsinglech);
  mH1F_EpdPhPi0Mult->Fill(ngoodbothph);
  mH1F_EpdChPi0Mult->Fill(ngoodbothch);
  for( short i=0; i<5; ++i ){
    //std::cout << " - |emtrig["<<i<<"]:"<<emtrig[i] <<"|n:"<<ngoodpi0s << std::endl;
    if( emtrig[i]>=0 ){
      ((TH1*) mH1F_Pi0MultAllCuts->UncheckedAt(i))->Fill(ngoodpi0s);
    }
  }
  //std::cout << "NGoodPi0s:"<<ngoodpi0s << std::endl;

  return kStOk;
}

void StMuFcsPi0TreeMaker::Clear(Option_t* option)
{
  //Reset Variables
  mFoundVertex = 0;
  mUseVertex = -999.0;
  
  mValidTrigFound = false;
  mEmTrigFound = false;
  mTrigEm0 = -1;
  mTrigEm1 = -1;
  mTrigEm2 = -1;
  mTrigEm3 = -1;
  
  mEvtInfo->Clear("C");
  mNTrig = 0; //Since ROOT only writes up to the size of mNTrig then only need to reset this back to zero and next loop will overwrite array as neccessary
  mPhArr->Clear("C");
  mPi0Arr->Clear("C");
  
  return;
}

/*void StMuFcsPi0TreeMaker::AnalyzePi0s()
{

  return;

  }*/

void StMuFcsPi0TreeMaker::Print(Option_t* opt) const
{
  TString option(opt);
  option.ToLower();
  if( option.Contains("a") ){ option = "etgp"; }
  if( mEvtInfo!=0 && option.Contains("e") ){ mEvtInfo->Print(); }
  if( option.Contains("t") ){
    std::cout << "## Trigger Information|NTrig:"<<mNTrig << std::endl;
    for( int i=0; i<mNTrig; ++i ){
      std::cout << " + |TrigId:"<<mTriggers[i];
      if( mFcsTrigMap!=0 && mEvtInfo!=0 ){ std::cout << "|TrigName:"<< mFcsTrigMap->nameFromId(mTriggers[i],mEvtInfo->mRunNum); }
      std::cout << std::endl;
    }
  }
  if( option.Contains("g") ){
    std::cout << "## Photon Information|Size:"<<mPhArr->GetEntriesFast() << std::endl;
    for( int i=0; i<mPhArr->GetEntriesFast(); ++i ){
      std::cout << " + ";
      mPhArr->At(i)->Print();
    }
  }
  if( option.Contains("p") ){
    std::cout << "## Pi0 Information|Size:"<<mPi0Arr->GetEntriesFast() << std::endl;
    for( int i=0; i<mPi0Arr->GetEntriesFast(); ++i ){
      std::cout << " + ";
      mPi0Arr->At(i)->Print();
    }
  }
}

std::vector<Double_t> StMuFcsPi0TreeMaker::ProjectToEpd(Double_t xfcs, Double_t yfcs, Double_t zfcs, Double_t zvertex)
{
  //Assume x,y=0 at vertex so only need zvertex as origin and is the initial point for the direction
  Double_t linedirection[3] = {xfcs,yfcs,zfcs-zvertex}; //This is the direction vector from the origin to the fcs point (FcsXYZ-Vertex)
  Double_t EpdZ = 375.0;   //Need a point on epd plane so picking x,y=0 is valid and formula below reflects this. Also only care about West EPD which is along positive z-axis
  //Tiles are parallel to z-axis so normal vector is {0,0,1}. The formula below reflects this
  //Solution of intersection of line and plane where line has direction {xdir,ydir,zdir}*t and starts at {0,0,zvertex} and a plane with a normal that points along the z-axis and has a point on the plane at {0,0,EpdZ}; "t" is the free parameter in the parametric equation of the line.
  double tintersection = (EpdZ-zvertex) / (linedirection[2]);
  std::vector<Double_t> intersection;
  intersection.emplace_back(linedirection[0]*tintersection);
  intersection.emplace_back(linedirection[1]*tintersection);
  intersection.emplace_back(linedirection[2]*tintersection+zvertex);
  return intersection;
}

void StMuFcsPi0TreeMaker::PaintEventQa(TCanvas* canv,  const char* savename) const
{
  canv->Clear();

  canv->Divide(2,2);
  canv->cd(1);
  mH1D_Entries->Draw("hist e");
  canv->cd(2)->SetLogy();
  mH1F_Triggers->Draw("hist e");
  canv->cd(3)->SetLogz();
  mH2F_foundVvertex->Draw("colz");
  canv->cd(4);
  mH1F_RndmSpin->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintPolarization(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(2,2);
  canv->cd(1);
  mH1D_BluePol->Draw("hist e");
  canv->cd(2);
  mH1D_BluePolErr->Draw("hist e");
  canv->cd(3);
  mH1D_YellowPol->Draw("hist e");
  canv->cd(4);
  mH1D_YellowPolErr->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintPhotonQa(TCanvas* canv, const char* savename)   const
{
  canv->Clear();
  
  canv->Divide(3,3);
  canv->cd(1)->SetLogz();
  mH2F_PhotonHeatMap->Draw("colz");
  canv->cd(2)->SetLogz();
  mH2F_EpdProjHitMap->Draw("colz");
  canv->cd(3)->SetLogz();
  mH2F_EpdProjHitMap_Vcut->Draw("colz");
  canv->cd(4)->SetLogz();
  mH2F_EpdNmip->Draw("colz");
  canv->cd(5)->SetLogy();
  mH1F_ClusterEnergy->Draw("hist e");
  canv->cd(6)->SetLogy();
  mH1F_PointMult->Draw("hist e");
  canv->cd(7)->SetLogy();
  mH1F_PointEnergy->Draw("hist e");
  canv->cd(8)->SetLogz();
  mH2F_Energy_ph1Vph2->Draw("colz");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintAllPi0(TCanvas* canv,  const char* savename)  const
{
  canv->Clear();
  
  canv->Divide(3,2);
  canv->cd(1)->SetLogy();
  mH1F_AllPi0Mult->Draw("hist e");
  canv->cd(2);
  mH1F_AllPi0Zgg->Draw("hist e");
  canv->cd(3);
  //mH1F_AllPi0Phi->Draw("hist e");
  //canv->cd(4);
  //mH1F_AllPi0Eta->Draw("hist e");
  mH2F_AllPi0_etaVphi->Draw("colz");
  canv->cd(4)->SetLogy();
  mH1F_AllPi0En->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_AllPi0Pt->Draw("hist e");
  canv->cd(6);
  mH1F_AllPi0Mass->Draw("hist e");
  //canv->cd(8);
  //mH1F_AllPointPairMass->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintNoEpdCut(TCanvas* canv,  const char* savename)  const
{
  canv->Clear();
  
  canv->Divide(3,2);
  canv->cd(1)->SetLogy();
  mH1F_NoEpdCutPi0Mult->Draw("hist e");
  canv->cd(2);
  mH1F_NoEpdCutZgg->Draw("hist e");
  canv->cd(3);
  mH2F_NoEpdCut_etaVphi->Draw("colz");
  canv->cd(4)->SetLogy();
  mH1F_NoEpdCutEn->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_NoEpdCutPt->Draw("hist e");
  canv->cd(6);
  mH1F_NoEpdCutAllMass->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdPhPi0(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,2);
  canv->cd(1)->SetLogy();
  mH1F_EpdPhPi0Mult->Draw("hist e");
  canv->cd(2);
  mH1F_EpdPhZgg->Draw("hist e");
  canv->cd(3);
  //mH1F_EpdPhPhi->Draw("hist e");
  //canv->cd(4);
  //mH1F_EpdPhEta->Draw("hist e");
  mH2F_EpdPh_etaVphi->Draw("colz");
  canv->cd(4)->SetLogy();
  mH1F_EpdPhEn->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_EpdPhPt->Draw("hist e");
  canv->cd(6);
  mH1F_EpdPhAllMass->Draw("hist e");
  //canv->cd(8);
  //mH1F_EpdPhAllPoints->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdChPi0(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,2);
  canv->cd(1)->SetLogy();
  mH1F_EpdChPi0Mult->Draw("hist e");
  canv->cd(2);
  mH1F_EpdChZgg->Draw("hist e");
  canv->cd(3);
  //mH1F_EpdChPhi->Draw("hist e");
  //canv->cd(4);
  //mH1F_EpdChEta->Draw("hist e");
  mH2F_EpdCh_etaVphi->Draw("colz");
  canv->cd(4)->SetLogy();
  mH1F_EpdChEn->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_EpdChPt->Draw("hist e");
  canv->cd(6);
  mH1F_EpdChAllMass->Draw("hist e");
  //canv->cd(8);
  //mH1F_EpdChAllPoints->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdSinglePh(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,2);
  canv->cd(1)->SetLogy();
  mH1F_EpdSinglePhPi0Mult->Draw("hist e");
  canv->cd(2);
  mH1F_EpdSinglePhZgg->Draw("hist e");
  canv->cd(3);
  mH2F_EpdSinglePh_etaVphi->Draw("colz");
  canv->cd(4)->SetLogy();
  mH1F_EpdSinglePhEn->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_EpdSinglePhPt->Draw("hist e");
  canv->cd(6);
  mH1F_EpdSinglePhAllMass->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdSingleCh(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,2);
  canv->cd(1)->SetLogy();
  mH1F_EpdSingleChPi0Mult->Draw("hist e");
  canv->cd(2);
  mH1F_EpdSingleChZgg->Draw("hist e");
  canv->cd(3);
  mH2F_EpdSingleCh_etaVphi->Draw("colz");
  canv->cd(4)->SetLogy();
  mH1F_EpdSingleChEn->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_EpdSingleChPt->Draw("hist e");
  canv->cd(6);
  mH1F_EpdSingleChAllMass->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintPi0Overlap(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,3);

  //canv->cd(1);
  canv->cd(1)->SetLogy();
  TLegend* legpad1 = new TLegend(0.5,0.5,0.93,0.93,"","nbNDC");
  TH1* h1pi0mult = mH1F_NoEpdCutPi0Mult->DrawCopy("hist e");
  h1pi0mult->SetStats(0);
  h1pi0mult->SetLineColor(kBlack);
  legpad1->AddEntry(h1pi0mult,"NoEpdCut","fle");
  //AddHistStatsOneline(legpad1,h1pi0mult,"NoEpdCut");
  TH1* h1phmult = mH1F_EpdPhPi0Mult->DrawCopy("hist e same");
  h1phmult->SetLineColor(kBlue);
  std::stringstream ss_legname;
  ss_legname << "EpdNmip<"<<mEpdNmipCut;
  //AddHistStatsOneline(legpad1,h1phmult,ss_legname.str().c_str());
  legpad1->AddEntry(h1phmult,ss_legname.str().c_str(),"fle");
  TH1* h1singlephmult = mH1F_EpdSinglePhPi0Mult->DrawCopy("hist e same");
  h1singlephmult->SetLineColor(kRed);
  ss_legname.str("");
  ss_legname << "EpdNmip<"<<mEpdNmipCut << " single point";
  //AddHistStatsOneline(legpad1,h1singlephmult,ss_legname.str().c_str());
  legpad1->AddEntry(h1singlephmult,ss_legname.str().c_str(),"fle");
  TH1* h1singlechmult = mH1F_EpdSingleChPi0Mult->DrawCopy("hist e same");
  h1singlechmult->SetLineColor(kGreen+2);
  ss_legname.str("");
  ss_legname << "EpdNmip>="<<mEpdNmipCut << " single point";
  //AddHistStatsOneline(legpad1,h1singlechmult,ss_legname.str().c_str());
  legpad1->AddEntry(h1singlechmult,ss_legname.str().c_str(),"fle");
  legpad1->Draw();
  
  //canv->cd(2);
  canv->cd(2);//->SetLogy();
  TH1* h1pi0zgg = mH1F_NoEpdCutZgg->DrawCopy("hist e");
  h1pi0zgg->SetLineColor(kBlack);
  TH1* h1phzgg = mH1F_EpdPhZgg->DrawCopy("hist e same");
  h1phzgg->SetLineColor(kBlue);
  TH1* h1singlephzgg = mH1F_EpdSinglePhZgg->DrawCopy("hist e same");
  h1singlephzgg->SetLineColor(kRed);
  TH1* h1singlechzgg  = mH1F_EpdSingleChZgg->DrawCopy("hist e same");
  h1singlechzgg->SetLineColor(kGreen+2);
  
  canv->cd(3);
  //canv->cd(3)->SetLogy();
  TH1* h1_pi0_phi = ((TH2*)mH2F_NoEpdCut_etaVphi)->ProjectionX("h1_pi0_phi");
  //TH1* h1_allpi0_phi_norm = h1_allpi0_phi->DrawNormalized("hist e");
  h1_pi0_phi->SetLineColor(kBlack);
  h1_pi0_phi->Draw("hist e");
  TH1* h1_epdph_phi = ((TH2*)mH2F_EpdPh_etaVphi)->ProjectionX("h1_epdph_phi");
  //TH1* h1_epdph_phi_norm = h1_epdph_phi->DrawNormalized("hist e same");
  h1_epdph_phi->SetLineColor(kBlue);
  h1_epdph_phi->Draw("hist e same");
  TH1* h1_epdsingleph_phi = ((TH2*)mH2F_EpdSinglePh_etaVphi)->ProjectionX("h1_epdsingleph_phi");
  //TH1* h1_epdph_phi_norm = h1_epdph_phi->DrawNormalized("hist e same");
  h1_epdsingleph_phi->SetLineColor(kRed);
  h1_epdsingleph_phi->Draw("hist e same");
  TH1* h1_epdsinglech_phi = ((TH2*)mH2F_EpdSingleCh_etaVphi)->ProjectionX("h1_epdsinglech_phi");
  //TH1* h1_epdch_phi_norm = h1_epdch_phi->DrawNormalized("hist e same");
  h1_epdsinglech_phi->SetLineColor(kGreen+2);
  h1_epdsinglech_phi->Draw("hist e same");
  //Cleanup unneeded histograms that are copies of the normalized ones and not drawn on the canvas and therefor won't be cleaned up
  //delete h1_allpi0_phi; h1_allpi0_phi=0;
  //delete h1_epdph_phi;  h1_epdph_phi=0;
  //delete h1_epdch_phi;  h1_epdch_phi=0;
  
  canv->cd(4);
  //canv->cd(4)->SetLogy();
  TH1* h1_pi0_eta = ((TH2*)mH2F_NoEpdCut_etaVphi)->ProjectionY("h1_pi0_eta");
  //TH1* h1_pi0_eta_norm = h1_allpi0_eta->DrawNormalized("hist e");
  h1_pi0_eta->SetLineColor(kBlack);
  h1_pi0_eta->Draw("hist e");
  TH1* h1_epdph_eta = ((TH2*)mH2F_EpdPh_etaVphi)->ProjectionY("h1_epdph_eta");
  //TH1* h1_epdph_eta_norm = h1_epdph_eta->DrawNormalized("hist e same");
  h1_epdph_eta->SetLineColor(kRed);
  h1_epdph_eta->Draw("hist e same");
  TH1* h1_epdsingleph_eta = ((TH2*)mH2F_EpdPh_etaVphi)->ProjectionY("h1_epdsingleph_eta");
  //TH1* h1_epdch_eta_norm = h1_epdch_eta->DrawNormalized("hist e same");
  h1_epdsingleph_eta->SetLineColor(kRed);
  h1_epdsingleph_eta->Draw("hist e same");
  TH1* h1_epdsinglech_eta = ((TH2*)mH2F_EpdSingleCh_etaVphi)->ProjectionY("h1_epdsinglech_eta");
  //TH1* h1_epdch_eta_norm = h1_epdch_eta->DrawNormalized("hist e same");
  h1_epdsinglech_eta->SetLineColor(kGreen+2);
  h1_epdsinglech_eta->Draw("hist e same");
  //Cleanup unneeded histograms that are copies of the normalized ones and not drawn on the canvas and therefor won't be cleaned up
    //delete h1_allpi0_eta; h1_allpi0_eta=0;
    //delete h1_epdph_eta;  h1_epdph_eta=0;
    //delete h1_epdch_eta;  h1_epdch_eta=0;

  //canv->cd(5);
  canv->cd(5)->SetLogy();
  TH1* h1pi0en = mH1F_NoEpdCutEn->DrawCopy("hist e");
  h1pi0en->SetLineColor(kBlack);
  TH1* h1phen = mH1F_EpdPhEn->DrawCopy("hist e same");
  h1phen->SetLineColor(kBlue);
  TH1* h1singlephen = mH1F_EpdSinglePhEn->DrawCopy("hist e same");
  h1singlephen->SetLineColor(kRed);
  TH1* h1singlechen = mH1F_EpdSingleChEn->DrawCopy("hist e same");
  h1singlechen->SetLineColor(kGreen+2);
  
  canv->cd(6)->SetLogy();
  TH1* h1pi0pt = mH1F_NoEpdCutPt->DrawCopy("hist e");
  h1pi0pt->SetLineColor(kBlack);
  TH1* h1phpt = mH1F_EpdPhPt->DrawCopy("hist e same");
  h1phpt->SetLineColor(kBlue);
  TH1* h1singlephpt = mH1F_EpdSinglePhPt->DrawCopy("hist e same");
  h1singlephpt->SetLineColor(kRed);
  TH1* h1singlechpt = mH1F_EpdSingleChPt->DrawCopy("hist e same");
  h1singlechpt->SetLineColor(kGreen+2);

  // canv->cd(7);
  // //canv->cd(7)->SetLogy();
  // mH1F_AllPi0Mass->SetLineColor(kBlack);
  // mH1F_AllPi0Mass->Draw("hist e");
  // mH1F_EpdPhAllMass->SetLineColor(kBlue);
  // mH1F_EpdPhAllMass->Draw("hist e same");
  // mH1F_EpdChAllMass->SetLineColor(kGreen+2);
  // mH1F_EpdChAllMass->Draw("hist e same");
  
  // canv->cd(8);
  // //canv->cd(8)->SetLogy();
  // TH1* h1pmass = mH1F_AllPi0Mass->DrawNormalized("hist e");
  // h1pmass->GetXaxis()->SetRangeUser(0,0.3);
  // h1pmass->SetLineColor(kBlack);
  // //h1pmass->Rebin(2);
  // TH1* h1phmass = mH1F_EpdPhAllMass->DrawNormalized("hist e same");
  // h1phmass->SetLineColor(kBlue);
  // //h1phmass->Rebin(2);
  // TH1* h1chmass = mH1F_EpdChAllMass->DrawNormalized("hist e same");
  // h1chmass->SetLineColor(kGreen+2);

  //canv->cd(9);
  //mH1F_AllPi0Mass->Draw("hist e");
  //mH1F_AllPi0Mass->GetXaxis()->SetRangeUser(0,0.3);
  //mH1F_AllPi0Mass->SetLineColor(kBlack);
  //mH1F_EpdPhAllMass->Draw("hist e same");
  //mH1F_EpdPhAllMass->SetLineColor(kBlue);
  //mH1F_EpdChAllMass->Draw("hist e same");
  //mH1F_EpdChAllMass->SetLineColor(kGreen+2);

  //canv->cd(10);//->SetLogy();
  //TH1* h1allmass = mH1F_AllPi0Mass->DrawNormalized("hist e");
  //  h1allmass->GetXaxis()->SetRangeUser(0,0.3);
  //h1allmass->SetLineColor(kBlack);
  //  h1allmass->Rebin(2);
  //TH1* h1phallmass = mH1F_EpdPhAllMass->DrawNormalized("hist e same");
  //h1phallmass->SetLineColor(kBlue);
  //  h1phallmass->Rebin(2);
  //TH1* h1challmass = mH1F_EpdChAllMass->DrawNormalized("hist e same");
  //h1challmass->SetLineColor(kGreen+2);

  canv->cd(7);
  TLegend* legendpad7 = new TLegend(0.5,0.5,0.93,0.93,"","nbNDC");
  mH1F_NoEpdCutAllMass->SetTitle("Invariant Mass distributions after most cuts");
  mH1F_NoEpdCutAllMass->Draw("hist e");
  mH1F_NoEpdCutAllMass->SetStats(0);
  mH1F_NoEpdCutAllMass->SetLineColor(kBlack);
  legendpad7->AddEntry(mH1F_NoEpdCutAllMass,"No Epd nmip cut","fle");
  //mH1F_InvMassAllCuts->UncheckedAt(0)->Draw("hist e same");
  //((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->SetStats(0);
  //((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->SetLineColor(kBlue);
  mH1F_EpdPhAllMass->Draw("hist e same");
  mH1F_EpdPhAllMass->SetStats(0);
  mH1F_EpdPhAllMass->SetLineColor(kBlue);
  ss_legname.str("");
  ss_legname << "Epd nmip<"<<mEpdNmipCut;
  legendpad7->AddEntry(mH1F_EpdPhAllMass,ss_legname.str().c_str(),"fle");
  mH1F_EpdSinglePhAllMass->Draw("hist e same");
  mH1F_EpdSinglePhAllMass->SetStats(0);
  mH1F_EpdSinglePhAllMass->SetLineColor(kRed);
  ss_legname.str("");
  ss_legname << "Epd nmip<"<<mEpdNmipCut << " single point";
  legendpad7->AddEntry(mH1F_EpdSinglePhAllMass,ss_legname.str().c_str(),"fle");
  mH1F_EpdSingleChAllMass->Draw("hist e same");
  mH1F_EpdSingleChAllMass->SetStats(0);
  mH1F_EpdSingleChAllMass->SetLineColor(kGreen+2);
  ss_legname.str("");
  ss_legname << "Epd nmip>="<<mEpdNmipCut << "  single point";
  legendpad7->AddEntry(mH1F_EpdSingleChAllMass,ss_legname.str().c_str(),"fle");  
  legendpad7->Draw();

  // canv->cd(12);
  // TLegend* legpad12 = new TLegend(0.5,0.5,0.93,0.93,"","nbNDC");
  // TH1* h1allcutmass = ((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->DrawNormalized("hist e"); //Draw this first as it has the largest y-value
  // h1allcutmass->SetTitle("Normalized Invariant Mass distributions different cuts");
  // h1allcutmass->SetStats(0);
  // h1allcutmass->SetLineColor(kBlue);
  // TH1* h1allcutbutepd = mH1F_NoEpdCutAllMass->DrawNormalized("hist e same");
  // h1allcutbutepd->SetStats(0);
  // h1allcutbutepd->SetLineColor(kBlack);
  // TH1* h1allcutepdch = mH1F_EpdChAllMass->DrawNormalized("hist e same");
  // h1allcutepdch->SetStats(0);
  // h1allcutepdch->SetLineColor(kGreen+2);
  // legpad12->AddEntry(h1allcutmass,"EPD nmip<0.7","fle");
  // legpad12->AddEntry(h1allcutbutepd,"No EPD nmip cut","fle");
  // legpad12->AddEntry(h1allcutepdch,"EPD nmip>=0.7","fle");
  // legpad12->Draw();
  
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdQa(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,3);

  //canv->cd(1);
  canv->cd(1)->SetLogy();
  ((TH1*)mH1F_Pi0MultAllCuts->UncheckedAt(0))->SetLineColor(kBlack);
  mH1F_EpdSinglePhPi0Mult->SetLineColor(kBlue);
  mH1F_EpdSingleChPi0Mult->SetLineColor(kGreen+2);
  mH1F_Pi0MultAllCuts->UncheckedAt(0)->Draw("hist e");
  mH1F_EpdSinglePhPi0Mult->Draw("hist e same");
  mH1F_EpdSingleChPi0Mult->Draw("hist e same");
  
  canv->cd(2);//->SetLogy();
  ((TH1*)mH1F_AllCuts_Zgg->UncheckedAt(0))->SetLineColor(kBlack);
  mH1F_EpdSinglePhZgg->SetLineColor(kBlue);
  mH1F_EpdSingleChZgg->SetLineColor(kGreen+2);
  mH1F_AllCuts_Zgg->UncheckedAt(0)->Draw("hist e");
  mH1F_EpdSinglePhZgg->Draw("hist e same");
  mH1F_EpdSingleChZgg->Draw("hist e same");
  
  //canv->cd(3);
  canv->cd(3)->SetLogy();
  ((TH1*)mH1F_AllCuts_Pi0En->UncheckedAt(0))->SetLineColor(kBlack);
  mH1F_EpdSinglePhEn->SetLineColor(kBlue);
  mH1F_EpdSingleChEn->SetLineColor(kGreen+2);
  mH1F_AllCuts_Pi0En->UncheckedAt(0)->Draw("hist e");
  mH1F_EpdSinglePhEn->Draw("hist e same");
  mH1F_EpdSingleChEn->Draw("hist e same");

  canv->cd(4);
  TH1* h1_pi0cut_phi = ((TH2*)mH2F_AllCuts_Pi0_etaVphi->UncheckedAt(0))->ProjectionX("h1_picut_phi");
  TH1* h1_singleph_phi = ((TH2*)mH2F_EpdSinglePh_etaVphi)->ProjectionX("h1_singleph_phi");
  TH1* h1_singlech_phi = ((TH2*)mH2F_EpdSinglePh_etaVphi)->ProjectionX("h1_singlech_phi");
  h1_pi0cut_phi->SetLineColor(kBlack);
  h1_singleph_phi->SetLineColor(kBlue);
  h1_singlech_phi->SetLineColor(kGreen+2);
  h1_pi0cut_phi->Draw("hist e");
  h1_singleph_phi->Draw("hist e same");
  h1_singlech_phi->Draw("hist e same");  

  canv->cd(5);
  TH1* h1_pi0cut_eta = ((TH2*)mH2F_AllCuts_Pi0_etaVphi->UncheckedAt(0))->ProjectionY("h1_picut_eta");
  TH1* h1_singleph_eta = ((TH2*)mH2F_EpdSinglePh_etaVphi)->ProjectionY("h1_singleph_eta");
  TH1* h1_singlech_eta = ((TH2*)mH2F_EpdSinglePh_etaVphi)->ProjectionY("h1_singlech_eta");
  h1_pi0cut_eta->SetLineColor(kBlack);
  h1_singleph_eta->SetLineColor(kBlue);
  h1_singlech_eta->SetLineColor(kGreen+2);
  h1_pi0cut_eta->Draw("hist e");
  h1_singleph_eta->Draw("hist e same");
  h1_singlech_eta->Draw("hist e same");
  
  //canv->cd(6);
  canv->cd(6)->SetLogy();
  TH1* h1_pi0cut_pt = ((TH2*)mH2F_AllCuts_Pi0_ptVeta->UncheckedAt(0))->ProjectionY("h1_pi0cut_pt");
  h1_pi0cut_pt->SetLineColor(kBlack);
  mH1F_EpdSinglePhPt->SetLineColor(kBlue);
  mH1F_EpdSingleChPt->SetLineColor(kGreen+2);
  h1_pi0cut_pt->Draw("hist e");
  mH1F_EpdSinglePhPt->Draw("hist e same");
  mH1F_EpdSingleChPt->Draw("hist e same");

  canv->cd(7);
  ((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->SetLineColor(kBlack);
  mH1F_EpdSinglePhAllMass->SetLineColor(kBlue);
  mH1F_EpdSingleChAllMass->SetLineColor(kGreen+2);  
  mH1F_InvMassAllCuts->UncheckedAt(0)->Draw("hist e");
  mH1F_EpdSinglePhAllMass->Draw("hist e same");
  mH1F_EpdSingleChAllMass->Draw("hist e same");
  
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintInvMassEpdQa(TCanvas* canv, const char* savename ) const
{
  canv->Clear();

  //canv->Divide(2,1);
  //canv->cd(1)->SetLogy();
  canv->cd();
  TLegend* legpad1 = new TLegend(0.7,0.7,0.93,0.93,"","nbNDC");
  /*
  TH1* h1_allmass = mH1F_AllPi0Mass->DrawCopy("hist e");
  h1_allmass->SetLineColor(kBlack);
  h1_allmass->SetTitle("Invariant Mass distributions with different cuts");
  h1_allmass->SetStats(0);
  //h1allmass->GetYaxis()->SetRangeUser(0,0.022);
  TH1* h1_allbutepdcut = mH1F_InvMassAllButEpdCut->DrawCopy("hist e same");
  h1_allbutepdcut->SetStats(0);
  h1_allbutepdcut->SetLineColor(kViolet);
  TH1* h1_allcutepdph = mH1F_EpdPhAllMass->DrawCopy("hist e same");
  h1_allcutepdph->SetStats(0);
  h1_allcutepdph->SetLineColor(kBlue);
  TH1* h1_allcutepdch = mH1F_EpdChAllMass->DrawCopy("hist e same");
  h1_allcutepdch->SetStats(0);
  h1_allcutepdch->SetLineColor(kGreen+2);
  TH1* h1_singleph = mH1F_EpdSinglePhAllMass->DrawCopy("hist e same");
  h1_singleph->SetStats(0);
  h1_singleph->SetLineColor(kRed);
  TH1* h1_singlech = mH1F_EpdSingleChAllMass->DrawCopy("hist e same");
  h1_singlech->SetLineColor(kOrange);
  h1_singlech->SetStats(0);
  */
  //h1allmass->GetYaxis()->SetRangeUser(0,0.022);
  TH1* h1_allbutepdcut = (TH1*)mH1F_NoEpdCutAllMass->Clone("h1_allbutepdcut");
  //h1_allbutepdcut->SetStats(0);
  //h1_allbutepdcut->SetLineColor(kViolet);
  TH1* h1_allcutepdph = (TH1*)mH1F_EpdPhAllMass->Clone("h1_allcutepdph");
  h1_allcutepdph->Divide(h1_allbutepdcut);
  h1_allcutepdph->SetStats(0);
  h1_allcutepdph->SetLineColor(kBlue);
  TH1* h1_allcutepdch = (TH1*)mH1F_EpdChAllMass->Clone("h1_allcutepdch");
  h1_allcutepdch->Divide(h1_allbutepdcut);
  h1_allcutepdch->SetStats(0);
  h1_allcutepdch->SetLineColor(kGreen+2);
  TH1* h1_singleph = (TH1*)mH1F_EpdSinglePhAllMass->Clone("h1_singleph");
  h1_singleph->Divide(h1_allbutepdcut);
  h1_singleph->SetStats(0);
  h1_singleph->SetLineColor(kRed);
  TH1* h1_singlech = (TH1*)mH1F_EpdSingleChAllMass->Clone("h1_singlech");
  h1_singlech->Divide(h1_allbutepdcut);
  h1_singlech->SetLineColor(kOrange);
  h1_singlech->SetStats(0);

  h1_allcutepdph->Draw("hist e");
  h1_allcutepdph->GetYaxis()->SetRangeUser(0,1);
  h1_allcutepdch->Draw("hist e same");
  h1_singleph->Draw("hist e same");
  h1_singlech->Draw("hist e same");
  //legpad1->AddEntry(h1_allmass,"All point pair","fle");
  //legpad1->AddEntry(h1_allbutepdcut,"All but EPD mip cut","fle");
  std::stringstream ss_legname;
  ss_legname << "EpdNmip<"<<mEpdNmipCut << "/CutMass";
  legpad1->AddEntry(h1_allcutepdph,ss_legname.str().c_str(),"fle");
  ss_legname.str("");
  ss_legname << "EpdNmip>="<<mEpdNmipCut << "/CutMass";
  legpad1->AddEntry(h1_allcutepdch,ss_legname.str().c_str(),"fle");
  ss_legname.str("");
  ss_legname << "EpdNmip<"<<mEpdNmipCut << " one point/CutMass";
  legpad1->AddEntry(h1_singleph,ss_legname.str().c_str(),"fle");
  ss_legname.str("");
  ss_legname << "EpdNmip>="<<mEpdNmipCut << " one point/CutMass";
  legpad1->AddEntry(h1_singlech,ss_legname.str().c_str(),"fle");
  legpad1->Draw();
  /*
  canv->cd(2);
  TLegend* legpad2 = new TLegend(0.5,0.5,0.93,0.93,"","nbNDC");
  TH1* h1_allmass_norm = mH1F_AllPi0Mass->DrawNormalized("hist e");
  h1_allmass_norm->SetLineColor(kBlack);
  h1_allmass_norm->SetTitle("Normalized Invariant Mass distributions with different cuts");
  h1_allmass_norm->SetStats(0);
  h1_allmass_norm->GetXaxis()->SetRangeUser(0,0.6);
  h1_allmass_norm->GetYaxis()->SetRangeUser(0,0.015);
  TH1* h1_allbutepdcut_norm = mH1F_NoEpdCutAllMass->DrawNormalized("hist e same");
  h1_allbutepdcut_norm->SetStats(0);
  h1_allbutepdcut_norm->SetLineColor(kViolet);
  TH1* h1_allcutepdph_norm = mH1F_EpdPhAllMass->DrawNormalized("hist e same");
  h1_allcutepdph_norm->SetStats(0);
  h1_allcutepdph_norm->SetLineColor(kBlue);
  TH1* h1_allcutepdch_norm = mH1F_EpdChAllMass->DrawNormalized("hist e same");
  h1_allcutepdch_norm->SetStats(0);
  h1_allcutepdch_norm->SetLineColor(kGreen+2);
  TH1* h1_singleph_norm = mH1F_EpdSinglePhAllMass->DrawNormalized("hist e same");
  h1_singleph_norm->SetStats(0);
  h1_singleph_norm->SetLineColor(kRed);
  TH1* h1_singlech_norm = mH1F_EpdSingleChAllMass->DrawNormalized("hist e same");
  h1_singlech_norm->SetLineColor(kOrange);
  h1_singlech_norm->SetStats(0);

  //legpad2->AddEntry(h1_allmass_norm,"All point pair","fle");
  legpad2->AddEntry(h1_allbutepdcut_norm,"All but EPD nmip cut","fle");
  legpad2->AddEntry(h1_allcutepdph_norm,"EPD nmip<0.7","fle");
  legpad2->AddEntry(h1_allcutepdch_norm,"EPD nmip>=0.7","fle");
  legpad2->AddEntry(h1_singleph_norm,"EPD nmip<0.7 one point","fle");
  legpad2->AddEntry(h1_singlech_norm,"EPD nmip>=0.7 one point","fle");
  legpad2->Draw();
  */
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::AddHistStatsOneline( TLegend* HistLeg, const TH1* h1, const std::string &title )
{
  //This function is good for when many histograms are plotted
  if( HistLeg==0 || h1==0 ){return;}
  if( h1->GetDimension()==1 ){
    std::stringstream ss_entry;
    if( title.size()==0 ){ ss_entry << h1->GetName(); }
    else{ ss_entry << title; }
    ss_entry << "|E:"<< h1->GetEntries();
    ss_entry << "|M:"<< h1->GetMean();
    ss_entry << "|R:"<< h1->GetRMS();
    //ss_entry << "|U:"<< h1->GetBinContent(0);
    //ss_entry << "|O:"<< h1->GetBinContent(h1->GetNbinsX()+1);
    HistLeg->AddEntry(h1, ss_entry.str().c_str(),"fle" );
  }
}

void StMuFcsPi0TreeMaker::PaintEnergyZoom(TCanvas* canv, const char* savename) const
{
  canv->Clear();

  canv->Divide(2,2);

  canv->cd(1);
  TH1* h1_encopy = (TH1*)mH1F_ClusterEnergy->DrawClone("hist e");
  h1_encopy->GetXaxis()->SetRangeUser(75,85);
  h1_encopy->SetLineColor(kBlack);

  canv->cd(2)->SetLogz();
  mH2F_PhotonHeatMapG->SetStats(0);
  mH2F_PhotonHeatMapG->Draw("colz");

  canv->cd(3)->SetLogz();
  mH2F_PhotonHeatMapB->SetStats(0);
  mH2F_PhotonHeatMapB->Draw("colz");
  
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdNmipCuts(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(2,1);
  for( UInt_t i=0; i<2; ++i ){
    canv->cd(i+1);
    for( Int_t icut=0, ipad=1; icut<NEPDCUTS; ++icut,++ipad ){
      ((TH1*)mH1F_InvMassEpdCuts[i]->UncheckedAt(icut))->SetLineColor(icut+1);//Hack to get rainbow colors since I know I only have 8 histograms
      if( icut==0 ){ ((TH1*)mH1F_InvMassEpdCuts[i]->UncheckedAt(icut))->Draw("hist e"); }
      else{  ((TH1*)mH1F_InvMassEpdCuts[i]->UncheckedAt(icut))->Draw("hist e same"); }
    }
  }
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintPi0Cuts(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(3,3);
  canv->cd(1)->SetLogy();
  //mH1F_NFoundPhiBin->Draw("hist e");
  mH1F_AllCuts_xF->UncheckedAt(0)->Draw("hist e");
  canv->cd(2)->SetLogy();
  mH1F_AllCuts_xFZoom->UncheckedAt(0)->Draw("hist e");  
  //canv->cd(3)->SetLogy();
  //mH1F_AllCuts_Pi0En->Draw("hist e");
  canv->cd(3);
  mH2F_AllCuts_Pi0_massVen->UncheckedAt(0)->Draw("colz");
  canv->cd(4);
  //mH1F_AllCuts_Pi0Phi->Draw("hist e");
  mH2F_AllCuts_Pi0_etaVphi->UncheckedAt(0)->Draw("colz");
  canv->cd(5);
  mH1F_InvMassAllCuts->UncheckedAt(0)->Draw("hist e");
  canv->cd(6)->SetLogy();
  mH1F_NBadEpdProj->Draw("hist e");
  mH1F_NBadEpdProjVcut->Draw("hist e same");
  canv->cd(7);
  mH2F_AllCuts_Pi0_yVx->UncheckedAt(0)->Draw("colz");
  canv->cd(8);
  mH1F_Pi0MultAllCuts->UncheckedAt(0)->Draw("hist e");
  canv->cd(9);
  TH1* h1f_multcut_clone = (TH1*)mH1F_Pi0MultAllCuts->UncheckedAt(0)->Clone("h1f_multcut_clone");
  h1f_multcut_clone->SetBinContent(1,0); //Artifically delete all entries in zero bin to zoom in on the higher regions
  h1f_multcut_clone->Draw("hist e");
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintInvMassCuts(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  for( short ixbin=0; ixbin<NXFBIN; ++ixbin ){
    canv->DivideSquare(NPHIBIN);
    for( short phibin=0; phibin<NPHIBIN; ++phibin ){
      canv->cd(phibin+1);
      std::stringstream histname;
      histname << "H1F_InvMass_xf"<<ixbin << "_phi"<<phibin;
      TH1D* hist_proj = ((TH3*)mH3F_AllCutsInvMass_xfVphi)->ProjectionZ( histname.str().c_str(), phibin+1,phibin+1, ixbin+1,ixbin+1 );
      hist_proj->SetTitle( histname.str().c_str() );
      hist_proj->Draw("hist e");
      //mH1F_InvMassAllCutsByEnByPhi[ebin][phibin]->Draw("hist e");
    }
    canv->Print(savename);
    canv->Clear();  //Should hopefully delete the projection histograms 
  }
}

void StMuFcsPi0TreeMaker::PaintNpi0Inc(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(2,2);
  int ipad = 1;
  for( int ibeam=0; ibeam<2; ++ibeam ){
    for( int ispin=0; ispin<2; ++ispin ){
      canv->cd(ipad++);
      mH2F_NPi0Inc_xfVphi[ibeam][ispin]->SetStats(0);
      mH2F_NPi0Inc_xfVphi[ibeam][ispin]->Draw("colz");
    }
  }
  canv->Print(savename);
  /*
  for( short ebin=0; ebin<NENERGYBIN; ++ebin ){
    canv->DivideSquare(NPHIBIN);
    for( short phibin=0; phibin<NPHIBIN; ++phibin ){
      canv->cd(phibin+1);
      //mH1F_NPi0ByEnByPhi[ebin][phibin]->Draw();
    }
    canv->Print(savename);
    canv->Clear();
  }
  */
}

void StMuFcsPi0TreeMaker::PaintNpi0Bg1(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(2,2);
  int ipad = 1;
  for( int ibeam=0; ibeam<2; ++ibeam ){
    for( int ispin=0; ispin<2; ++ispin ){
      canv->cd(ipad++);
      mH2F_NPi0Bg1_xfVphi[ibeam][ispin]->SetStats(0);
      mH2F_NPi0Bg1_xfVphi[ibeam][ispin]->Draw("colz");
    }
  }
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintNpi0Bg2(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(2,2);
  int ipad = 1;
  for( int ibeam=0; ibeam<2; ++ibeam ){
    for( int ispin=0; ispin<2; ++ispin ){
      canv->cd(ipad++);
      mH2F_NPi0Bg2_xfVphi[ibeam][ispin]->SetStats(0);
      mH2F_NPi0Bg2_xfVphi[ibeam][ispin]->Draw("colz");
    }
  }
  canv->Print(savename);
}

  
Int_t StMuFcsPi0TreeMaker::LoadGraphsFromFile(TFile* file, TObjArray* graphs )
{
  Int_t gloaded = 0;
  gloaded += StMuFcsRun22QaMaker::MakeGraph(file,graphs,mGE_AllCuts_InvMass,"GE_AllCuts_InvMass","Mean of Invariant Mass (Err=RMS) vs. Run Index");
  gloaded += StMuFcsRun22QaMaker::MakeGraph(file,graphs,mGE_AllCuts_Pi0En,"GE_AllCuts_Pi0En","Mean Pi0 Energy (Err=RMS) vs. Run Index");
  return gloaded;
}


void StMuFcsPi0TreeMaker::FillGraphs(Int_t irun)
{
  //TF1* pi0gausfit = new TF1("pi0gausfit","gaus(0)",0.1,0.2);
  //mH1F_InvMassAllCuts->Fit(pi0gausfit,"RQN");
  mGE_AllCuts_InvMass->SetPoint(irun,irun,((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->GetMean());
  mGE_AllCuts_InvMass->SetPointError(irun,0,((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->GetRMS());
  //delete pi0gausfit;
  //pi0gausfit = 0;

  mGE_AllCuts_Pi0En->SetPoint(irun,irun,((TH1*)mH1F_AllCuts_Pi0En->UncheckedAt(0))->GetMean());
  mGE_AllCuts_Pi0En->SetPointError(irun,0,((TH1*)mH1F_AllCuts_Pi0En->UncheckedAt(0))->GetRMS());

  return;
}

void StMuFcsPi0TreeMaker::DrawQaGraphs(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->cd();
  canv->Divide(1,2);

  canv->cd(1);
  mGE_AllCuts_InvMass->Draw("AL");
  canv->cd(2);
  mGE_AllCuts_Pi0En->Draw("AL");
  
  canv->SaveAs(savename);
}


void StMuFcsPi0TreeMaker::MergeForTssa( TH1* totalhistinc[][2], TH1* totalhistbg1[][2], TH1* totalhistbg2[][2], TH3* mergedinvmass, TH1* mergedpolblue, TH1* mergedpolyell, TH1* mergedpolblueerr, TH1* mergedpolyellerr )
{
  if( mH1F_RndmSpin->GetBinContent(1)>0.1 ){ std::cout << "  + RandomSpinFound" << std::endl; return; } //Don't merge histograms from files with random spin patterns
  for( int ibeam=0; ibeam<2; ++ibeam ){
    for( int ispin=0; ispin<2; ++ispin ){
      totalhistinc[ibeam][ispin]->Add(mH2F_NPi0Inc_xfVphi[ibeam][ispin]);
      totalhistbg1[ibeam][ispin]->Add(mH2F_NPi0Bg1_xfVphi[ibeam][ispin]);
      totalhistbg2[ibeam][ispin]->Add(mH2F_NPi0Bg2_xfVphi[ibeam][ispin]);
    }
  }
  mergedinvmass->Add(mH3F_AllCutsInvMass_xfVphi);
  mergedpolblue->Add(mH1D_BluePol);
  mergedpolblueerr->Add(mH1D_BluePolErr);
  mergedpolyell->Add(mH1D_YellowPol);
  mergedpolyellerr->Add(mH1D_YellowPolErr);
  //mergedpoldata->Add(mH1D_Entries);
}

void StMuFcsPi0TreeMaker::DoTssaAna( TH1* npi0[][2], TH1* h1_rawasymVphi[][StMuFcsPi0TreeMaker::NXFBIN] )
{
  //Double_t pi = TMath::Pi();
  //TH1* h1_rawasymVphi[2][NENERGYBIN];  //Histograms of computed raw asymmetries for blue and yellow beams for each energy bin
  //memset(h1_rawasymVphi,0,sizeof(h1_rawasymVphi));
  
  //for( int ibeam = 0; ibeam<2; ++ibeam ){
  for( int ixbin = 0; ixbin<NXFBIN; ++ixbin ){
    int ibeam = 0; //Blue beam first since that uses the normal STAR direction convention
    for( int iphibin=0; iphibin<NPHIBIN/2 ; ++iphibin ){
      Double_t nupatphi         = npi0[ibeam][0]->GetBinContent(npi0[ibeam][0]->GetBin(iphibin+1,ixbin+1));
      Double_t ndownatphipluspi = npi0[ibeam][1]->GetBinContent(npi0[ibeam][1]->GetBin(iphibin+1+NPHIBIN/2,ixbin+1));
      Double_t ndownatphi       = npi0[ibeam][1]->GetBinContent(npi0[ibeam][1]->GetBin(iphibin+1,ixbin+1));
      Double_t nupatphipluspi   = npi0[ibeam][0]->GetBinContent(npi0[ibeam][0]->GetBin(iphibin+1+NPHIBIN/2,ixbin+1));
      
      Double_t nupdown = sqrt( nupatphi * ndownatphipluspi );
      Double_t ndownup = sqrt( ndownatphi * nupatphipluspi );
      Double_t numerator = nupdown - ndownup;
      Double_t denominator = nupdown + ndownup;
      //std::cout << "|iphibin:"<<iphibin << "|iebin:"<<iebin << "|ibeam:"<<0 << std::endl;
      //std::cout << " + |nu0:"<< nupatphi << "|ndp:"<<ndownatphipluspi << "|nd0:"<<ndownatphi << "|nup:"<<nupatphipluspi << "|numer:"<<numerator << "|denom:"<<denominator << std::endl;
	
      if( denominator != 0 ){
	// Propagate errors for A_N
	Double_t AnErr = sqrt(ndownatphi*nupatphi*nupatphipluspi + ndownatphipluspi*nupatphi*nupatphipluspi + ndownatphi*ndownatphipluspi*nupatphi + ndownatphi*ndownatphipluspi*nupatphipluspi)/(denominator*denominator);
	//std::cout << " + |An:"<<numerator/denominator<< "|nu0_err:"<<nupatphi_err << "|ndp_err:"<<ndownatphipluspi_err << "|nd0_err:"<<ndownatphi_err << "|nup_err:"<< nupatphipluspi_err << "|AnErr:"<<AnErr << std::endl;
	  
	h1_rawasymVphi[ibeam][ixbin]->SetBinContent(iphibin+1,numerator/denominator);
	h1_rawasymVphi[ibeam][ixbin]->SetBinError(iphibin+1,AnErr);
      }
      else{
	//If denomnator is zero set bin to zero to avoid bad division and set large error
	h1_rawasymVphi[ibeam][ixbin]->SetBinContent(iphibin+1,0);
	h1_rawasymVphi[ibeam][ixbin]->SetBinError(iphibin+1,1);
      }
    }
    ibeam = 1; //Yellow beam next since need to flip the order because "left" for yellow is the opposite of the "left" for blue
    for( int iphibin=NPHIBIN/2; iphibin<NPHIBIN; ++iphibin ){
      Double_t nupatphi         = npi0[ibeam][0]->GetBinContent(npi0[ibeam][0]->GetBin(iphibin+1,ixbin+1));
      Double_t ndownatphipluspi = npi0[ibeam][1]->GetBinContent(npi0[ibeam][1]->GetBin(iphibin+1-NPHIBIN/2,ixbin+1));
      Double_t ndownatphi       = npi0[ibeam][1]->GetBinContent(npi0[ibeam][1]->GetBin(iphibin+1,ixbin+1));
      Double_t nupatphipluspi   = npi0[ibeam][0]->GetBinContent(npi0[ibeam][0]->GetBin(iphibin+1-NPHIBIN/2,ixbin+1));
      
      Double_t nupdown = sqrt( nupatphi * ndownatphipluspi );
      Double_t ndownup = sqrt( ndownatphi * nupatphipluspi );
      Double_t numerator = nupdown - ndownup;
      Double_t denominator = nupdown + ndownup;
      //std::cout << "|iphibin:"<<iphibin << "|iebin:"<<iebin << "|ibeam:"<<1 << std::endl;
      //std::cout << " + |nu0:"<< nupatphi << "|ndp:"<<ndownatphipluspi << "|nd0:"<<ndownatphi << "|nup:"<<nupatphipluspi << "|numer:"<<numerator << "|denom:"<<denominator << std::endl;
	
      if( denominator != 0 ){
	// Propagate errors for A_N
	Double_t AnErr = sqrt(ndownatphi*nupatphi*nupatphipluspi + ndownatphipluspi*nupatphi*nupatphipluspi + ndownatphi*ndownatphipluspi*nupatphi + ndownatphi*ndownatphipluspi*nupatphipluspi)/(denominator*denominator);
	//std::cout << " + |An:"<<numerator/denominator << "|AnErr:"<<AnErr << std::endl;
	h1_rawasymVphi[ibeam][ixbin]->SetBinContent(iphibin+1-NPHIBIN/2,numerator/denominator); //instead of flipping can also multiply by -1
	h1_rawasymVphi[ibeam][ixbin]->SetBinError(iphibin+1-NPHIBIN/2,AnErr);
      }
      else{
	//If denomnator is zero set bin to zero to avoid bad division and set large error
	h1_rawasymVphi[ibeam][ixbin]->SetBinContent(iphibin+1-NPHIBIN/2,0);
	h1_rawasymVphi[ibeam][ixbin]->SetBinError(iphibin+1-NPHIBIN/2,1);
      }
    }
  }
}

void StMuFcsPi0TreeMaker::DoTssaFit( TH1* h1_rawasymVphi[][StMuFcsPi0TreeMaker::NXFBIN], TH1* h1_bluepoldata, TH1* h1_yellowpoldata, TH1* h1_AnResult[] )
{
  Double_t pi = TMath::Pi();
  //Double_t ebins[NENERGYBIN+1] = {0, 15, 20, 25, 30, 40, 55, 70, 100};  //Copied from above
  //h1_AnResult[0] = new TH1D("H1D_AnResultBlue","A_N vs. Energy;Energy;A_N",StMuFcsPi0TreeMaker::NENERGYBIN,ebins);
  //h1_AnResult[1] = new TH1D("H1D_AnResultYellow","A_N vs. Energy;Energy;A_N",StMuFcsPi0TreeMaker::NENERGYBIN,ebins);
  //Double_t nentries = h1_poldata->GetBinContent(1);
  for( int ibeam = 0; ibeam<2; ++ibeam ){
    for( int ixbin = 0; ixbin<NXFBIN; ++ixbin ){
      //std::stringstream ss_fname;
      //ss_fname << "AnFit_" << ibeam << "_" << iebin;
      //TF1* AnFit = new TF1("AnFit","[0]*cos(x+[1])+[2]",-pi/2.0,pi/2.0);
      TF1* AnFit = new TF1("AnFit","[0]*cos(x+[1])",-pi/2.0,pi/2.0);
      AnFit->SetParName(0,"A");
      AnFit->SetParName(1,"#phi_{0}");
      //AnFit->SetParName(2,"yoff");
      h1_rawasymVphi[ibeam][ixbin]->Fit(AnFit,"R");
      Double_t rawasym = AnFit->GetParameter(0);
      //Double_t totalpolpercent = h1_poldata->GetBinContent(ibeam+2);  //Since blue beam [0] was stored in bin 2 and yellow beam [1] was stored in bin 3
      //Double_t lumcorrectedpol = totalpolpercent/nentries;            //luminosity corrected polarization
      Double_t lumcorrectedpol = (ibeam==0)?(h1_bluepoldata->GetMean()/100.0):(h1_yellowpoldata->GetMean()/100.0); //Divide by 100 since value is a percent
      h1_AnResult[ibeam]->SetBinContent( ixbin+1, rawasym/lumcorrectedpol );
      h1_AnResult[ibeam]->SetBinError( ixbin+1, AnFit->GetParError(0)/lumcorrectedpol );
    }
  }
}

void StMuFcsPi0TreeMaker::DoPi0Fits(TH3* mH3F_invmass, TH1* hist_proj[] )
{
  for( short ebin=0; ebin<NXFBIN; ++ebin ){
    if( hist_proj[ebin]==0 ){
      std::stringstream histname;
      histname << "H1F_InvMass_xf"<<ebin;
      hist_proj[ebin] = (TH1*)((TH3F*)mH3F_invmass)->ProjectionZ( histname.str().c_str(), 1,NPHIBIN, ebin+1,ebin+1 );
      hist_proj[ebin]->SetTitle( histname.str().c_str() );
    }
    TF1* PeakFit = new TF1("PeakFit",skewgaus,0.1,0.2,4);
    PeakFit->SetParameter(0,hist_proj[ebin]->GetBinContent(hist_proj[ebin]->GetMaximumBin()));
    PeakFit->SetParameter(1,0.135);
    PeakFit->SetParameter(2,0.03);
    PeakFit->SetParameter(3,0);
    //PeakFit->SetParLimits(3,0,10); //Force positive skew
    //PeakFit->FixParameter(3,0); //debuging for gaussian
    //TF1* Bg1Fit = new TF1("Bg1Fit","[0] +[1]*x +[2]*x*x +[3]*x*x*x + [4]*x*x*x*x",0.3,0.4);
    //TF1* Bg1Fit = new TF1("Bg1Fit","1 +[0]*x +[1]*(2*x*x-1) +[2]*(4*x*x*x-3*x)",0.3,0.4); //Fit to Chebyshev
    //TF1* Bg1Fit = new TF1("Bg1Fit","[0] +[1]*x*x +[2]*x*x*x*x",0.3,0.4); //Fit to even
    //TF1* Bg1Fit = new TF1("Bg1Fit","[0] +[1]*x +[2]*x*x*x",0.3,0.4);     //Fit to odd
    TF1* Bg1Fit = new TF1("Bg1Fit",pol4bg,0,0.9,5);
    TF1* Bg2Fit = new TF1("Bg2Fit","[0] +[1]*x +[2]*x*x +[3]*x*x*x + [4]*x*x*x*x",0.7,0.9);
    hist_proj[ebin]->Fit(PeakFit,"QR");
    hist_proj[ebin]->Fit(Bg1Fit,"QR+");
    hist_proj[ebin]->Fit(Bg2Fit,"QR+");
    //TF1* GlobalFit = new TF1("GlobalFit","[0] +[1]*x +[2]*x*x +[3]*x*x*x +[4]*x*x*x*x + gaus(5)",0,0.45);
    //TF1* GlobalFit = new TF1("GlobalFit","1 +[0]*x +[1]*(2*x*x-1) +[2]*(4*x*x*x-3*x) + gaus(3)",0,0.4);
    //TF1* GlobalFit = new TF1("GlobalFit",pol4skewgaus,0,0.45,9);
    TF1* GlobalFit = new TF1("GlobalFit",pol4skewgaus,0.0,1.0,9);
    GlobalFit->SetParameter(0,Bg1Fit->GetParameter(0));
    GlobalFit->SetParameter(1,Bg1Fit->GetParameter(1));
    GlobalFit->SetParameter(2,Bg1Fit->GetParameter(2));
    GlobalFit->SetParameter(3,Bg1Fit->GetParameter(3));
    GlobalFit->SetParameter(4,Bg1Fit->GetParameter(4));
    //GlobalFit->SetParameter(5,PeakFit->GetParameter(0)); //DOn't use from first fit
    GlobalFit->SetParameter(5,hist_proj[ebin]->GetBinContent(hist_proj[ebin]->GetMaximumBin()));
    GlobalFit->SetParName(5,"amplitude");
    //GlobalFit->SetParameter(6,PeakFit->GetParameter(1));
    GlobalFit->SetParameter(6,0.135);
    GlobalFit->SetParName(6,"mean");
    GlobalFit->SetParameter(7,0.03);
    GlobalFit->SetParName(7,"sigma");
    GlobalFit->SetParameter(8,0);  //Fix the skew parameter in global fit since skew should only be affected by peak region
    GlobalFit->SetParName(8,"skew");
    //GlobalFit->SetParameter(3,PeakFit->GetParameter(0));
    //GlobalFit->SetParameter(4,PeakFit->GetParameter(1));
    //GlobalFit->SetParameter(5,PeakFit->GetParameter(2));
    hist_proj[ebin]->Fit(GlobalFit,"R+");
    TF1* BgFit = new TF1("BgFit","[0] +[1]*x +[2]*x*x +[3]*x*x*x +[4]*x*x*x*x", 0,1.0);
    //TF1* BgFit = new TF1("BgFit","1 +[0]*x +[1]*(2*x*x-1) +[2]*(4*x*x*x-3*x)",0,0.5); //Fit to Chebyshev
    BgFit->FixParameter(0,Bg1Fit->GetParameter(0));
    BgFit->FixParameter(1,Bg1Fit->GetParameter(1));
    BgFit->FixParameter(2,Bg1Fit->GetParameter(2));
    BgFit->FixParameter(3,Bg1Fit->GetParameter(3));
    BgFit->FixParameter(4,Bg1Fit->GetParameter(4));
    hist_proj[ebin]->Fit(BgFit,"RQ+");
    TF1* BgGlobalFit = new TF1("BgGlobalFit","[0] +[1]*x +[2]*x*x +[3]*x*x*x +[4]*x*x*x*x", 0,1.0);
    //TF1* BgFit = new TF1("BgGlobalFit","1 +[0]*x +[1]*(2*x*x-1) +[2]*(4*x*x*x-3*x)",0,0.5); //Fit to Chebyshev
    BgGlobalFit->FixParameter(0,GlobalFit->GetParameter(0));
    BgGlobalFit->FixParameter(1,GlobalFit->GetParameter(1));
    BgGlobalFit->FixParameter(2,GlobalFit->GetParameter(2));
    BgGlobalFit->FixParameter(3,GlobalFit->GetParameter(3));
    BgGlobalFit->FixParameter(4,GlobalFit->GetParameter(4));
    hist_proj[ebin]->Fit(BgGlobalFit,"RQ+");
  }
}

void StMuFcsPi0TreeMaker::DoBgCorrectedAn(TH1* h1_invmass_xf[], TH1* h1_an_inc, TH1* h1_an_bg, TH1* h1_anresult )
{
  if( h1_anresult==0 ){ return; }
  if( h1_invmass_xf==0 ){ return; }
  if( h1_an_inc==0 ){ return; }
  if( h1_an_bg==0 ){ return; }
  for( Int_t ixbin=0; ixbin<NXFBIN; ++ixbin ){
    TF1* bgfunc = 0;
    if( ixbin<=(NXFBIN-1) ){ bgfunc = h1_invmass_xf[ixbin]->GetFunction("BgGlobalFit"); }
    //if( ixbin<4 ){ bgfunc = h1_invmass_en[ixbin]->GetFunction("BgGlobalFit"); }
    else{ bgfunc = h1_invmass_xf[ixbin]->GetFunction("Bg2Fit"); }  //For all but the last energy bin this function had the most reasonable background shape(*/
    if( bgfunc==0 ){ continue; }
    Double_t xlowbin = h1_invmass_xf[ixbin]->FindBin(0.1);
    Double_t xhighbin = h1_invmass_xf[ixbin]->FindBin(0.2);
    Double_t npi0_inc = h1_invmass_xf[ixbin]->Integral(xlowbin,xhighbin,"width");
    Double_t npi0_bg = bgfunc->Integral(0.1,0.2);
    Double_t ratio = npi0_bg/npi0_inc;
    Double_t ansignal = (h1_an_inc->GetBinContent(ixbin+1) - ratio*h1_an_bg->GetBinContent(ixbin+1)) / (1.0-ratio);
    Double_t aninc_err = h1_an_inc->GetBinError(ixbin+1);
    Double_t anbg_err = h1_an_bg->GetBinError(ixbin+1);
    Double_t ansig_err = (1.0/(1.0-ratio))*sqrt(aninc_err*aninc_err + ratio*ratio*anbg_err*anbg_err);
    h1_anresult->SetBinContent(ixbin+1, ansignal);
    h1_anresult->SetBinError(ixbin+1, ansig_err );
  }
}

void StMuFcsPi0TreeMaker::PolData::Print() const
{
  std::cout << " * |fill:"<<mFillNum
	    << "|En:"<<mBeamEn
	    << "|StartTime:"<<mStartTime
	    << "|BP0:"<<mBlueP0
	    << "|BErrP0:"<<mBlueErrP0
	    << "|BdPdT:"<<mBluedPdT
	    << "|BErrdPdT:"<<mBlueErrdPdT
    	    << "|YP0:"<<mYellowP0
	    << "|YErrP0:"<<mYellowErrP0
	    << "|YdPdT:"<<mYellowdPdT
	    << "|YErrdPdT:"<<mYellowErrdPdT
	    << std::endl;
}

Int_t StMuFcsPi0TreeMaker::ReadPolFile(const char* filename)
{
  std::ifstream in_polfile(filename);       //input (in_) polarization (pol) file
  if( !in_polfile.is_open() ){ return 0; }
  Int_t fillnum = 0;
  Int_t value = 0;
  Double_t fvalue = 0;
  //[Don't use eof in while loop since bit only gets set after reading file so it will loop twice on last line](https://stackoverflow.com/questions/5605125/why-is-iostreameof-inside-a-loop-condition-i-e-while-stream-eof-cons)
  while( in_polfile >> fillnum ){
    if( fillnum<=32000 ){ break; }
    PolData* temp = new PolData();
    temp->mFillNum = fillnum;                             //First entry is fill number
    in_polfile >> value; temp->mBeamEn = value;           //Second entry is beam energy
    in_polfile >> value; temp->mStartTime = value;        //Third entry is start time
    in_polfile >> value;                                  //Fourth entry is stop time (drop)
    in_polfile >> fvalue; temp->mBlueP0 = fvalue;         //Fifth entry is blue beam polarization
    in_polfile >> fvalue; temp->mBlueErrP0 = fvalue;      //Sixth entry is blue beam polarization error
    in_polfile >> fvalue; temp->mBluedPdT = fvalue;       //Seventh entry is blue beam polarization decay
    in_polfile >> fvalue; temp->mBlueErrdPdT = fvalue;    //Eighth entry is the blue beam polarization decay error
    in_polfile >> fvalue; temp->mYellowP0 = fvalue;       //Ninth entry is yellow beam polarization
    in_polfile >> fvalue; temp->mYellowErrP0 = fvalue;    //Tenth entry is yellow beam polarization error
    in_polfile >> fvalue; temp->mYellowdPdT = fvalue;     //Eleventh entry is yellow beam polarization decay
    in_polfile >> fvalue; temp->mYellowErrdPdT = fvalue;  //Twelfth entry is the yellow beam polarization decay error
    mPolarizationData[fillnum] = temp;
  }
  in_polfile.close();
  //for( std::map<Int_t,PolData*>::iterator itr=mPolarizationData.begin(); itr!=mPolarizationData.end(); ++itr ){ itr->second->Print(); }
  return mPolarizationData.size();
}

Double_t StMuFcsPi0TreeMaker::pol4bg(Double_t* x, Double_t* par)
{
  //Rejecting these points gives good background at both the high and low end
  if( x[0]<0.045 ){ TF1::RejectPoint(); return 0; }
  if( 0.09<x[0] && x[0]<0.21 ){ TF1::RejectPoint(); return 0; }
  if( 0.4<x[0] && x[0]<0.6 ){ TF1::RejectPoint(); return 0; }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t StMuFcsPi0TreeMaker::skewgaus(Double_t*x, Double_t* par)
{
  if( par[2]==0 ){ return 1.e30; }
  Double_t xarg = (x[0]-par[1])/par[2];
  Double_t gaus = (1.0/(TMath::Sqrt(2.0*TMath::Pi())))*TMath::Exp(-0.5*xarg*xarg);
  Double_t skew = 0.5*(1.0+TMath::Erf( (par[3]*xarg)/TMath::Sqrt2() ));
  return (2.0/par[2])*par[0]*gaus*skew;
}

Double_t StMuFcsPi0TreeMaker::pol4skewgaus(Double_t*x, Double_t* par)
{
  //if( x[0]<0.045 ){ TF1::RejectPoint(); return 0; }
  if( x[0]<0.045 || (0.4<x[0] && x[0]<0.6) ){ TF1::RejectPoint(); return 0; } //Avoid low mass region and 0.4 to 0.6 to skip eta meson mass
  Double_t pol4 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
  if( par[7]==0 ){ return 1.e30; }
  Double_t xarg = (x[0]-par[6])/par[7];
  Double_t gaus = (1.0/(TMath::Sqrt(2.0*TMath::Pi())))*TMath::Exp(-0.5*xarg*xarg);
  Double_t skew = 0.5*(1.0+TMath::Erf( (par[8]*xarg)/TMath::Sqrt2() ));
  return pol4 + (2.0/par[7])*par[5]*gaus*skew;
}

void StMuFcsPi0TreeMaker::PaintPhotonQaForDefense(TCanvas* canv, const char* savename)   const
{
  canv->Clear();
  
  canv->Divide(2,2);
  canv->cd(1)->SetLogz();
  mH2F_PhotonHeatMap->SetStats(0);
  mH2F_PhotonHeatMap->Draw("colz");
  canv->cd(2)->SetLogz();
  mH2F_EpdProjHitMap_Vcut->SetStats(0);
  mH2F_EpdProjHitMap_Vcut->Draw("colz");
  canv->cd(3)->SetLogy();
  TH1* h1_nmippoint = ((TH2*)mH2F_EpdNmip)->ProjectionY("h1_nmippoint",2,2);
  h1_nmippoint->SetTitle("EPD nMIP for all found points;nMIP;");
  h1_nmippoint->Draw("hist e");
  canv->cd(4)->SetLogy();
  mH1F_PointEnergy->Draw("hist e");

  canv->Print(savename);
  mH2F_PhotonHeatMap->SetStats(1);
  mH2F_EpdProjHitMap_Vcut->SetStats(1);
}


void StMuFcsPi0TreeMaker::PaintPi0QaForDefense(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  canv->Divide(3,2);

  canv->cd(1)->SetLogy();
  mH1F_AllCuts_Pi0En->UncheckedAt(0)->Draw("hist e");

  canv->cd(2);
  mH2F_AllCuts_Pi0_etaVphi->UncheckedAt(0)->Draw("colz");

  canv->cd(3);
  mH1F_InvMassAllCuts->UncheckedAt(0)->Draw("hist e");

  canv->cd(4)->SetLogz();
  ((TH1*)mH2F_AllCuts_Pi0_yVx->UncheckedAt(0))->SetStats(0);
  mH2F_AllCuts_Pi0_yVx->UncheckedAt(0)->Draw("colz");

  canv->cd(5);//->SetLogy();
  TH1* h1phzgg = (mH1F_EpdPhZgg->GetEntries()>0) ? mH1F_EpdPhZgg->DrawNormalized("hist e") : mH1F_EpdPhZgg->DrawCopy("hist e");
  h1phzgg->SetLineColor(kBlue);
  
  canv->cd(6);
  TLegend* legpad12 = new TLegend(0.5,0.5,0.93,0.93,"","nbNDC");
  TH1* h1allcutmass = ((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->GetEntries()>0 ? ((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->DrawNormalized("hist e") : ((TH1*)mH1F_InvMassAllCuts->UncheckedAt(0))->DrawCopy("hist e"); //Draw this first as it has the largest y-value
  h1allcutmass->SetTitle("Normalized Invariant Mass distributions after most cuts;M_{inv} (GeV/c^{2});");
  h1allcutmass->SetStats(0);
  h1allcutmass->SetLineColor(kBlue);
  TH1* h1allcutbutepd = ((TH1*)mH1F_NoEpdCutAllMass)->GetEntries()>0 ? ((TH1*)mH1F_NoEpdCutAllMass)->DrawNormalized("hist e same") : ((TH1*)mH1F_NoEpdCutAllMass)->DrawCopy("hist e same");
  h1allcutbutepd->SetStats(0);
  h1allcutbutepd->SetLineColor(kBlack);
  TH1* h1allcutepdch = ((TH1*)mH1F_EpdChAllMass)->GetEntries()>0 ? ((TH1*)mH1F_EpdChAllMass)->DrawNormalized("hist e same") : ((TH1*)mH1F_EpdChAllMass)->DrawCopy("hist e same");
  h1allcutepdch->SetStats(0);
  h1allcutepdch->SetLineColor(kGreen+2);
  std::stringstream ss_legname;
  ss_legname << "EpdNmip<"<<mEpdNmipCut;
  legpad12->AddEntry(h1allcutmass,ss_legname.str().c_str(),"fle");
  legpad12->AddEntry(h1allcutbutepd,"No EPD nmip cut","fle");
  ss_legname.str("");
  ss_legname << "EpdNmip>="<<mEpdNmipCut;
  legpad12->AddEntry(h1allcutepdch,ss_legname.str().c_str(),"fle");
  legpad12->Draw();
  
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintAllHistOneTrigger(TCanvas* canv, int trigidx, const char* savename) const
{
  canv->Clear();
  canv->Divide(4,3);

  canv->cd(1);
  mH1F_InvMassAllCuts->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(2);
  mH1F_Pi0MultAllCuts->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(3);
  mH1F_AllCuts_xF->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(4);
  mH1F_AllCuts_xFZoom->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(5);
  mH1F_AllCuts_Zgg->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(6);
  mH1F_AllCuts_Dgg->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(7);
  mH1F_AllCuts_Pi0En->UncheckedAt(trigidx)->Draw("hist e");
  canv->cd(8);
  mH2F_AllCuts_Pi0_massVen->UncheckedAt(trigidx)->Draw("colz");
  canv->cd(9);
  mH2F_AllCuts_Pi0_xfVen->UncheckedAt(trigidx)->Draw("colz");
  canv->cd(10);
  mH2F_AllCuts_Pi0_ptVeta->UncheckedAt(trigidx)->Draw("colz");
  canv->cd(11);
  mH2F_AllCuts_Pi0_etaVphi->UncheckedAt(trigidx)->Draw("colz");
  canv->cd(12);
  mH2F_AllCuts_Pi0_yVx->UncheckedAt(trigidx)->Draw("colz");
  
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintOneHistAllTrigger(TCanvas* canv, TObjArray* histarr, const char* drawoption, const char* savename) const
{
  canv->Clear();
  canv->Divide(3,2);

  TString opt(drawoption);
  canv->cd(1);
  histarr->UncheckedAt(0)->Draw(opt.Data());
  canv->cd(2);
  histarr->UncheckedAt(1)->Draw(opt.Data());
  canv->cd(3);
  histarr->UncheckedAt(2)->Draw(opt.Data());
  canv->cd(4);
  histarr->UncheckedAt(3)->Draw(opt.Data());
  canv->cd(5);
  histarr->UncheckedAt(4)->Draw(opt.Data());
  
  canv->Print(savename);
}

Int_t StMuFcsPi0TreeMaker::DrawEpdProjection(TCanvas* canvas, const char* savename)
{
  if( mEpdGeo==0 ){ return 0; }
  if( canvas==0 ) { return 0; }
  canvas->Clear();
  canvas->cd();
  canvas->DrawFrame(-100,-100,100,100);
  if( mEpdTileMap.size()>0 ){
    for( std::map<Int_t,TPolyLine*>::iterator polyitr=mEpdTileMap.begin(); polyitr!=mEpdTileMap.end(); ++polyitr ){
      polyitr->second->SetLineWidth(1);
      polyitr->second->SetLineColor(kBlack);
      polyitr->second->SetFillColorAlpha(kWhite,0);
    }
  }
  else{
    //loop over all west epd tiles and make the polygon need to do this every time otherwise root does not draw the fill correctly
    for(int i_pp=1; i_pp<=12; ++i_pp){     //Supersector runs [1,12]
      for( int i_tt=1; i_tt<=31; ++i_tt ){ //Tile number [1,31]
	TPolyLine* polyline = EpdTilePoly(i_pp,i_tt);
	if( polyline!=0 ){
	  polyline->SetLineColor(kBlack);
	  polyline->SetFillColorAlpha(kWhite,0);
	  std::pair<std::map<Int_t,TPolyLine*>::iterator,bool> itr = mEpdTileMap.emplace(100*i_pp+i_tt,polyline);
	  //std::cout << "|i_pp:"<<i_pp<<"|i_tt:"<<i_tt<<"|key:"<<100*i_pp+i_tt << std::endl;
	  if( ! (itr.second) ){ std::cout << "Could not create polyline map element" << std::endl; }
	}
      }
    }
  }

  //Create boxes for all FCS towers
  std::vector<TBox*> FcsTowers;
  double zepd = 375.0;
  //double zfcs=710.0+13.90+15.0;
  //double zr = zepd/zfcs;
  //float yoff[kFcsNDet]={100.0,100.0,-100.0,-100.0,100.0,100.0};
  //float xoff[kFcsNDet]={  0.0,  0.0,   0.0,   0.0,  0.0,  0.0};
  for( int idet=0; idet<kFcsNDet; idet++ ){
    if(idet<=kFcsEcalSouthDetId){ //ecal & hcal
      for(int id=0; id<mFcsDb->maxId(idet); id++){
	StThreeVectorF xyz = mFcsDb->getStarXYZ(idet,id);
	double wx=mFcsDb->getXWidth(idet);
	double wy=mFcsDb->getYWidth(idet);
	double zr = zepd/xyz.z();
	TBox* cell=new TBox(zr*(xyz.x()-wx/2.0), zr*(xyz.y()-wy/2.0), 
			    zr*(xyz.x()+wx/2.0), zr*(xyz.y()+wy/2.0));
	cell->SetFillColorAlpha(kWhite,0);
	cell->SetLineColorAlpha(kGray,0.5);
	cell->SetLineWidth(1);
	FcsTowers.push_back(cell);
      }
    }
  }

  /*std::map<Int_t,bool> fcsepdhits;
  for( unsigned int ihit=mMuFcsColl->indexOfFirstHit(4); ihit<mMuFcsColl->numberOfHits(); ++ihit ){
    StMuFcsHit* fcshit = mMuFcsColl->getHit(ihit);
    unsigned short det = fcshit->detectorId();
    unsigned short id = fcshit->id();
    int pp=0;
    int tt=0;
    mFcsDb->getEPDfromId(det,id,pp,tt);
    int fcsepdkey = 100*pp+tt;
    //std::cout << "|det:"<<det <<"|id:"<<id << "|fcsepdkey:"<<fcsepdkey << std::endl;
    if( fcsepdkey!=0 ){ fcsepdhits[fcsepdkey] = false; }
    //fcshit->energy()
  }*/
  //loop over all epd hits
  unsigned int nepdhits = 0;
  StSPtrVecEpdHit* epdhits = 0;
  if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
  else if( mEpdColl!=0 ){
    epdhits = &(mEpdColl->epdHits());
    nepdhits = epdhits->size();
  }
  else{ LOG_ERROR << "StMuEpdRun22QaMaker::FillEpdinfo() - If you see this error then there is a bug that is setting EPD hits improperly" << endm; return 0; }
  StMuEpdHit* muepdhit = 0;
  StEpdHit* epdhit = 0;
  //std::cout << "|nepdhits:"<<nepdhits << std::endl;
  for(unsigned int i=0; i<nepdhits; ++i ){
    if( mMuEpdHits!=0 ){ muepdhit = (StMuEpdHit*)mMuEpdHits->UncheckedAt(i); } //To match similar in StMuDstMaker->epdHit(int i)
    else if( epdhits!=0 ){ epdhit = (StEpdHit*)((*epdhits)[i]); }
    else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A VERY SERIOUS BUG IN THE CODE" << endm; return 0; } 
    //std::cout << "|i:"<<i << "|muepdhit:"<<muepdhit << "|epdhit:"<<epdhit << std::endl;
    int ew    = muepdhit!=0 ? muepdhit->side()    : epdhit->side();      //east=-1, west=1
    if( ew==-1 ){ continue; }
    int epdpp = muepdhit!=0 ? muepdhit->position(): epdhit->position();  //Supersector runs [1,12]
    int epdtt = muepdhit!=0 ? muepdhit->tile()    : epdhit->tile();      //Tile number [1,31]
    float nmip = muepdhit!=0 ? muepdhit->nMIP(): epdhit->nMIP();         //The ADC value of the hit divided by the MIP peak position; e.g. if nmip==1 then adc value sits at the MIP peak
    std::map<Int_t,TPolyLine*>::iterator epdhitit = mEpdTileMap.find(100*epdpp+epdtt);
    //std::map<Int_t,bool>::iterator fcsepdhitit = fcsepdhits.find(100*epdpp+epdtt);
    //if( fcsepdhitit!=fcsepdhits.end() ){ fcsepdhitit->second = true; }
    //else{ std::cout << "|FoundEpdHit not in Fcs Hit:"<<100*epdpp+epdtt << std::endl; }
    if( epdhitit!=mEpdTileMap.end() ){
      if( 0<nmip && nmip<mEpdNmipCut ){
	epdhitit->second->SetFillColorAlpha(kBlue,0.2);
	//std::cout << "  <nmip|ihit:"<<i << "|epdpp:"<<epdpp << "|epdtt:"<<epdtt <<"|nmip:"<<nmip << std::endl;
      }
      else if( nmip>mEpdNmipCut ){
	//Int_t color = kOrange+nmip;
	//if( nmip>10 ){ color = kOrange+10; }
	//epdhitit->second->SetFillColorAlpha(kRed,0.2);
	Int_t color = kOrange;
	epdhitit->second->SetFillColor(color);
	//std::cout << "  >nmip|ihit:"<<i << "|epdpp:"<<epdpp << "|epdtt:"<<epdtt <<"|nmip:"<<nmip << std::endl;
      }
      else{
	epdhitit->second->SetFillColor(kGray);
      }
    }
    else{
      std::cout << "Not in map!" << std::endl;
    }
  }

  //loop over all photon candidates
  std::vector<TMarker*> FcsMarkers;
  std::vector<TEllipse*> FcsPointEllipse;
  for( Int_t iph = 0; iph<mPhArr->GetEntriesFast(); ++iph ){
    if( !(mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh) ){ continue; }
    //std::cout << "|iph:"<<iph << "|iphnew:"<<iph-noldhits << std::endl;
    FcsPhotonCandidate* ph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(iph);
    if( ph==0 ){ std::cout << "==========I=CANNOT=BE=ZERO==========" << std::endl; return 0; }
    std::vector<Double_t> epdproj = ProjectToEpd(ph->mX,ph->mY,ph->mZ,mUseVertex);
    //std::cout << "|i:"<<iph <<"|x:"<<epdproj.at(0) << "|y:"<<epdproj.at(1) << std::endl;
    TMarker* fcsmarker = 0;
    TEllipse* ellipse = 0;
    if( ph->mFromCluster ){
      fcsmarker = new TMarker(epdproj.at(0),epdproj.at(1),33);
      fcsmarker->SetMarkerSize(1.5);
    }
    else{
      fcsmarker = new TMarker(epdproj.at(0),epdproj.at(1),30);
      fcsmarker->SetMarkerSize(2);
      ellipse = new TEllipse(epdproj.at(0),epdproj.at(1),10,10);
      int linesize=1;
      //A bit hacky but should work since size and energy is linearly increasing
      if( ph->mEn>10 ){ linesize=2; }
      if( ph->mEn>20 ){ linesize=3; }
      if( ph->mEn>30 ){ linesize=4; }
      if( ph->mEn>40 ){ linesize=5; }
      if( ph->mEn>50 ){ linesize=6; }
      if( ph->mEn>60 ){ linesize=7; }
      if( ph->mEn>70 ){ linesize=8; }
      if( ph->mEn>80 ){ linesize=9; }
      if( ph->mEn>90 ){ linesize=10; }
      ellipse->SetLineWidth(linesize);
      FcsPointEllipse.push_back(ellipse);
    }
    fcsmarker->SetMarkerColor(kBlack);
    
    //loop over all west epd tiles so that even if no hit recorded can use as a veto
    //Supersector runs [1,12]
    //Tile number [1,31]
    for(int i_pp=1; i_pp<=12; ++i_pp){
      for( int i_tt=1; i_tt<=31; ++i_tt ){
	if( ph->mEpdMatch[0] == (100*i_pp+i_tt) ){
	  //if( mEpdGeo->IsInTile(i_pp,i_tt,1,epdproj.at(0),epdproj.at(1)) ){
	  //std::cout << " + |iph:"<<iph << "|nmip:"<<ph->mEpdHitNmip[0] << "|epdkey:"<<100*i_pp+i_tt << std::endl;
	  std::map<Int_t,TPolyLine*>::iterator epdhitit = mEpdTileMap.find(100*i_pp+i_tt);
	  if( epdhitit!=mEpdTileMap.end() ){
	    if( ph->mEpdHitNmip[0]==0 ){
	      epdhitit->second->SetLineColor(kGreen+2);
	      epdhitit->second->SetFillColorAlpha(kGreen,0.5);
	      epdhitit->second->SetLineWidth(2);
	      fcsmarker->SetMarkerColor(kGreen+2);
	      if( ellipse!=0 ){ ellipse->SetLineColor(kGreen+2); }
	    }
	    else if( 0<ph->mEpdHitNmip[0] && ph->mEpdHitNmip[0]<mEpdNmipCut ){
	      epdhitit->second->SetLineColor(kBlue);
	      epdhitit->second->SetFillColorAlpha(kBlue,0.5);
	      epdhitit->second->SetLineWidth(2);
	      fcsmarker->SetMarkerColor(kBlue);
	      if( ellipse!=0 ){ ellipse->SetLineColor(kBlue); }
	      //std::cout << "<nmip|iph:"<<iph << "|i_pp:"<<i_pp << "|i_tt:"<<i_tt << "|nmip"<<ph->mEpdHitNmip[0] << std::endl;
	    }
	    else if( ph->mEpdHitNmip[0]>=mEpdNmipCut ){
	      epdhitit->second->SetLineColor(kRed);
	      epdhitit->second->SetFillColorAlpha(kRed,0.5);
	      epdhitit->second->SetLineWidth(2);
	      fcsmarker->SetMarkerColor(kRed);
	      if( ellipse!=0 ){ ellipse->SetLineColor(kRed); }
	    }
	    else{
	      epdhitit->second->SetLineColor(kGray);
	      fcsmarker->SetMarkerColor(kGray);
	      if( ellipse!=0 ){ ellipse->SetLineColor(kGray); }
	    }
	  }
	  else{
	    std::cout << "Not in map!" << std::endl;
	  }
	}
      }
    }
    FcsMarkers.push_back(fcsmarker);
  }

  //Draw Epd tiles first so they appear in the background
  for( std::map<Int_t,TPolyLine*>::iterator polyitr=mEpdTileMap.begin(); polyitr!=mEpdTileMap.end(); ++polyitr ){
    //std::cout << "|id:"<<polyitr->first << "|polyline:"<<polyitr->second << std::endl;
    polyitr->second->Draw("f");
    polyitr->second->Draw();
  }
  //Draw photon candidates next on top of the epd tiles
  for( unsigned int i=0; i<FcsMarkers.size(); ++i ){
    FcsMarkers.at(i)->Draw();
  }
  for( unsigned int i=0; i<FcsPointEllipse.size(); ++i ){
    FcsPointEllipse.at(i)->SetFillColorAlpha(kWhite,0);
    FcsPointEllipse.at(i)->Draw();
  }
  for( unsigned int i=0; i<FcsTowers.size(); ++i ){
    FcsTowers.at(i)->Draw("l");
  }

  /*for( std::map<Int_t,bool>::iterator fcsepditr=fcsepdhits.begin(); fcsepditr!=fcsepdhits.end(); ++fcsepditr ){
    if( fcsepditr->second==false ){ std::cout << "|NOHIT:"<<fcsepditr->first << std::endl; }
    }*/

  canvas->Print(savename);
  //Clean memory
  for( unsigned int i=0; i<FcsMarkers.size(); ++i ){
    delete FcsMarkers.at(i);
  }
  for( unsigned int i=0; i<FcsPointEllipse.size(); ++i ){
    delete FcsPointEllipse.at(i);
  }
  for( unsigned int i=0; i<FcsTowers.size(); ++i ){
    delete FcsTowers.at(i);
  }
  return 1;
}

void StMuFcsPi0TreeMaker::CheckInsideEpdTile(FcsPhotonCandidate* photon, Double_t projx, Double_t projy ) const
{
  //loop over all west epd tiles so that even if no hit recorded can use as a veto
  for(int i_pp=1; i_pp<=12; ++i_pp){     //Supersector runs [1,12]
    for( int i_tt=1; i_tt<=31; ++i_tt ){ //Tile number [1,31]
      if( mEpdGeo->IsInTile(i_pp,i_tt, 1, projx,projy) ){ //Only care about west EPD tiles; hence the '1'
	//TPolyLine* epd_ccw = EpdCCWOuterCorner(i_pp,i_tt);
	photon->mEpdHitNmip[0] = 0;
	photon->mEpdMatch[0] = 100*i_pp + i_tt;
	//std::cout << " + |projx:"<<projx << "|projy:"<<projy << "|nmip:"<< photon->mEpdHitNmip[0] << "|epdkey:"<<photon->mEpdMatch[0] << std::endl;
	break; //Inside match should be unique
      }
    }
  }
  if( photon->mEpdMatch[0]==0 ){ //If no intersection found it would be -1 so now check all the CCW adjacencies
    //std::cout << "   - |projx:"<<projx << "|projy:"<<projy << "|nmip:"<< photon->mEpdHitNmip[0] << "|epdkey:"<<photon->mEpdMatch[0] << std::endl;
    int ccwcounter = 1; //Should be 1 but Hack to check the drawing of the projections
    for(int i_pp=1; i_pp<=12; ++i_pp){     //Supersector runs [1,12]
      for( int i_tt=1; i_tt<=31; ++i_tt ){ //Tile number [1,31]
	TPolyLine* epd_outerccw = EpdCCWOuterCorner(i_pp,i_tt);
	if( TMath::IsInside( projx, projy, epd_outerccw->GetN(), epd_outerccw->GetX(), epd_outerccw->GetY() ) ){
	  //TPolyLine* epd_ccw = EpdCCWOuterCorner(i_pp,i_tt);
	  photon->mEpdHitNmip[ccwcounter] = 0;
	  photon->mEpdMatch[ccwcounter] = 100*i_pp + i_tt;
	  //std::cout << "     - |ccwcounter:"<<ccwcounter << "|nmip:"<< photon->mEpdHitNmip[ccwcounter] << "|epdkey:"<<photon->mEpdMatch[ccwcounter] << std::endl;
	  ++ccwcounter;
	}
      }
    }
  }
  
  int ncorners = 0;
  for( int icorner=0; icorner<5; ++icorner ){
    //std::cout << " + |projx:"<<projx << "|projy:"<<projy << "|epdkey:"<<photon->mEpdHitNmip[icorner] << std::endl;
    if( photon->mEpdMatch[icorner]!=0 ){ ++ncorners; }
  }
  
  if( ncorners>1 ){
    int bestcorner = 0;
    Double_t mindist = 999; //Pick some large distance so that the minimum will get set with first loop
    for( int icorner=0; icorner<5; ++icorner ){
      //Pick the best corner and set it to 0 value since the algorithm above only cares about the match in 0
      //std::cout << " + |projx:"<<projx << "|projy:"<<projy << "|nmip:"<< photon->mEpdHitNmip[icorner] << "|epdkey:"<<photon->mEpdMatch[icorner];
      if( photon->mEpdMatch[icorner]!=0 ){
	int epdpp = photon->mEpdMatch[icorner]/100;
	int epdtt = photon->mEpdMatch[icorner] - epdpp*100;
	TVector3 epdhitxyz = mEpdGeo->TileCenter(epdpp,epdtt,1);//1 for west
	Double_t distx = projx-epdhitxyz.x();
	Double_t disty = projy-epdhitxyz.y();
	Double_t dist = TMath::Sqrt(distx*distx+disty*disty);
	//std::cout << "|("<<epdhitxyz.x() << ","<<epdhitxyz.y() <<")|dx:"<<distx << "|dy:"<< disty << "|dist:"<<dist;
	if( dist<mindist ){ bestcorner = icorner; }
      }
      //else{std::cout << "|("<<0 << ","<<0 <<")"; }
      //std::cout << std::endl;
    }
    //std::cout << "   + |ncorners:"<<ncorners << std::endl;
    photon->mEpdMatch[0] = photon->mEpdMatch[bestcorner];
    photon->mEpdHitNmip[0] = 0;
  }
  
    /*
    //For all tiles except 1, 2, or 3 only check Outer CCW since this will cover the gap for all tiles except 1, 2, or 3
    TPolyLine* epd_outerccw = EpdCCWOuterCorner(i_pp,i_tt);
      if( TMath::IsInside( projx, projy, epd_outerccw->GetN(), epd_outerccw->GetX(), epd_outerccw->GetY() ) ){
	photon->mEpdHitNmip[1] = 0;
	photon->mEpdMatch[1] = 100*i_pp + i_tt;
      }
      //else if( i_tt==1 || i_tt==2 || i_tt==3 ){
      //For tiles 1, 2, and 3 need to check other than the outer CCW because of the pentagonal structure of tile 1 means that inner CCW, inner CW, and outer CW have a weird overlap that needs checking
      TPolyLine* epd_innerccw = EpdCCWInnerCorner(i_pp,i_tt);
      if( TMath::IsInside( projx, projy, epd_innerccw->GetN(), epd_innerccw->GetX(), epd_innerccw->GetY() ) ){
	photon->mEpdHitNmip[2] = 0;
	photon->mEpdMatch[2] = 100*i_pp + i_tt;
      }
      TPolyLine* epd_innercw = EpdCWInnerCorner(i_pp,i_tt);
      if( TMath::IsInside( projx, projy, epd_innercw->GetN(), epd_innercw->GetX(), epd_innercw->GetY() ) ){
	photon->mEpdHitNmip[3] = 0;
	photon->mEpdMatch[3] = 100*i_pp + i_tt;
      }
      TPolyLine* epd_outercw = EpdCWOuterCorner(i_pp,i_tt);
      if( TMath::IsInside( projx, projy, epd_outercw->GetN(), epd_outercw->GetX(), epd_outercw->GetY() ) ){
	photon->mEpdHitNmip[4] = 0;
	photon->mEpdMatch[4] = 100*i_pp + i_tt;
      }
    } //for i_tt
  }   //for i_pp
    */
}

TPolyLine* StMuFcsPi0TreeMaker::EpdTilePoly(short pp, short tt) const
{
  if( mEpdGeo==0 ){ return 0; }
  double x[5] = {0};
  double y[5] = {0};
  int ncorners = 0;
  short eastwest = 1;   //1 means EPD west
  mEpdGeo->GetCorners(pp,tt,eastwest,&ncorners,x,y);
  if( ncorners==0 ){ return 0; }
  std::vector<double> xvals;
  std::vector<double> yvals;
  for( int i=0; i<ncorners; ++i ){
    xvals.emplace_back(x[i]);
    yvals.emplace_back(y[i]);
    if( i==ncorners-1 ){
      xvals.emplace_back(x[0]);
      yvals.emplace_back(y[0]);
    }
  }
  //std::cout << "|pp:"<<pp << "|tt:"<<tt << "|n:"<<ncorners << "|";
  for( unsigned int j=0; j<xvals.size(); ++j ){
    //std::cout << "("<<xvals.at(j) << ","<<yvals.at(j) << ")|";
  }
  //std::cout << std::endl;
  TPolyLine* polyline = new TPolyLine(xvals.size(),xvals.data(),yvals.data());  //Equal sizes so shouldn't matter which is used
  return polyline;
}


//Makes the 2x2 
TPolyLine* StMuFcsPi0TreeMaker::EpdCCWOuterCorner(short pp, short tt) const
{
  //std::vector<Int_t> adjaenttiles = StMuEpdRun22QaMaker::GetAdjacentEpdIds(pp,tt);
  //std::vector<TPolyLine*> alllines;
  //alllines.push_back(EpdTilePoly(pp,tt)); //Start with center tile
  TPolyLine* main = EpdTilePoly(pp,tt);
  Int_t nmain = main->GetN();
  Double_t* xmain = main->GetX();
  Double_t* ymain = main->GetY();
  std::list<double> xvals;
  std::list<double> yvals;
  for( int i=0; i<nmain-1; ++i ){
    xvals.push_back(xmain[i]);
    yvals.push_back(ymain[i]);
  }
  std::list<double>::iterator xitr = xvals.begin();
  std::list<double>::iterator yitr = yvals.begin();
  /*
  for( ; xitr!=xvals.end() && yitr!=yvals.end(); ++xitr, ++yitr ){
    std::cout << " + |x:"<< *xitr << "|y:"<<*yitr << std::endl;
    }*/
  bool CCWIstt1 = false;
  xitr = xvals.begin();
  yitr = yvals.begin();
  Int_t adjpp = 0;
  Int_t adjtt = 0;
  StMuEpdRun22QaMaker::GetEpdTileCCW(pp,tt,adjpp,adjtt);
  //std::cout << "CCW|adjpp:"<<adjpp << "|adjtt:"<<adjtt << std::endl;
  if( adjpp!=0 && adjtt!=0 ){
    //Know there is something in the counter clockwise direction
    TPolyLine* adjline = EpdTilePoly(adjpp,adjtt);
    Int_t nadj = adjline->GetN();
    Double_t* adjxvals = adjline->GetX();
    Double_t* adjyvals = adjline->GetY();
    int lastcorner = 3;  //For rectangular tiles this is the index to use
    //std::cout << "|nadj:"<<nadj << std::endl;
    if( nadj==6 ){
      //For tile pp1
      CCWIstt1 = true;
      lastcorner = 4;  //For tt1 which is pentagonal this is the index to use
    }
    //if( nadj==5 ){  //5 is actually 4 corners since last element is the same as the start
      //std::cout << "HERECCW:"<<adjxvals[0] <<"|"<<adjyvals[0] << std::endl;
      //Since "insert" will insert before the iterator advance to second corner which is the outer CCW corner
      std::advance(xitr,1);
      std::advance(yitr,1);
      //Want to add CW inner edge first and then inner CCW edge
      xvals.insert(xitr,adjxvals[lastcorner]);
      yvals.insert(yitr,adjyvals[lastcorner]);
      xvals.insert(xitr,adjxvals[0]);
      yvals.insert(yitr,adjyvals[0]);
      xvals.insert(xitr,adjxvals[1]);
      yvals.insert(yitr,adjyvals[1]);
      xvals.insert(xitr,adjxvals[2]);
      yvals.insert(yitr,adjyvals[2]);
      if( CCWIstt1 ){
	//Add extra corner for tt1
	xvals.insert(xitr,adjxvals[3]);
	yvals.insert(yitr,adjyvals[3]);
      }
      //}
    delete adjline;
    
    /*for( std::list<double>::iterator xit=xvals.begin(), yit=yvals.begin(); xit!=xvals.end() && yit!=yvals.end(); ++xit, ++yit ){
      std::cout << " + |x:"<< *xit << "|y:"<<*yit << std::endl;
      }*/
  }
  StMuEpdRun22QaMaker::GetEpdTileOuter(pp,tt,adjpp,adjtt);
  //std::cout << "Outer|adjpp:"<<adjpp << "|adjtt:"<<adjtt << std::endl;
  if( adjpp!=0 && adjtt!=0 ){
    //Erase point on this corner
    //std::cout << "  - |x:"<<*xitr << "|y:"<<*yitr << std::endl;
    xitr = xvals.erase(xitr); //Gets set to next element (corner)
    yitr = yvals.erase(yitr);
    //std::cout << "AFTERERASE" << std::endl;
    /*for( std::list<double>::iterator xit=xvals.begin(), yit=yvals.begin(); xit!=xvals.end() && yit!=yvals.end(); ++xit, ++yit ){
      std::cout << " + |x:"<< *xit << "|y:"<<*yit << std::endl;
      }*/
    --xitr;   //Go back to previous corner
    --yitr;
    //std::cout << "  - |x:"<<*xitr << "|y:"<<*yitr << std::endl;
    //Know there is something above so add the points accordingly
    TPolyLine* adjline = EpdTilePoly(adjpp,adjtt);
    Int_t nadj = adjline->GetN();
    Double_t* adjxvals = adjline->GetX();
    Double_t* adjyvals = adjline->GetY();
    if( nadj==5 ){
      //std::cout << "HEREOUTER:"<<adjxvals[0] <<"|"<<adjyvals[0] << std::endl;
      //Since "insert" will insert before the iterator advance to third corner of original tile
      std::advance(xitr,1);
      std::advance(yitr,1);
      xvals.insert(xitr,adjxvals[0]);
      yvals.insert(yitr,adjyvals[0]);
      xvals.insert(xitr,adjxvals[1]);
      yvals.insert(yitr,adjyvals[1]);
      xvals.insert(xitr,adjxvals[2]);
      yvals.insert(yitr,adjyvals[2]);
      xvals.insert(xitr,adjxvals[3]);
      yvals.insert(yitr,adjyvals[3]);
    }
    else{ std::cout << "MAJOR ERROR:TILE1 CANNOT BE AN OUTER TILE" << std::endl; }
    delete adjline;
    /*
    for( std::list<double>::iterator xit=xvals.begin(), yit=yvals.begin(); xit!=xvals.end() && yit!=yvals.end(); ++xit, ++yit ){
      std::cout << " + |x:"<< *xit << "|y:"<<*yit << std::endl;
      }*/
  }
  StMuEpdRun22QaMaker::GetEpdTileOuterCCW(pp,tt,adjpp,adjtt);
  //std::cout << "OuterCCW|adjpp:"<<adjpp << "|adjtt:"<<adjtt << std::endl;
  if( adjpp!=0 && adjtt!=0 ){
    //Know there is something in the outer CCW position
    TPolyLine* adjline = EpdTilePoly(adjpp,adjtt);
    Int_t nadj = adjline->GetN();
    Double_t* adjxvals = adjline->GetX();
    Double_t* adjyvals = adjline->GetY();
    if( nadj==5 ){
      //std::cout << "HEREOUTERCCW:"<<adjxvals[0] <<"|"<<adjyvals[0] << std::endl;
      xitr = xvals.begin();
      yitr = yvals.begin();
       //Get to the points related to the inner corner and remove them
      std::advance(xitr,4);
      std::advance(yitr,4);
      //std::cout << "  - |x:"<<*xitr << "|y:"<<*yitr << std::endl;
      xitr = xvals.erase(xitr); //Gets set to next element (corner)
      yitr = yvals.erase(yitr);
      //std::cout << "AFTERERASE1" << std::endl;
      //std::cout << "  - |x:"<<*xitr << "|y:"<<*yitr << std::endl;
      xitr = xvals.erase(xitr);
      yitr = yvals.erase(yitr);
      if( CCWIstt1 ){
	//Delete extra corner when tt1 was added CCW
	xitr = xvals.erase(xitr);
	yitr = yvals.erase(yitr);
      }
      //std::cout << "AFTERERASE2" << std::endl;
      //std::cout << "  - |x:"<<*xitr << "|y:"<<*yitr << std::endl;
      /*for( std::list<double>::iterator xit=xvals.begin(), yit=yvals.begin(); xit!=xvals.end() && yit!=yvals.end(); ++xit, ++yit ){
	std::cout << " + |x:"<< *xit << "|y:"<<*yit << std::endl;
	}*/
      //Iterator has moved to next element which also needs to be deleted
      //Since "insert" will insert before the iterator it is now pointing to correct "top left corner"
      //Want to add clockwise edge first
      xvals.insert(xitr,adjxvals[0]);
      yvals.insert(yitr,adjyvals[0]);
      xvals.insert(xitr,adjxvals[1]);
      yvals.insert(yitr,adjyvals[1]);
      xvals.insert(xitr,adjxvals[2]);
      yvals.insert(yitr,adjyvals[2]);
      //xvals.push_back(adjxvals[3]); //Don't insert last point which is now inside the shape
      //yvals.push_back(adjyvals[3]);
    }
    delete adjline;
  }

  std::vector<double> newxvals;
  std::vector<double> newyvals;
  xitr=xvals.begin();
  yitr=yvals.begin();
  //std::cout << "=====final=====" << std::endl;
  for( ; xitr!=xvals.end() && yitr!=yvals.end(); ++xitr, ++yitr ){
    newxvals.push_back(*xitr);
    newyvals.push_back(*yitr);
    //std::cout << " + |x:"<< *xitr << "|y:"<<*yitr << std::endl;
  }
  //Add first element back in to close the polygon
  xitr = xvals.begin();
  yitr = yvals.begin();
  newxvals.push_back(*xitr);
  newyvals.push_back(*yitr);
    
  delete main;
  TPolyLine* newline = new TPolyLine(newxvals.size(),newxvals.data(),newyvals.data());  //Equal sizes so shouldn't matter which is used
  return newline;
}

TPolyLine* StMuFcsPi0TreeMaker::EpdCWOuterCorner(short pp, short tt) const
{
  int newtt = 0;
  int newpp = 0;
  StMuEpdRun22QaMaker::GetEpdTileCW(pp,tt,newpp,newtt);
  //std::cout << "|newpp:"<<newpp << "|newtt:"<<newtt << std::endl;
  if( newtt!=0 && newpp!=0 ){ return EpdCCWOuterCorner(newpp,newtt); }
  else{ return 0; }
}

TPolyLine* StMuFcsPi0TreeMaker::EpdCCWInnerCorner(short pp, short tt) const
{
  int newtt = 0;
  int newpp = 0;
  StMuEpdRun22QaMaker::GetEpdTileInner(pp,tt,newpp,newtt);
  //std::cout << "|pp:"<<pp << "|tt:"<<tt << "|newpp:"<<newpp << "|newtt:"<<newtt << std::endl;
  if( newpp!=0 && newtt!=0 ){ return EpdCCWOuterCorner(newpp,newtt); }
  else{
    if( tt==3 ){
      //Special for tile 3 to group tile 1, 2 and 3 together
      TPolyLine* tt1 = EpdTilePoly(pp,1);
      TPolyLine* tt2 = EpdTilePoly(pp,2);
      TPolyLine* tt3 = EpdTilePoly(pp,3);
      Double_t* xtt1 = tt1->GetX();
      Double_t* ytt1 = tt1->GetY();
      Double_t* xtt2 = tt2->GetX();
      Double_t* ytt2 = tt2->GetY();
      Double_t* xtt3 = tt3->GetX();
      Double_t* ytt3 = tt3->GetY();
      std::vector<double> xvals;
      std::vector<double> yvals;
      //Manually adding corners
      xvals.push_back(xtt1[0]);
      yvals.push_back(ytt1[0]);
      xvals.push_back(xtt1[1]);
      yvals.push_back(ytt1[1]);
      xvals.push_back(xtt2[0]);
      yvals.push_back(ytt2[0]);
      xvals.push_back(xtt2[1]);
      yvals.push_back(ytt2[1]);
      xvals.push_back(xtt2[2]);
      yvals.push_back(ytt2[2]);
      xvals.push_back(xtt3[1]);
      yvals.push_back(ytt3[1]);
      xvals.push_back(xtt3[2]);
      yvals.push_back(ytt3[2]);
      xvals.push_back(xtt3[3]);
      yvals.push_back(ytt3[3]);
      xvals.push_back(xtt1[3]);
      yvals.push_back(ytt1[3]);
      xvals.push_back(xtt1[4]);
      yvals.push_back(ytt1[4]);
      //Add last points back in to close the polygon
      xvals.push_back(xtt1[0]);
      yvals.push_back(ytt1[0]);
      
      delete tt1;
      delete tt2;
      delete tt3;
      TPolyLine* newline = new TPolyLine(xvals.size(),xvals.data(),yvals.data());  //Equal sizes so shouldn't matter which is used
      return newline;      
    }
    if( tt==1 ){
      //For tile 1 only add CCW tile
      TPolyLine* main = EpdTilePoly(pp,tt);
      Int_t nmain = main->GetN();
      Double_t* xmain = main->GetX();
      Double_t* ymain = main->GetY();
      std::list<double> xvals;
      std::list<double> yvals;
      for( int i=0; i<nmain-1; ++i ){
	xvals.push_back(xmain[i]);
	yvals.push_back(ymain[i]);
      }
      std::list<double>::iterator xitr = xvals.begin();
      std::list<double>::iterator yitr = yvals.begin();
      xitr = xvals.begin();
      yitr = yvals.begin();
      bool CCWIstt1 = false;
      Int_t adjpp = 0;
      Int_t adjtt = 0;
      StMuEpdRun22QaMaker::GetEpdTileCCW(pp,tt,adjpp,adjtt);
      if( adjpp!=0 && adjtt!=0 ){
	//Know there is something in the counter clockwise direction
	TPolyLine* adjline = EpdTilePoly(adjpp,adjtt);
	Int_t nadj = adjline->GetN();
	Double_t* adjxvals = adjline->GetX();
	Double_t* adjyvals = adjline->GetY();
	int lastcorner = 3;  //For rectangular tiles this is the index to use
	//std::cout << "|nadj:"<<nadj << std::endl;
	if( nadj==6 ){
	  //For tile pp1
	  CCWIstt1 = true;
	  lastcorner = 4;  //For tt1 which is pentagonal this is the index to use
	}
	std::advance(xitr,1);
	std::advance(yitr,1);
	//Want to add CW inner edge first and then inner CCW edge
	xvals.insert(xitr,adjxvals[lastcorner]);
	yvals.insert(yitr,adjyvals[lastcorner]);
	xvals.insert(xitr,adjxvals[0]);
	yvals.insert(yitr,adjyvals[0]);
	xvals.insert(xitr,adjxvals[1]);
	yvals.insert(yitr,adjyvals[1]);
	xvals.insert(xitr,adjxvals[2]);
	yvals.insert(yitr,adjyvals[2]);
	if( CCWIstt1 ){
	  //Add extra corner for tt1
	  xvals.insert(xitr,adjxvals[3]);
	  yvals.insert(yitr,adjyvals[3]);
	}
	delete adjline;
      }
      std::vector<double> newxvals;
      std::vector<double> newyvals;
      xitr=xvals.begin();
      yitr=yvals.begin();
      //std::cout << "=====final=====" << std::endl;
      for( ; xitr!=xvals.end() && yitr!=yvals.end(); ++xitr, ++yitr ){
	newxvals.push_back(*xitr);
	newyvals.push_back(*yitr);
	//std::cout << " + |x:"<< *xitr << "|y:"<<*yitr << std::endl;
      }
      //Add first element back in to close the polygon
      xitr = xvals.begin();
      yitr = yvals.begin();
      newxvals.push_back(*xitr);
      newyvals.push_back(*yitr);
    
      delete main;
      TPolyLine* newline = new TPolyLine(newxvals.size(),newxvals.data(),newyvals.data());  //Equal sizes so shouldn't matter which is used
      return newline;
    }
  }
  return 0;
}

TPolyLine* StMuFcsPi0TreeMaker::EpdCWInnerCorner(short pp, short tt) const
{
  int newtt = 0;
  int newpp = 0;
  StMuEpdRun22QaMaker::GetEpdTileInnerCW(pp,tt,newpp,newtt);
  //std::cout << "|pp:"<<pp << "|tt:"<<tt << "|newpp:"<<newpp << "|newtt:"<<newtt << std::endl;
  if( newpp!=0 && newtt!=0 ){ return EpdCCWOuterCorner(newpp,newtt); }
  else{
    if( tt==2 ){
      //Special for tile 2 to group tile 1, 2 and 3 together
      TPolyLine* tt1 = EpdTilePoly(pp,1);
      TPolyLine* tt2 = EpdTilePoly(pp,2);
      TPolyLine* tt3 = EpdTilePoly(pp,3);
      Double_t* xtt1 = tt1->GetX();
      Double_t* ytt1 = tt1->GetY();
      Double_t* xtt2 = tt2->GetX();
      Double_t* ytt2 = tt2->GetY();
      Double_t* xtt3 = tt3->GetX();
      Double_t* ytt3 = tt3->GetY();
      std::vector<double> xvals;
      std::vector<double> yvals;
      //Manually adding corners
      xvals.push_back(xtt1[0]);
      yvals.push_back(ytt1[0]);
      xvals.push_back(xtt1[1]);
      yvals.push_back(ytt1[1]);
      xvals.push_back(xtt2[0]);
      yvals.push_back(ytt2[0]);
      xvals.push_back(xtt2[1]);
      yvals.push_back(ytt2[1]);
      xvals.push_back(xtt2[2]);
      yvals.push_back(ytt2[2]);
      xvals.push_back(xtt3[1]);
      yvals.push_back(ytt3[1]);
      xvals.push_back(xtt3[2]);
      yvals.push_back(ytt3[2]);
      xvals.push_back(xtt3[3]);
      yvals.push_back(ytt3[3]);
      xvals.push_back(xtt1[3]);
      yvals.push_back(ytt1[3]);
      xvals.push_back(xtt1[4]);
      yvals.push_back(ytt1[4]);
      //Add last points back in to close the polygon
      xvals.push_back(xtt1[0]);
      yvals.push_back(ytt1[0]);
      
      delete tt1;
      delete tt2;
      delete tt3;
      TPolyLine* newline = new TPolyLine(xvals.size(),xvals.data(),yvals.data());  //Equal sizes so shouldn't matter which is used
      return newline;      
    }
    if( tt==1 ){
      //For tile 1 only add CC tile
      TPolyLine* main = EpdTilePoly(pp,tt);
      Int_t nmain = main->GetN();
      Double_t* xmain = main->GetX();
      Double_t* ymain = main->GetY();
      std::list<double> xvals;
      std::list<double> yvals;
      for( int i=0; i<nmain-1; ++i ){
	xvals.push_back(xmain[i]);
	yvals.push_back(ymain[i]);
      }
      std::list<double>::iterator xitr = xvals.end();
      std::list<double>::iterator yitr = yvals.end();
      xitr = xvals.end();
      yitr = yvals.end();
      Int_t adjpp = 0;
      Int_t adjtt = 0;
      StMuEpdRun22QaMaker::GetEpdTileCW(pp,tt,adjpp,adjtt);
      if( adjpp!=0 && adjtt!=0 ){
	//Know there is something in the clockwise direction and it is tt 1 with 5 corners based on adjacent mapping
	TPolyLine* adjline = EpdTilePoly(adjpp,adjtt);
	//Int_t nadj = adjline->GetN();
	Double_t* adjxvals = adjline->GetX();
	Double_t* adjyvals = adjline->GetY();
	//advance to last corner
	--xitr;
	--yitr;
	//Want to add CCW outer edge first and then outer CW edge and so on
	xvals.insert(xitr,adjxvals[1]);
	yvals.insert(yitr,adjyvals[1]);
	xvals.insert(xitr,adjxvals[2]);
	yvals.insert(yitr,adjyvals[2]);
	xvals.insert(xitr,adjxvals[3]);
	yvals.insert(yitr,adjyvals[3]);
	xvals.insert(xitr,adjxvals[4]);
	yvals.insert(yitr,adjyvals[4]);
	xvals.insert(xitr,adjxvals[0]);
	yvals.insert(yitr,adjyvals[0]);
	delete adjline;
      }
      std::vector<double> newxvals;
      std::vector<double> newyvals;
      xitr=xvals.begin();
      yitr=yvals.begin();
      //std::cout << "=====final=====" << std::endl;
      for( ; xitr!=xvals.end() && yitr!=yvals.end(); ++xitr, ++yitr ){
	newxvals.push_back(*xitr);
	newyvals.push_back(*yitr);
	//std::cout << " + |x:"<< *xitr << "|y:"<<*yitr << std::endl;
      }
      //Add first element back in to close the polygon
      xitr = xvals.begin();
      yitr = yvals.begin();
      newxvals.push_back(*xitr);
      newyvals.push_back(*yitr);
    
      delete main;
      TPolyLine* newline = new TPolyLine(newxvals.size(),newxvals.data(),newyvals.data());  //Equal sizes so shouldn't matter which is used
      return newline;
    }
  }
  return 0;
}

Int_t StMuFcsPi0TreeMaker::GetColor(Double_t Value, Double_t MinVal, Double_t MaxVal)
{
  Double_t MinHue = 275.0;
  Double_t MaxHue = 0.0;
  Float_t Red,Green,Blue;
  if( Value < MinVal ){ TColor::HSV2RGB(MinHue,1.0,1.0, Red,Green,Blue); }
  else if( Value > MaxVal ){ TColor::HSV2RGB(MaxHue,1.0,1.0, Red,Green,Blue); }
  else
    {
      Double_t percent = (Value-MinVal)/(MaxVal-MinVal);
      Float_t Hue = ((MaxHue-MinHue)*percent)+MinHue;
      TColor::HSV2RGB(Hue,1.0,1.0, Red,Green,Blue);
    }
    return TColor::GetColor( Red, Green, Blue );
}

void StMuFcsPi0TreeMaker::DrawSelectEpdAdjTiles(TCanvas* canv) const
{
  if( mEpdGeo==0 ){ return; }
  if( canv==0 ) { return; }
  std::vector<TPolyLine*> lines;
  //std::vector<TPolyLine*> extralines;

  canv->Clear();
  canv->cd();
  canv->DrawFrame(-30,0,30,50);
  TPolyLine* line_pp1tt8     = EpdTilePoly(1,8);
  TPolyLine* line_pp1tt7     = EpdTilePoly(1,7);
  TPolyLine* outccw_pp1tt8   = EpdCCWOuterCorner(1,8);
  TPolyLine* outcw_pp1tt8    = EpdCWOuterCorner(1,8);
  TPolyLine* innerccw_pp1tt8 = EpdCCWInnerCorner(1,8);
  TPolyLine* innercw_pp1tt8  = EpdCWInnerCorner(1,8);
  line_pp1tt8->SetLineColor(kBlack);
  line_pp1tt8->SetLineWidth(2);
  line_pp1tt7->SetLineColor(kBlack);
  outccw_pp1tt8->SetLineColor(kBlue);
  outcw_pp1tt8->SetLineColor(kGreen+2);
  innerccw_pp1tt8->SetLineColor(kRed);
  innercw_pp1tt8->SetLineColor(kMagenta);
  line_pp1tt8->Draw();
  line_pp1tt7->Draw();
  outccw_pp1tt8->Draw();
  outcw_pp1tt8->Draw();
  innerccw_pp1tt8->Draw();
  innercw_pp1tt8->Draw();
  canv->Print("EpdAdj_pp1tt8.png");

  canv->Clear();
  canv->cd();
  canv->DrawFrame(-15,0,10,15);
  TPolyLine* line_pp1tt1     = EpdTilePoly(1,1);
  TPolyLine* line_pp1tt2     = EpdTilePoly(1,2);
  TPolyLine* line_pp1tt3     = EpdTilePoly(1,3);
  TPolyLine* outccw_pp1tt1   = EpdCCWOuterCorner(1,1);
  TPolyLine* outcw_pp1tt1    = EpdCWOuterCorner (1,1);
  TPolyLine* innerccw_pp1tt1 = EpdCCWInnerCorner(1,1);
  TPolyLine* innercw_pp1tt1  = EpdCWInnerCorner (1,1);
  TPolyLine* outccw_pp1tt2   = EpdCCWOuterCorner(1,2);
  TPolyLine* outcw_pp1tt2    = EpdCWOuterCorner (1,2);
  TPolyLine* innerccw_pp1tt2 = EpdCCWInnerCorner(1,2);
  TPolyLine* innercw_pp1tt2  = EpdCWInnerCorner (1,2);
  TPolyLine* outccw_pp1tt3   = EpdCCWOuterCorner(1,3);
  TPolyLine* outcw_pp1tt3    = EpdCWOuterCorner (1,3);
  TPolyLine* innerccw_pp1tt3 = EpdCCWInnerCorner(1,3);
  TPolyLine* innercw_pp1tt3  = EpdCWInnerCorner (1,3);
  line_pp1tt1->SetLineColor(kBlack);
  line_pp1tt2->SetLineColor(kBlack);
  line_pp1tt3->SetLineColor(kBlack);
  outccw_pp1tt1->SetLineColor(kBlue);
  outcw_pp1tt1->SetLineColor(kGreen+2);
  innerccw_pp1tt1->SetLineColor(kRed);
  innercw_pp1tt1->SetLineColor(kMagenta);
  outccw_pp1tt2->SetLineColor(kBlue);
  outcw_pp1tt2->SetLineColor(kGreen+2);
  innerccw_pp1tt2->SetLineColor(kRed);
  innercw_pp1tt2->SetLineColor(kMagenta);
  outccw_pp1tt3->SetLineColor(kBlue);
  outcw_pp1tt3->SetLineColor(kGreen+2);
  innerccw_pp1tt3->SetLineColor(kRed);
  innercw_pp1tt3->SetLineColor(kMagenta);

  line_pp1tt1->SetLineWidth(2);
  line_pp1tt1->Draw();
  line_pp1tt2->SetLineWidth(1);
  line_pp1tt2->Draw();
  line_pp1tt3->SetLineWidth(1);
  line_pp1tt3->Draw();
  outccw_pp1tt1->Draw();
  outcw_pp1tt1->Draw();
  innerccw_pp1tt1->Draw();
  innercw_pp1tt1->Draw();
  canv->Print("EpdAdj_pp1tt1.png");

  canv->Clear();
  canv->cd();
  canv->DrawFrame(-15,0,10,25);
  line_pp1tt1->SetLineWidth(1);
  line_pp1tt1->Draw();
  line_pp1tt2->SetLineWidth(2);
  line_pp1tt2->Draw();
  line_pp1tt3->SetLineWidth(1);
  line_pp1tt3->Draw();
  outccw_pp1tt2->Draw();
  outcw_pp1tt2->Draw();
  innerccw_pp1tt2->Draw();
  innercw_pp1tt2->Draw();
  canv->Print("EpdAdj_pp1tt2.png");

  canv->Clear();
  canv->cd();
  canv->DrawFrame(-15,0,10,25);
  line_pp1tt1->SetLineWidth(1);
  line_pp1tt1->Draw();
  line_pp1tt2->SetLineWidth(1);
  line_pp1tt2->Draw();
  line_pp1tt3->SetLineWidth(2);
  line_pp1tt3->Draw();
  outccw_pp1tt3->Draw();
  outcw_pp1tt3->Draw();
  innerccw_pp1tt3->Draw();
  innercw_pp1tt3->Draw();
  canv->Print("EpdAdj_pp1tt3.png");

  //std::cout << "|outccw_pp1tt3:"<<outccw_pp1tt1 << "|outcw_pp1tt3:"<<outcw_pp1tt1 << "|innerccw_pp1tt3:"<<innerccw_pp1tt1 << "|innercw_pp1tt3:"<<innercw_pp1tt1 << std::endl;
  canv->Clear();
  canv->cd();
  canv->DrawFrame(-100,-100,100,100);
  Int_t itile=0;
  for(int i_pp=1; i_pp<=12; ++i_pp){     //Supersector runs [1,12]
    for( int i_tt=1; i_tt<=31; ++i_tt ){ //Tile number [1,31]
      TPolyLine* lineouterccw = EpdCCWOuterCorner(i_pp,i_tt);
      if( lineouterccw!=0 ){
	Int_t color = GetColor(itile++,0,12*31);
	lineouterccw->SetLineColor(color);
	lineouterccw->SetFillColorAlpha(kWhite,0);
	lines.push_back(lineouterccw);
      }
    }
  }
  //std::cout << " ===== |itile:"<<itile << " ===== " << std::endl;
  canv->SaveAs("EpdAllAdjMap.png");
  delete line_pp1tt8;
  delete line_pp1tt7;
  delete outccw_pp1tt8;
  delete outcw_pp1tt8;
  delete innerccw_pp1tt8;
  delete innercw_pp1tt8;
  delete line_pp1tt1;
  delete line_pp1tt2;
  delete line_pp1tt3;
  delete outccw_pp1tt1;
  delete outcw_pp1tt1;
  delete innerccw_pp1tt1;
  delete innercw_pp1tt1;
  delete outccw_pp1tt2;
  delete outcw_pp1tt2;
  delete innerccw_pp1tt2;
  delete innercw_pp1tt2;
  delete outccw_pp1tt3;
  delete outcw_pp1tt3;
  delete innerccw_pp1tt3;
  delete innercw_pp1tt3;
  for( unsigned int i=0; i<lines.size(); ++i ){
    delete lines.at(i);
  }
}

void StMuFcsPi0TreeMaker::DrawAllEpdAdjTiles(TCanvas* canv, Int_t tt) const
{
  if( tt<1 || tt>31 ){ return; }
  canv->Clear();
  canv->cd();
  canv->Divide(4,3);
  std::vector<TPolyLine*> lines;
  for( int i_pp=1; i_pp<=12; ++i_pp ){
    TPolyLine* tile_main  = EpdTilePoly(i_pp,tt);
    TPolyLine* outerccw   = EpdCCWOuterCorner(i_pp,tt);
    TPolyLine* outercw    = EpdCWOuterCorner(i_pp,tt);
    TPolyLine* innerccw   = EpdCCWInnerCorner(i_pp,tt);
    TPolyLine* innercw    = EpdCWInnerCorner(i_pp,tt);

    //To get right canvas dimesnosns
    double xlow=999;
    double ylow=999;
    double xhigh=-999;
    double yhigh=-999;
    if( tile_main!=0 ){
      Int_t n=tile_main->GetN();
      double* x=tile_main->GetX();
      double* y=tile_main->GetY();
      for( int i=0; i<n; ++i ){
	if( xlow>x[i] ) { xlow =x[i]; }
	if( ylow>y[i] ) { ylow =y[i]; }
	if( xhigh<x[i] ){ xhigh=x[i]; }
	if( yhigh<y[i] ){ yhigh=y[i]; }
      }
    }
    if( outerccw!=0 ){
      Int_t n=outerccw->GetN();
      double* x=outerccw->GetX();
      double* y=outerccw->GetY();
      for( int i=0; i<n; ++i ){
	if( xlow>x[i] ) { xlow =x[i]; }
	if( ylow>y[i] ) { ylow =y[i]; }
	if( xhigh<x[i] ){ xhigh=x[i]; }
	if( yhigh<y[i] ){ yhigh=y[i]; }
      }
    }
    if( outercw!=0 ){
      Int_t n=outercw->GetN();
      double* x=outercw->GetX();
      double* y=outercw->GetY();
      for( int i=0; i<n; ++i ){
	if( xlow>x[i] ) { xlow =x[i]; }
	if( ylow>y[i] ) { ylow =y[i]; }
	if( xhigh<x[i] ){ xhigh=x[i]; }
	if( yhigh<y[i] ){ yhigh=y[i]; }
      }
    }
    if( innerccw!=0 ){
      Int_t n=innerccw->GetN();
      double* x=innerccw->GetX();
      double* y=innerccw->GetY();
      for( int i=0; i<n; ++i ){
	if( xlow>x[i] ) { xlow =x[i]; }
	if( ylow>y[i] ) { ylow =y[i]; }
	if( xhigh<x[i] ){ xhigh=x[i]; }
	if( yhigh<y[i] ){ yhigh=y[i]; }
      }
    }
    if( innercw!=0 ){
      Int_t n=innercw->GetN();
      double* x=innercw->GetX();
      double* y=innercw->GetY();
      for( int i=0; i<n; ++i ){
	if( xlow>x[i] ) { xlow =x[i]; }
	if( ylow>y[i] ) { ylow =y[i]; }
	if( xhigh<x[i] ){ xhigh=x[i]; }
	if( yhigh<y[i] ){ yhigh=y[i]; }
      }
    }
    
    canv->cd(i_pp)->DrawFrame(xlow-4,ylow-4,xhigh+4,yhigh+4);
    tile_main->SetLineColor(kBlack);
    tile_main->SetLineWidth(2);
    tile_main->Draw();
    if( outerccw!=0 ){
      outerccw->SetLineColor(kBlue);
      outerccw->Draw();
    }
    if( outercw!=0 ){
      outercw->SetLineColor(kGreen+2);
      outercw->Draw();
    }
    if( innerccw!=0 ){
      innerccw->SetLineColor(kRed);
      innerccw->Draw();
    }
    if( innercw!=0 ){
      innercw->SetLineColor(kMagenta);
      innercw->Draw();
    }
    lines.push_back(tile_main);
    lines.push_back(outerccw);
    lines.push_back(outercw);
    lines.push_back(innerccw);
    lines.push_back(innercw);
  }
  std::stringstream ss_savename;
  ss_savename << "EpdAdjAllPpTt"<< tt <<".png";
  canv->Print(ss_savename.str().c_str());
  for( unsigned int i=0; i<lines.size(); ++i ){
    delete lines.at(i);
  }  
}

