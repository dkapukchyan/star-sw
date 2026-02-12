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

StMuFcsPi0TreeMaker::StMuFcsPi0TreeMaker(const Char_t* name) : StMuFcsTreeMaker(name)
{
  memset(mH2F_NPi0Inc_xfVphi, 0,sizeof(mH2F_NPi0Inc_xfVphi));
  memset(mH2F_NPi0Bg1_xfVphi, 0,sizeof(mH2F_NPi0Bg1_xfVphi));
  memset(mH2F_NPi0Bg2_xfVphi, 0,sizeof(mH2F_NPi0Bg2_xfVphi));
}

StMuFcsPi0TreeMaker::~StMuFcsPi0TreeMaker()
{
}

UInt_t StMuFcsPi0TreeMaker::LoadHists( TFile* file )
{
  if( mHists==0 ){ return 0; }
  /*
  if( othermkr!=0 ){
    this->CopyHists(othermkr);
    }*/
  UInt_t loaded = StMuFcsTreeMaker::LoadHists(file);
  //UInt_t loaded = 0;

  //mH2F_foundVvertex = (TH1*)mHists->FindObject("H2F_foundVvertex");

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

Int_t StMuFcsPi0TreeMaker::Make()
{
  return StMuFcsTreeMaker::Make();
}

Int_t StMuFcsPi0TreeMaker::Make_TssaAna(TClonesArray* pointpairs)
{
  std::cout << this->ClassName() << "|Start:Make_TssaAna()" << std::endl;
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
  for( int i=0; i<mPhPairArr->GetEntriesFast(); ++i ){
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

/*  
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
*/

/*void StMuFcsPi0TreeMaker::MergeForTssa( TH1* totalhistinc[][2], TH1* totalhistbg1[][2], TH1* totalhistbg2[][2], TH3* mergedinvmass, TH1* mergedpolblue, TH1* mergedpolyell, TH1* mergedpolblueerr, TH1* mergedpolyellerr )
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
  }*/

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

/*
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
    
}
*/


