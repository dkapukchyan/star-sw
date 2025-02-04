#include "StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEventTypes.h"
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

StMuFcsPi0TreeMaker::StMuFcsPi0TreeMaker(const Char_t* name) : StMaker(name)
{
  mSpinRndm.SetSeed(0);
  memset(mTriggers,0,sizeof(mTriggers));
  memset(mH1F_InvMassEpdCuts,0,sizeof(mH1F_InvMassEpdCuts));
  //memset(mH1F_InvMassZggCuts,0,sizeof(mH1F_InvMassZggCuts));
  //memset(mH1F_InvMassPtCuts, 0,sizeof(mH1F_InvMassPtCuts));
  //memset(mH1F_InvMassAllCutsByEnByPhi, 0,sizeof(mH1F_mH1F_InvMassAllCutsByEnByPhi));
  //memset(mH1F_NPi0ByEnByPhi, 0,sizeof(mH1F_NPi0ByEnByPhi));
  memset(mH2F_NPi0_enVphi, 0,sizeof(mH2F_NPi0_enVphi));
}

StMuFcsPi0TreeMaker::~StMuFcsPi0TreeMaker()
{
  delete mEvtInfo;
  delete mPhArr;
  delete mPi0Arr;
  delete mPi0Tree;
  if( mInternalHists ){ delete mHists; } //This deletes the file too
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
  if( mHists!=0 ){ LOG_WARN << "StMuFcsRun22QaMaker::setHistManager() - HistManager exists and is external - no changes made" << endm; return; }
  else{ mHists = hm; }
}

UInt_t StMuFcsPi0TreeMaker::LoadHists( TFile* file )
{
  if( mHists==0 ){ return 0; }
  UInt_t loaded = 0;
  loaded += mHists->AddH1F(file,mH1F_Entries,"H1_Entries", "Number of entries;", 3,-0.5, 2.5);
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
  
  loaded += mHists->AddH1F(file,mH1F_BestPi0Mass,"H1F_BestPi0Mass","Pi0s with photon combination of two highest energy pairs;Invariant Mass (GeV);", 500,0,1);
  //loaded += mHists->AddH2F(file,mH2F_Pi0HeatMap,"H2F_Pi0HeatMap","STAR x,y locations for Pi0s;x (cm);y (cm)", 500,-250,250, 500,-250,250);
  loaded += mHists->AddH1F(file,mH1F_PointMult,"H1F_PointMult","Point Multiplicity with only an energy cut;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_BestPi0Zgg,"H1F_BestPi0Zgg","Zgg of Pi0s with photon combination using highest energy pairs;Zgg;", 100,0,1);
  loaded += mHists->AddH1F(file,mH1F_BestPi0Phi,"H1F_BestPi0Phi","Phi of Pi0s with photon combination using highest energy pairs;Phi;", 100,-TMath::Pi(),TMath::Pi());
  loaded += mHists->AddH1F(file,mH1F_BestPi0Eta,"H1F_BestPi0Eta","Eta of Pi0s with photon combination using highest energy pairs;Eta;", 100,1,5.5);
  loaded += mHists->AddH1F(file,mH1F_BestPi0En,"H1F_BestPi0En","Energy of Pi0s with photon combination using highest energy pairs;Energy (GeV)", 1000,0,100);
  loaded += mHists->AddH1F(file,mH1F_BestPi0Pt,"H1F_BestPi0Pt","Pt of Pi0s with photon combination using highest energy pairs;Pt (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_AllPointPairMass,"H1F_AllPointsPi0Mass","Invariant mass of all point pair combinations;Invariant Mass (GeV);", 500,0,1);
  
  loaded += mHists->AddH1F(file,mH1F_EpdPhInvMass,"H1F_EpdPhInvMass","Pi0s with highest energy pairs from EPD cut photons;Invariant Mass (GeV);", 500,0,1);
  //loaded += mHists->AddH2F(file,mH2F_EpdPhHeatMap,"H2F_EpdPhHeatMap","STAR x,y locations for Pi0s with EPD cut photons;x (cm);y (cm)", 500,-250,250, 500,-250,250);
  loaded += mHists->AddH1F(file,mH1F_EpdPhPointMult,"H1F_EpdPhPointMult","Point Multiplicity with an energy cut and EPD cut on photon;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_EpdPhZgg,"H1F_EpdPhZgg","Zgg of Pi0s using highest energy pairs and Epd Cut Photons;Zgg;", 100,0,1);
  loaded += mHists->AddH1F(file,mH1F_EpdPhPhi,"H1F_EpdPhPhi","Phi of Pi0s using highest energy pairs and Epd Cut Photons;Phi;", 100,-TMath::Pi(),TMath::Pi());
  loaded += mHists->AddH1F(file,mH1F_EpdPhEta,"H1F_EpdPhEta","Eta of Pi0s using highest energy pairs and Epd Cut Photons;Eta;", 100,1,5.5);
  loaded += mHists->AddH1F(file,mH1F_EpdPhEn,"H1F_EpdPhEn","Energy of Pi0s using highest energy pairs and Epd Cut Photons;Energy (GeV)", 1000,0,100);
  loaded += mHists->AddH1F(file,mH1F_EpdPhPt,"H1F_EpdPhPt","Pt of Pi0s using highest energy pairs and Epd Cut Photons;Pt (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_EpdPhAllPoints,"H1F_EpdPhAllPoints","Invariant mass of all point pair combinations with Epd Cut Photons;Invariant Mass (GeV);", 500,0,1);


  loaded += mHists->AddH1F(file,mH1F_EpdChInvMass,"H1F_EpdChInvMass","Pi0s with highest energy pairs from EPD cut photons;Invariant Mass (GeV);", 500,0,1);
  //loaded += mHists->AddH2F(file,mH2F_EpdChHeatMap,"H2F_EpdChHeatMap","STAR x,y locations for Pi0s with EPD cut photons;x (cm);y (cm)", 500,-250,250, 500,-250,250);
  loaded += mHists->AddH1F(file,mH1F_EpdChPointMult,"H1F_EpdChPointMult","Point Multiplicity with an energy cut and EPD cut on photon;Point Multiplicity", 30,0,30);
  loaded += mHists->AddH1F(file,mH1F_EpdChZgg,"H1F_EpdChZgg","Zgg of Pi0s using highest energy pairs and Epd Cut Charged;Zgg;", 100,0,1);
  loaded += mHists->AddH1F(file,mH1F_EpdChPhi,"H1F_EpdChPhi","Phi of Pi0s using highest energy pairs and Epd Cut Charged;Phi;", 100,-TMath::Pi(),TMath::Pi());
  loaded += mHists->AddH1F(file,mH1F_EpdChEta,"H1F_EpdChEta","Eta of Pi0s using highest energy pairs and Epd Cut Charged;Eta;", 100,1,5.5);
  loaded += mHists->AddH1F(file,mH1F_EpdChEn,"H1F_EpdChEn","Energy of Pi0s using highest energy pairs and Epd Cut Charged;Energy (GeV)", 1000,0,100);
  loaded += mHists->AddH1F(file,mH1F_EpdChPt,"H1F_EpdChPt","Pt of Pi0s using highest energy pairs and Epd Cut Charged;Pt (GeV)", 100,0,10);
  loaded += mHists->AddH1F(file,mH1F_EpdChAllPoints,"H1F_EpdChAllPoints","Invariant mass of all point pair combinations with Epd Cut Charged;Invariant Mass (GeV);", 500,0,1); //This makes it such that this bin size is twice that of the 0,1 range with 500 bins

  //loaded += mHists->AddH1F(file,mH1F_2ndBestPi0Mass,"H1F_2ndBestPi0Mass", "Second best candidate for pi0 mass;Invariant Mass (GeV)", 500,0,1 );
  //loaded += mHists->AddH1F(file,mH1F_3rdBestPi0Mass,"H1F_3rdBestPi0Mass", "Third best candidate for pi0 mass;Invariant Mass (GeV)", 500,0,1 );
  //loaded += mHists->AddH1F(file,mH1F_LowPointMult,"H1F_LowPointMult", "Using best pi0 mass with point multiplicity (<4);Invariant Mass (GeV)", 500,0,1);

  //loaded += mHists->AddH1F(file,mH1F_GoodPointMult,"H1F_GoodPointMult","Point Multiplicity used for pi0 candidate pairs from photons with EPD cut;Point Multiplicity", 30,0,30);
  //loaded += mHists->AddH1F(file,mH1F_PointsEpd,"H1F_PointsEpd","Number of points with EPD nmip cut;Point Multiplicty",30,0,30);


  if( mH1F_InvMassEpdCuts[0]==0 ){ mH1F_InvMassEpdCuts[0] = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_InvMassEpdCuts[0],NEPDCUTS,"H1F_InvMassEpdCuts_AllTrig","Different EPD NMIP cuts all triggers",500,0,1);
  if( mH1F_InvMassEpdCuts[1]==0 ){ mH1F_InvMassEpdCuts[1] = new TObjArray(); }
  loaded += mHists->AddH1FArr(file,mH1F_InvMassEpdCuts[1],NEPDCUTS,"H1F_InvMassEpdCuts_EmTrig","Different EPD NMIP cuts EM triggers",500,0,1);
  //loaded += mHists->AddH1FArr(file,mH1F_InvMassZggCuts,8,"H1F_InvMassZggCuts",500,0,1);
  //loaded += mHists->AddH1FArr(file,mH1F_InvMassPtCuts,8,"H1F_InvMassPtCuts",500,0,1);
  loaded += mHists->AddH1F(file,mH1F_InvMassAllCuts,"H1F_AllInvMassCuts","Invariant Mass of two photons after all cuts applied", 500,0,1);
  loaded += mHists->AddH1F(file,mH1F_Pi0MultAllCuts,"H1F_Pi0MultAllCuts","Number of potential pi0s per event after all cuts", 20,0,20);
  //loaded += mHists->AddH2F(file,mH2F_AllCuts_Poi_yVx,"H2F_AllCuts_Poi_yVx","Point distribution y vs. x after all cuts", 400,-200,200, 300,-150,150);
  loaded += mHists->AddH2F(file,mH2F_AllCuts_Pi0_yVx,"H2F_AllCuts_Pi0_yVx","Pi0 distribution projected to FCS y vs. x after all cuts", 400,-200,200, 300,-150,150);

  loaded += mHists->AddH1F(file,mH1F_NFoundPhiBin,"H1F_NFoundPhiBin","Number of found phi bins",3,0,3);
  loaded += mHists->AddH1F(file,mH1F_AllCuts_xF,"H1F_AllCuts_xF","xF of pi0s after all cuts applied", 100,0,1);
  loaded += mHists->AddH1F(file,mH1F_AllCuts_Pi0En,"H1F_AllCuts_Pi0En","Energy of pi0s after all cuts applied", 100,0,100);
  //loaded += mHists->AddH1F(file,mH1F_AllCuts_Pi0Phi,"H1F_AllCuts_Pi0Phi","Phi of pi0s after all cuts applied", NPHIBIN,0,TMath::TwoPi());
  loaded += mHists->AddH2F(file,mH2F_AllCuts_Pi0_etaVphi,"H2F_AllCuts_Pi0_etaVphi","Eta vs. Phi distriubtion of pi0s", NPHIBIN,0,TMath::TwoPi(), 70,0,7);
  //loaded += mHists->AddH2F(file,mH2F_EpdNmip,"H2F_EpdNmip","EpdNmip;cluster;nmip", 2,0,2, 50,0,5);
			   
  //TString entitletext[NENERGYBIN] = { "En<=10", "10<En<=30", "30<En<=50", "50<En<=70", "70<En<=100", "En>100" };
  //TString phititletext[NPHIBIN] = { "0<=#phi<#frac{#pi}{4}", "#frac{#pi}{4}<=#phi<#frac{#pi}{2}", "#frac{#pi}{2}<=#phi<#frac{3#pi}{4}", "#frac{3#pi}{4}<=#phi<#pi" };
  double pi  = TMath::Pi();
  Double_t ebins[NENERGYBIN+1] = {0, 15, 20, 25, 30, 40, 55, 70, 100};

  loaded += mHists->AddH2F(file,mH2F_NPi0_enVphi[0][0],"H2F_NPi0_enVphi_blue_up","Number of Pi0s in a given energy and phi bin for blue beam spin up", NPHIBIN,-pi/2.0,3.0*pi/2.0, NENERGYBIN,ebins );
  loaded += mHists->AddH2F(file,mH2F_NPi0_enVphi[0][1],"H2F_NPi0_enVphi_blue_down","Number of Pi0s in a given energy and phi bin for blue beam spin down", NPHIBIN,-pi/2.0,3.0*pi/2.0, NENERGYBIN,ebins );
  loaded += mHists->AddH2F(file,mH2F_NPi0_enVphi[1][0],"H2F_NPi0_enVphi_yellow_up","Number of Pi0s in a given energy and phi bin for yellow beam spin up", NPHIBIN,-pi/2.0,3.0*pi/2.0, NENERGYBIN,ebins );
  loaded += mHists->AddH2F(file,mH2F_NPi0_enVphi[1][1],"H2F_NPi0_enVphi_yellow_down","Number of Pi0s in a given energy and phi bin for yellow beam spin down", NPHIBIN,-pi/2.0,3.0*pi/2.0, NENERGYBIN,ebins );

  //[January 29, 2025] > [How to create a TH3F with fixed and variable bin sizes](https://root-forum.cern.ch/t/how-to-make-th3-histograms-with-variable-bin-edges/38789)
  TAxis tmpphi(NPHIBIN,-pi/2.0,3.0*pi/2.0);
  Double_t phiedges[NPHIBIN+1] = {0};
  for( int i=0; i<25; ++i ){ phiedges[i] = tmpphi.GetBinUpEdge(i); }
  TAxis tmpmass(500,0,1);
  Double_t massedges[501] = {0};
  for( int i=0; i<501; ++i ){ massedges[i] = tmpmass.GetBinUpEdge(i); }
  loaded += mHists->AddH3F(file,mH3F_AllCutsInvMass_enVphi,"H3F_InvMass_envVphi","Invariant mass after all cuts by phi and energy binning", NPHIBIN,phiedges, NENERGYBIN,ebins, 500,massedges );
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

  UInt_t totalhists = this->LoadHists(0); //This is total of histograms loaded from a file not created. Don't use mFileOutput as you are not trying to load from #mFileOutput
  mHists->SetOwner(kTRUE);
  LOG_INFO << "StMuFcsPi0TreeMaker::Init() - Loaded " << totalhists << " histograms" << endm;

  return kStOk;
}

void StMuFcsPi0TreeMaker::LoadDataFromFile(TFile* file) //, TTree&* tree, FcsEventInfo&* evt,Int_t& ntrig, Int_t&* triggers,  TClonesArray&* pharr, TClonesArray&* pi0arr, TH1&* hist)
{
  if( file==0 || file->IsZombie() ){ std::cout << "LoadDataFromFile - ERROR:Unable to load from null file or zombie file" << std::endl; return; }
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
  if( mFcsTrigMap==0 ){ std::cout << "StMuFcsRun22QaMaker::LoadDataFromFile() - No Trigger Map found" << std::endl; }
  else{ if( mFcsTrigMap->sizeOfTriggers()<=0 ){ std::cout << "StMuFcsRun22QaMaker::LoadDataFromFile() - Trigger Map is empty" << std::endl; } }
  
  if( mHists==0 ){ mHists = new HistManager(); }
  UInt_t totalhists = this->LoadHists(file);
  //std::cout << "TotalHistsLoaded: "<<totalhists << std::endl;
  //std::cout << "DUMB CHECK "<<mFcsTrigMap << std::endl;
  //std::string name(mFcsTrigMap->nameFromId(45,22349011));
  //std::cout << name << std::endl;
  
  return;
}

Int_t StMuFcsPi0TreeMaker::InitRun(int runnumber) {
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
Int_t StMuFcsPi0TreeMaker::Make() {
  mMuDstMkr = (StMuDstMaker*)GetInputDS("MuDst");
  if( mMuDstMkr==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make - !MuDstMkr" <<endm; return kStErr; }
  mMuDst = mMuDstMkr->muDst();
  if( mMuDst==0 ){ LOG_ERROR << "StMuFcsPi0TreeMaker::Make - !MuDst" << endm; return kStErr; }
  mMuEvent = mMuDst->event();
  if( mMuEvent==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make - !MuEvent" <<endm; return kStErr; }
  mTrigData = mMuEvent->triggerData();
  if( mTrigData==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make - !TrigData" <<endm; return kStErr; }
  mRunInfo = &(mMuEvent->runInfo());
  if( mRunInfo==0 ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::Make - !RunInfo" <<endm; return kStErr; }

  mH1F_Entries->Fill(1); //This is just counting valid make calls
  
  //Filter with Trigger Information first
  bool validtrigfound = false;
  bool emtrigfound = false;
  if( mIgnoreTrig ){ validtrigfound = true; } //If ignoring triggers then set validtrigfound to true so event is not skipped
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
	    validtrigfound=true;
	    mH1F_Triggers->Fill( thistrig.c_str(), 1);
	  }
	}
      }
      if( thistrig=="fcsEM0" || thistrig=="fcsEM1" || thistrig=="fcsEM2" || thistrig=="fcsEM3"
	  || thistrig=="fcsEM0_tpc" || thistrig=="fcsEM1_tpc" || thistrig=="fcsEM2_tpc" || thistrig=="fcsEM3_tpc" )
	{
	  emtrigfound = true;
	}
    }
    else{ mNTrig = ntrig; }
  }
  if( !validtrigfound ){ mNTrig=0; return kStSkip; } //Reset trigger array size before going to next event

  //Get EPD collection and/or hits
  mMuEpdHits = 0;
  mEpdColl = 0;
  mMuEpdHits = mMuDst->epdHits();
  if( mMuEpdHits!=0 ){ if( mMuEpdHits->GetEntriesFast()==0 ){mMuEpdHits=0;} }//If mMuEpdHits is not zero but has no hits set it to zero so rest of code processes from StEpdHitMaker
  if( mMuEpdHits==0 ){ LOG_INFO << "StMuFcsPi0TreeMaker::Make - No MuEPD hits" << endm;
    mEpdHitMkr = (StEpdHitMaker*)GetMaker("epdHit");
    if( mEpdHitMkr==0 ){ LOG_WARN << "StMuFcsPi0TreeMaker::Make - No StEpdHitMaker(\"epdHit\")" << endm; }
    else{ mEpdColl = mEpdHitMkr->GetEpdCollection(); }
    if( mEpdColl==0 ){ LOG_WARN << "StMuEpdRun22QaMaker::FillFcsInfo - No Epd hit information found" << endm; mEpdHitMkr=0; }//Set the hit maker back to zero so it can be used as a check that the epd collection doesn't exist
  }

  //Fcs Collection
  mMuFcsColl = mMuDst->muFcsCollection();
  if (!mMuFcsColl) { LOG_ERROR << "StMuFcsPi0TreeMaker::Make did not find MuFcsCollection" << endm; return kStErr; }

  //FcsEventInfo* evtinfo = (FcsEventInfo*) mEvtInfo->ConstructedAt(0);  
  mEvtInfo->mRunTime         = mMuEvent->eventInfo().time();
  mEvtInfo->mRunNum          = mMuEvent->runNumber();
  mEvtInfo->mFill            = mRunInfo->beamFillNumber(StBeamDirection::east);//Yellow beam
  mEvtInfo->mEvent           = mMuEvent->eventId();
  mEvtInfo->mBx48Id          = mTrigData->bunchId48Bit();
  mEvtInfo->mBx7Id           = mTrigData->bunchId7Bit();
  mEvtInfo->mTofMultiplicity = mTrigData->tofMultiplicity();

  //std::cout << "|runtime:"<<mEvtInfo->mRunTime << "|runnum:"<<mEvtInfo->mRunNum << "|event:"<<mEvtInfo->mEvent << std::endl;
  
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

  short foundvertex = 0;     //Bit vector for knowing which vertex was used 0 means no vertex, 1=Vpd,2=Epd,4=Bbc
  Double_t usevertex = -999.0; //Vertex to use
  if( mEvtInfo->mVpdVz > -998 )     { foundvertex = 1; usevertex = mEvtInfo->mVpdVz; }
  else if( mEvtInfo->mEpdVz > -998 ){ foundvertex = 2; usevertex = mEvtInfo->mEpdVz; }
  else if( mEvtInfo->mBbcVz > -998 ){ foundvertex = 4; usevertex = mEvtInfo->mBbcVz; }
  else{ foundvertex = 0; usevertex = -999; } //Redundant but don't want a dangling "else if" statement
  mEvtInfo->mFoundVertex = foundvertex;
  mH2F_foundVvertex->Fill(usevertex,foundvertex);

  //TClonesArray* hits = mMuFcsColl->getHitArray();
  //if( hits==0 ){ LOG_INFO << "StMuFcsRun22QaMaker::FillFcsInfo - No FCS hits" << endm; }
  TClonesArray* clusters = mMuFcsColl->getClusterArray();
  if( clusters==0 ){ LOG_INFO << "StMuFcsRun22QaMaker::FillFcsInfo - No FCS clusters" << endm; }
  TClonesArray* points = mMuFcsColl->getPointArray();
  if( points==0 ){ LOG_INFO << "StMuFcsRun22QaMaker::FillFcsInfo - No FCS points" << endm; }
  
  //std::cout << "|hits:"<<hits << "|clusters:"<<clusters << "|points:"<<points << std::endl;
  Int_t ncandidates = 0;
  if( clusters!=0 ){
    for( UInt_t idet=0; idet<=kFcsEcalSouthDetId; ++idet ){
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

	if( foundvertex > 0 ){
	  StLorentzVectorD iclu_p_withz = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, usevertex );
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
      unsigned int np = mMuFcsColl->numberOfPoints(idet);
      unsigned int ipoint=mMuFcsColl->indexOfFirstPoint(idet);
      np += ipoint;
      for( ; ipoint<np; ++ipoint ){
	StMuFcsPoint* point = (StMuFcsPoint*)points->At(ipoint);
	float ipoi_x = point->x();
	float ipoi_y = point->y();
	float ipoi_energy = point->energy();
	if( ipoi_energy<mEnCut ){ continue; }

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

	if( foundvertex > 0 ){ 	StLorentzVectorD ipoi_p_withz = mFcsDb->getLorentzVector( ipoi_pos, ipoi_energy, usevertex );
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

  //Check photon candidates if they have any hits in the EPD. Use a separate loop so that this information could be used in the pi0 checking loop if needed. In future may also want to check against FCS preshower (EPD) hits
  for( Int_t iph = 0; iph<mPhArr->GetEntriesFast(); ++iph ){
    FcsPhotonCandidate* ph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(iph);
    std::vector<Double_t> epdproj(ProjectToEpd(ph->mX,ph->mY,ph->mZ,usevertex));
    if( !(ph->mFromCluster) ){
      mH2F_EpdProjHitMap->Fill( epdproj.at(0),epdproj.at(1) );
      if( -150<=usevertex && usevertex<=150 ){ mH2F_EpdProjHitMap_Vcut->Fill(epdproj.at(0),epdproj.at(1)); }
    }
    //std::cout << " + |phx:"<<ph->mX << "|phy:"<<ph->mY << "|phz:"<<ph->mZ << "|v:"<<usevertex << std::endl;
    //std::cout << " + |epdx:"<<epdproj.at(0) << "|epdy:"<<epdproj.at(1) << "|epdz:"<<epdproj.at(2) << std::endl;

    //loop over all west epd tiles so that even if no hit recorded can use as a veto
    for(int i_pp=1; i_pp<=12; ++i_pp){     //Supersector runs [1,12]
      for( int i_tt=1; i_tt<=31; ++i_tt ){ //Tile number [1,31]
	if( mEpdGeo->IsInTile(i_pp,i_tt, 1, epdproj.at(0),epdproj.at(1)) ){ //Only care about west EPD tiles; hence the '1'
	  ph->mEpdHitNmip = 0;
	}
      }
    }
    //loop over all hits and if an nmip value exists set for the point
    unsigned int nepdhits = 0;
    StSPtrVecEpdHit* epdhits = 0;
    if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
    else if( mEpdColl!=0 ){
      epdhits = &(mEpdColl->epdHits());
      nepdhits = epdhits->size();
    }
    else{ LOG_ERROR << "StMuEpdRun22QaMaker::FillEpdinfo() - If you see this error then there is a bug that is setting EPD hits improperly" << endm; return kStErr; }
    StMuEpdHit* muepdhit = 0;
    StEpdHit* epdhit = 0;
    for(unsigned int i=0; i<nepdhits; ++i ){
      if( mMuEpdHits!=0 ){ muepdhit = (StMuEpdHit*)mMuEpdHits->UncheckedAt(i); } //To match similar in StMuDstMaker->epdHit(int i)
      else if( epdhits!=0 ){ epdhit = (StEpdHit*)((*epdhits)[i]); }
      else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A VERY SERIOUS BUG IN THE CODE" << endm; return kStErr; } 
      //std::cout << "|i:"<<i << "|muepdhit:"<<muepdhit << "|epdhit:"<<epdhit << std::endl;
      int ew    = muepdhit!=0 ? muepdhit->side()    : epdhit->side();      //east=-1, west=1
      if( ew==-1 ){ continue; }
      int epdpp = muepdhit!=0 ? muepdhit->position(): epdhit->position();  //Supersector runs [1,12]
      int epdtt = muepdhit!=0 ? muepdhit->tile()    : epdhit->tile();      //Tile number [1,31]
      if( mEpdGeo->IsInTile(epdpp,epdtt,ew, epdproj.at(0),epdproj.at(1)) ){
	//int adc = muepdhit!=0 ? muepdhit->adc() : epdhit->adc();
	float nmip = muepdhit!=0 ? muepdhit->nMIP(): epdhit->nMIP();         //The ADC value of the hit divided by the MIP peak position; e.g. if nmip==1 then adc value sits at the MIP peak
	ph->mEpdHitNmip = nmip;
	//std::cout << "|epdpp:"<<epdpp <<"|epdtt:"<<epdtt <<"|nmip:"<<nmip << std::endl;
      }
    }
    mH2F_EpdNmip->Fill(ph->mFromCluster,ph->mEpdHitNmip);
  }


  //std::cout << "|ncandidates:"<<ncandidates <<"|clustersize:"<<clustersize <<"|Size:"<<mPhArr->GetEntriesFast() << std::endl;
  Int_t npoints = ncandidates - clustersize; //Don't need to add 1 since including clustersize but not ncandidates
  mH1F_PointMult->Fill(npoints);
  mPhArr->Sort(); //Since this is properly sorted with clusters showing up first clustersize is unchanged. Also sorts by energy

  //std::cout << "|clustersize:"<<clustersize << "|ncandidates:"<<ncandidates << "|npoints:"<<npoints << std::endl;
  
  Int_t npi0candidate = 0;
  //Filling cluster pi0s and cluster photon/elecron epd nmip cut
  std::vector<Int_t> goodclusphotonsidx;    //Here I really mean photon candidates below epd cut threshold
  std::vector<Int_t> goodcluselectronsidx;  //Here I really mean photon candidates above epd cut threshold
  for( Int_t ic = 0; ic<clustersize; ++ic ){
    FcsPhotonCandidate* iclus = (FcsPhotonCandidate*) mPhArr->UncheckedAt(ic);
    if( !(iclus->mFromCluster) ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
    //std::cout << "|ic:"<<ic << std::endl;
    //std::cout << "  + ";
    //iclus->Print();
    if( ic==(clustersize-1) ){ continue; }

    if( iclus->mEpdHitNmip>-0.1){ //Only include candidates who have their nmip value set
      if( iclus->mEpdHitNmip<mEpdNmipCut ){ goodclusphotonsidx.emplace_back(ic); }
      else{ goodcluselectronsidx.emplace_back(ic); }
    }    
    for( Int_t jc=ic+1; jc<clustersize; jc++ ){
      FcsPhotonCandidate* jclus = (FcsPhotonCandidate*) mPhArr->UncheckedAt(jc);
      if( !(jclus->mFromCluster) ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
      TLorentzVector pi0Vert_LV = iclus->lvVert() + jclus->lvVert();
      if( ic==0 ){ //Since we have a sorted photon array highest two energies are the first two entries
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
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

  //Making pi0s with cluster cut photon candidates
  for( UInt_t iidx=0; iidx<goodclusphotonsidx.size(); ++iidx ){
    if( iidx==(goodclusphotonsidx.size()-1) ){ continue; }
    FcsPhotonCandidate* iepdph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodclusphotonsidx.at(iidx));
    for( UInt_t jidx=iidx+1; jidx<goodclusphotonsidx.size(); ++jidx ){
      FcsPhotonCandidate* jepdph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodclusphotonsidx.at(jidx));
      TLorentzVector pi0Vert_LV = iepdph->lvVert() + jepdph->lvVert();
      if( iidx==0 ){ //Since we sorted array by index first two indices represent highest energy pairs with epd cut
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
	pi0c->mFromCluster = true;
	pi0c->mFromPh = -1;
	pi0c->mPhoton1Idx = goodclusphotonsidx.at(iidx);
	pi0c->mPhoton2Idx = goodclusphotonsidx.at(jidx);
	
	pi0c->mPx = pi0Vert_LV.Px();
	pi0c->mPy = pi0Vert_LV.Py();
	pi0c->mPz = pi0Vert_LV.Pz();
	pi0c->mEn = pi0Vert_LV.E();
	
	pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
	pi0c->mDgg     = FcsPi0Candidate::dgg(*iepdph,*jepdph);
	pi0c->mZgg     = FcsPi0Candidate::zgg(*iepdph,*jepdph);
	pi0c->mAlpha   = FcsPi0Candidate::alpha(*iepdph,*jepdph);
	pi0c->mInvMass = pi0Vert_LV.Mag();
	break;
      }
      /*
      else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_EpdPhAllPoints->Fill(pi0Vert_LV.Mag());
	}*/
    }
  }

  //Making pi0s with cluster cut electron candidates
  for( UInt_t iidx=0; iidx<goodcluselectronsidx.size(); ++iidx ){
    if( iidx==(goodcluselectronsidx.size()-1) ){ continue; }
    FcsPhotonCandidate* iepdch = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodcluselectronsidx.at(iidx));
    for( UInt_t jidx=iidx+1; jidx<goodcluselectronsidx.size(); ++jidx ){
      FcsPhotonCandidate* jepdch = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodcluselectronsidx.at(jidx));
      if( iidx==0 ){ //Since we sorted array by index first two indices represent highest energy pairs with epd cut
	TLorentzVector pi0Vert_LV = iepdch->lvVert() + jepdch->lvVert();
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
	pi0c->mFromCluster = true;
	pi0c->mFromPh = 1;
	pi0c->mPhoton1Idx = goodcluselectronsidx.at(iidx);
	pi0c->mPhoton2Idx = goodcluselectronsidx.at(jidx);
	
	pi0c->mPx = pi0Vert_LV.Px();
	pi0c->mPy = pi0Vert_LV.Py();
	pi0c->mPz = pi0Vert_LV.Pz();
	pi0c->mEn = pi0Vert_LV.E();
	
	pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
	pi0c->mDgg     = FcsPi0Candidate::dgg(*iepdch,*jepdch);
	pi0c->mZgg     = FcsPi0Candidate::zgg(*iepdch,*jepdch);
	pi0c->mAlpha   = FcsPi0Candidate::alpha(*iepdch,*jepdch);
	pi0c->mInvMass = pi0Vert_LV.Mag();
	break;
      }
      /*
      else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_EpdChAllPoints->Fill(pi0Vert_LV.Mag());
	}*/
    }
  }


  //Filling point pi0s and cluster photon/elecron epd nmip cut
  std::vector<Int_t> goodpointphotonsidx;      //Here I really mean photon candidates below epd cut threshold
  std::vector<Int_t> goodpointelectronsidx;    //Here I really mean photon candidates above epd cut threshold
  //std::cout << "===== EventId:"<< mEvtInfo->mEvent <<" =====" << std::endl;
  FcsPhotonCandidate* firstphotoncut[NEPDCUTS] = {0};
  FcsPhotonCandidate* secondphotoncut[NEPDCUTS] = {0};
  Int_t n_noepdproj = 0;
  Int_t n_noepdproj_vcut = 0;
  for( Int_t ip = clustersize; ip<ncandidates; ++ip ){
    FcsPhotonCandidate* ipoi = (FcsPhotonCandidate*) mPhArr->UncheckedAt(ip);
    if( ipoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
    //std::cout << "|ip:"<<ip << std::endl;
    //std::cout << "  + ";
    //ipoi->Print();

    if( ipoi->mEpdHitNmip>-0.1){ //Only include candidates who have their nmip value set
      if( ipoi->mEpdHitNmip<mEpdNmipCut ){ goodpointphotonsidx.emplace_back(ip); }
      else{ goodpointelectronsidx.emplace_back(ip); }
      for( short i=0; i<NEPDCUTS; ++i ){
	if( ipoi->mEpdHitNmip < 0.2+0.1*static_cast<double>(i) ){
	  if( firstphotoncut[i]==0 ){ firstphotoncut[i]=ipoi; }
	  else{ if( secondphotoncut[i]==0 ){ secondphotoncut[i]=ipoi; } }
	}
      }
    }
    else{
      ++n_noepdproj;
      if( -150<=usevertex && usevertex<=150 ){ ++n_noepdproj_vcut; }
    }
    
    if( ip==(mPhArr->GetEntriesFast()-1) ){ continue; }
    for( Int_t jp=ip+1; jp<ncandidates; ++jp ){
      FcsPhotonCandidate* jpoi = (FcsPhotonCandidate*) mPhArr->UncheckedAt(jp);
      if( jpoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
      TLorentzVector pi0Vert_LV = ipoi->lvVert() + jpoi->lvVert();
      if( ip==clustersize ){ //Since we have a sorted photon array highest two energies are the first two entries
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
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
	mH1F_BestPi0Mass->Fill(pi0c->mInvMass);
	//mH2F_BestPi0HeatMap->Fill();
	mH1F_BestPi0Zgg->Fill(pi0c->mZgg);
	mH1F_BestPi0Phi->Fill(pi0c->phi());
	mH1F_BestPi0Eta->Fill(pi0c->mEta);
	mH1F_BestPi0En->Fill(pi0c->mEn);
	mH1F_BestPi0Pt->Fill(pi0c->pt());
	mH1F_AllPointPairMass->Fill(pi0Vert_LV.Mag());
	mH2F_Energy_ph1Vph2->Fill(ipoi->mEn,jpoi->mEn);
      }
      else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_AllPointPairMass->Fill(pi0Vert_LV.Mag());
      }
    }
  }
  mH1F_NBadEpdProj->Fill(n_noepdproj);
  mH1F_NBadEpdProjVcut->Fill(n_noepdproj_vcut);

  //Handle the epd cuts
  for( short i=0; i<NEPDCUTS; ++i ){
    if( firstphotoncut[i]!=0 && secondphotoncut[i]!=0 ){
      TLorentzVector pi0Vert_LV = firstphotoncut[i]->lvVert() + secondphotoncut[i]->lvVert();
      if( emtrigfound ){ ((TH1*)mH1F_InvMassEpdCuts[1]->UncheckedAt(i))->Fill(pi0Vert_LV.Mag()); }
      ((TH1*)mH1F_InvMassEpdCuts[0]->UncheckedAt(i))->Fill(pi0Vert_LV.Mag());
    }
  }

  //Making pi0s with point cut photon candidates
  mH1F_EpdPhPointMult->Fill(goodpointphotonsidx.size());
  for( UInt_t iidx=0; iidx<goodpointphotonsidx.size(); ++iidx ){
    if( iidx==(goodpointphotonsidx.size()-1) ){ continue; }
    FcsPhotonCandidate* iepdph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodpointphotonsidx.at(iidx));
    for( UInt_t jidx=iidx+1; jidx<goodpointphotonsidx.size(); ++jidx ){
      FcsPhotonCandidate* jepdph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodpointphotonsidx.at(jidx));
      TLorentzVector pi0Vert_LV = iepdph->lvVert() + jepdph->lvVert();
      if( iidx==0 ){ //Since we sorted array by index first two indices represent highest energy pairs with epd cut
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
	pi0c->mFromCluster = false;
	pi0c->mFromPh = -1;
	pi0c->mPhoton1Idx = goodpointphotonsidx.at(iidx);
	pi0c->mPhoton2Idx = goodpointphotonsidx.at(jidx);
	
	pi0c->mPx = pi0Vert_LV.Px();
	pi0c->mPy = pi0Vert_LV.Py();
	pi0c->mPz = pi0Vert_LV.Pz();
	pi0c->mEn = pi0Vert_LV.E();
	
	pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
	pi0c->mDgg     = FcsPi0Candidate::dgg(*iepdph,*jepdph);
	pi0c->mZgg     = FcsPi0Candidate::zgg(*iepdph,*jepdph);
	pi0c->mAlpha   = FcsPi0Candidate::alpha(*iepdph,*jepdph);
	pi0c->mInvMass = pi0Vert_LV.Mag();
	mH1F_EpdPhInvMass->Fill(pi0c->mInvMass);
	//mH2F_BestPi0HeatMap->Fill();
	mH1F_EpdPhZgg->Fill(pi0c->mZgg);
	mH1F_EpdPhPhi->Fill(pi0c->phi());
	mH1F_EpdPhEta->Fill(pi0c->mEta);
	mH1F_EpdPhEn->Fill(pi0c->mEn);
	mH1F_EpdPhPt->Fill(pi0c->pt());
	mH1F_EpdPhAllPoints->Fill(pi0Vert_LV.Mag());
      }
      else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_EpdPhAllPoints->Fill(pi0Vert_LV.Mag());
      }
    }
  }


  //Making pi0s with point cut electron candidates
  mH1F_EpdChPointMult->Fill(goodpointelectronsidx.size());
  for( UInt_t iidx=0; iidx<goodpointelectronsidx.size(); ++iidx ){
    if( iidx==(goodpointelectronsidx.size()-1) ){ continue; }
    FcsPhotonCandidate* iepdch = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodpointelectronsidx.at(iidx));
    for( UInt_t jidx=iidx+1; jidx<goodpointelectronsidx.size(); ++jidx ){
      FcsPhotonCandidate* jepdch = (FcsPhotonCandidate*) mPhArr->UncheckedAt(goodpointelectronsidx.at(jidx));
      TLorentzVector pi0Vert_LV = iepdch->lvVert() + jepdch->lvVert();
      if( iidx==0 ){ //Since we sorted array by index first two indices represent highest energy pairs with epd cut
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
	pi0c->mFromCluster = false;
	pi0c->mFromPh = 1;
	pi0c->mPhoton1Idx = goodpointelectronsidx.at(iidx);
	pi0c->mPhoton2Idx = goodpointelectronsidx.at(jidx);
	
	pi0c->mPx = pi0Vert_LV.Px();
	pi0c->mPy = pi0Vert_LV.Py();
	pi0c->mPz = pi0Vert_LV.Pz();
	pi0c->mEn = pi0Vert_LV.E();
	
	pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
	pi0c->mDgg     = FcsPi0Candidate::dgg(*iepdch,*jepdch);
	pi0c->mZgg     = FcsPi0Candidate::zgg(*iepdch,*jepdch);
	pi0c->mAlpha   = FcsPi0Candidate::alpha(*iepdch,*jepdch);
	pi0c->mInvMass = pi0Vert_LV.Mag();
	mH1F_EpdChInvMass->Fill(pi0c->mInvMass);
	//mH2F_BestPi0HeatMap->Fill();
	mH1F_EpdChZgg->Fill(pi0c->mZgg);
	mH1F_EpdChPhi->Fill(pi0c->phi());
	mH1F_EpdChEta->Fill(pi0c->mEta);
	mH1F_EpdChEn->Fill(pi0c->mEn);
	mH1F_EpdChPt->Fill(pi0c->pt());
	mH1F_EpdChAllPoints->Fill(pi0Vert_LV.Mag());
      }
      else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_EpdChAllPoints->Fill(pi0Vert_LV.Mag());
      }
    }
  }

  //XF = pi0Vert_LV.Pz()/255;
  //Double_t PiOverN = TMath::Pi()/static_cast<Double_t>(NPHIBIN);
  //Double_t PhiBin[NPHIBIN*2+1];  //Last element is to make sure it comes full circle and checking if phi is between 7Pi/4 to 8Pi/4 is easier than 7Pi/4 and 0
  //for( short i=0; i<=2*NPHIBIN; ++i ){ PhiBin[i] = static_cast<Double_t>(i)*PiOverN; }
  //{0.0, PiOver4, 2.0*PiOver4, 3.0*PiOver4, 4.0*PiOver4, 5.0*PiOver4, 6.0*PiOver4, 7.0*PiOver4, 8.0*PiOver4};
  //Int_t NPi0UpAtPhi[NENERGYBIN][NPHIBIN]         = {0};
  //Int_t NPi0DownAtPhi[NENERGYBIN][NPHIBIN]       = {0};
  //Int_t NPi0UpAtPhiPlusPi[NENERGYBIN][NPHIBIN]   = {0};
  //Int_t NPi0DownAtPhiPlusPi[NENERGYBIN][NPHIBIN] = {0};
  Int_t ngoodpi0s = 0;
  //std::cout << "NPi0s:"<<mPi0Arr->GetEntriesFast() << std::endl;
  for( int i=0; i<mPi0Arr->GetEntriesFast(); ++i ){
    FcsPi0Candidate* pi0 = (FcsPi0Candidate*)mPi0Arr->At(i);
    if( pi0==0 ){ continue; }
    //pi0->Print();
    Double_t pi0en = pi0->mEn;
    //if( !emtrigfound )     { continue; }
    if( ! (-150<=usevertex && usevertex<=150) ){ /*std::cout << " StMuFcsPi0TreeMaker::Make() - Failed vertex cut:"<<usevertex << std::endl;*/ continue; }
    if( pi0->mFromCluster ){ /*std::cout << "StMuFcsPi0TreeMaker::Make() - Not a point - "<<pi0->mFromCluster<< std::endl;*/ continue; }
    if( pi0->mZgg>0.7 )    { /*std::cout << "StMuFcsPi0TreeMaker::Make() - Failed Zgg:"<< pi0->mZgg << std::endl;*/ continue; }
    if( pi0->mFromPh>=0 )  { /*std::cout << "StMuFcsPi0TreeMaker::Make() - Failed photon cut - "<<pi0->mFromPh<< std::endl;*/ continue; }  //Ensure pi0 is coming from photon cut
    //Add pt cut based on trigger
    bool exceedtrigpt = false;
    Float_t trigptthr = -1;
    std::string trigname = "";
    if( !mIgnoreTrig ){
      if( !emtrigfound ){ continue; }
      if( mFcsTrigMap!=0 ){
	for( Int_t itrig=0; itrig<mNTrig; ++itrig ){
	  Float_t ptthres = mFcsTrigMap->GetPtThr(mTriggers[itrig]);
	  std::string thistrig = mFcsTrigMap->nameFromId(mTriggers[itrig],mMuEvent->runNumber());
	  /*
	    if( thistrig=="fcsEM0" || thistrig=="fcsEM1" || thistrig=="fcsEM2" || thistrig=="fcsEM3"
	    || thistrig=="fcsEM0_tpc" || thistrig=="fcsEM1_tpc" || thistrig=="fcsEM2_tpc" || thistrig=="fcsEM3_tpc" )
	    {
	    }
	  */
	  if( pi0->pt()>=ptthres ){ exceedtrigpt=true; trigptthr=ptthres; trigname=thistrig; }
	}
      }
    }
    else{ exceedtrigpt = true; }
    if( !exceedtrigpt ){ continue; }
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
    TLorentzVector pi0_lv = pi0->lv();
    double pi0_momentum[3] = {pi0->mPx,pi0->mPy,pi0->mPz};
    double pi0_vertex[3] = {0,0,usevertex};
    int det = 0; //North side if negative px
    //South side for positive px, set px==0 to south side since fcs planes would intersect in that case
    if( pi0_momentum[2]>=0 && pi0_momentum[0]>=0 ){ det=1; }
    if( pi0_momentum[2]<0  && pi0_momentum[0]<0  ){ det=1; }
    StThreeVectorD pi0_xyz = mFcsDb->projectLine(det,pi0_momentum,pi0_vertex,0); //Project to front face of FCS
    mH2F_AllCuts_Pi0_yVx->Fill(pi0_xyz.x(),pi0_xyz.y());

    Double_t phi = pi0_lv.Phi(); //Range of this phi is -pi to pi
    Double_t mpi = -TMath::Pi();
    if( mpi<=phi && phi<mpi/2.0 ){ phi += TMath::TwoPi(); } //Since my binning goes from -pi/2 to 3pi/2 need to add 2pi to angles in the region from [-pi,-pi/2)
    //@[Jan 27, 2025] > (Keep energy binnng and check energy plot for binning) Make a phi histogram one for blue (up & down) and yellow (up & down) from -Pi/2 to 3/2 Pi with whatever binning. bin 1 which is bottom most on left, with nbins/2, TH1F* mhphi[2][2] = {0}; //[blue,yellow] [up,down], Move to Ana code [mhasym[2] //[blue,yellow beam]. mhphi[0][0]->bin(1)*mhphi[0][1]->bin(1+nbin/2)]
    mH1F_AllCuts_xF->Fill( pi0->mPz/255 );
    mH1F_AllCuts_Pi0En->Fill( pi0en );
    mH2F_AllCuts_Pi0_etaVphi->Fill( phi,pi0->eta() );
    //short nfoundphibin = 0;
    if( 0.1<=pi0->mass() && pi0->mass()<=0.2){
      if( mEvtInfo->BlueSpin()==1    ){ mH2F_NPi0_enVphi[0][0]->Fill(phi,pi0en); }
      if( mEvtInfo->BlueSpin()==-1   ){ mH2F_NPi0_enVphi[0][1]->Fill(phi,pi0en); }
      if( mEvtInfo->YellowSpin()==1  ){ mH2F_NPi0_enVphi[1][0]->Fill(phi,pi0en); }
      if( mEvtInfo->YellowSpin()==-1 ){ mH2F_NPi0_enVphi[1][1]->Fill(phi,pi0en); }
    }
    /*
    for( int iphi=0; iphi<NPHIBIN; ++iphi ){
      if( mEvtInfo->BlueSpin()==1  && PhiBin[iphi]<=phi   && phi<PhiBin[iphi+1] ){ ++NPi0UpAtPhi[enbin][iphi]; phibin=iphi; ++nfoundphibin; } //Spin up and on right side is N^{up}(phi)
      if( mEvtInfo->BlueSpin()==1  && PhiBin[iphi+4]<=phi && phi<PhiBin[iphi+5] ){ ++NPi0UpAtPhiPlusPi[enbin][iphi]; phibin=iphi; ++nfoundphibin; } //Spin up and on left side is N^{up}(phi+pi)
      if( mEvtInfo->BlueSpin()==-1 && PhiBin[iphi]<=phi   && phi<PhiBin[iphi+1] ){ ++NPi0DownAtPhi[enbin][iphi]; phibin=iphi; ++nfoundphibin; } //Spin down and on right side is N^{down}(phi)
      if( mEvtInfo->BlueSpin()==-1 && PhiBin[iphi+4]<=phi && phi<PhiBin[iphi+5] ){ ++NPi0DownAtPhiPlusPi[enbin][iphi]; phibin=iphi; ++nfoundphibin; } //Spin down and on left side is N^{down}(phi+pi)
    }
    if( phibin<0 ){ std::cout << "Couldn't find a valid phi bin" << std::endl; break; }
    mH1F_NFoundPhiBin->Fill(nfoundphibin);
    */
    //mH1F_InvMassAllCutsByEnByPhi[enbin][phibin]->Fill(pi0->mass());
    ((TH3*)mH3F_AllCutsInvMass_enVphi)->Fill(phi,pi0en,(Double_t)pi0->mass());
    mH1F_InvMassAllCuts->Fill(pi0->mass());
  }
  mH1F_Pi0MultAllCuts->Fill(ngoodpi0s);
  //std::cout << "NGoodPi0s:"<<ngoodpi0s << std::endl;

  /*
  for( short ebin=0; ebin<NENERGYBIN; ++ebin ){
    for( short phibin=0; phibin<NPHIBIN; ++phibin ){
      mH1F_NPi0ByEnByPhi[ebin][phibin]->SetBinContent(1,NPi0UpAtPhi[ebin][phibin]);
      mH1F_NPi0ByEnByPhi[ebin][phibin]->SetBinContent(2,NPi0UpAtPhiPlusPi[ebin][phibin]);
      mH1F_NPi0ByEnByPhi[ebin][phibin]->SetBinContent(3,NPi0DownAtPhi[ebin][phibin]);
      mH1F_NPi0ByEnByPhi[ebin][phibin]->SetBinContent(4,NPi0DownAtPhiPlusPi[ebin][phibin]);
    }
  }
  */
  if( mPi0Tree!=0 ){ mPi0Tree->Fill(); }

  mEvtInfo->Clear();
  mNTrig = 0; //Since ROOT only writes up to the size of mNTrig then only need to reset this back to zero and next loop will overwrite array as neccessary
  mPhArr->Clear();
  mPi0Arr->Clear();
  
  return kStOk;
}

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
  mH1F_Entries->Draw("hist e");
  canv->cd(2)->SetLogy();
  mH1F_Triggers->Draw("hist e");
  canv->cd(3)->SetLogz();
  mH2F_foundVvertex->Draw("colz");
  canv->cd(4);
  mH1F_RndmSpin->Draw("hist e");

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
  mH1F_PointEnergy->Draw("hist e");
  canv->cd(7)->SetLogz();
  mH2F_Energy_ph1Vph2->Draw("colz");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintBestPi0(TCanvas* canv,  const char* savename)  const
{
  canv->Clear();
  
  canv->Divide(3,3);
  canv->cd(1);
  mH1F_PointMult->Draw("hist e");
  canv->cd(2);
  mH1F_BestPi0Zgg->Draw("hist e");
  canv->cd(3);
  mH1F_BestPi0Phi->Draw("hist e");
  canv->cd(4);
  mH1F_BestPi0Eta->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_BestPi0En->Draw("hist e");
  canv->cd(6)->SetLogy();
  mH1F_BestPi0Pt->Draw("hist e");
  canv->cd(7);
  mH1F_BestPi0Mass->Draw("hist e");
  canv->cd(8);
  mH1F_AllPointPairMass->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdPhPi0(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,3);
  canv->cd(1);
  mH1F_EpdPhPointMult->Draw("hist e");
  canv->cd(2);
  mH1F_EpdPhZgg->Draw("hist e");
  canv->cd(3);
  mH1F_EpdPhPhi->Draw("hist e");
  canv->cd(4);
  mH1F_EpdPhEta->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_EpdPhEn->Draw("hist e");
  canv->cd(6)->SetLogy();
  mH1F_EpdPhPt->Draw("hist e");
  canv->cd(7);
  mH1F_EpdPhInvMass->Draw("hist e");
  canv->cd(8);
  mH1F_EpdPhAllPoints->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintEpdChPi0(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,3);
  canv->cd(1);
  mH1F_EpdChPointMult->Draw("hist e");
  canv->cd(2);
  mH1F_EpdChZgg->Draw("hist e");
  canv->cd(3);
  mH1F_EpdChPhi->Draw("hist e");
  canv->cd(4);
  mH1F_EpdChEta->Draw("hist e");
  canv->cd(5)->SetLogy();
  mH1F_EpdChEn->Draw("hist e");
  canv->cd(6)->SetLogy();
  mH1F_EpdChPt->Draw("hist e");
  canv->cd(7);
  mH1F_EpdChInvMass->Draw("hist e");
  canv->cd(8);
  mH1F_EpdChAllPoints->Draw("hist e");

  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintPi0Overlap(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  
  canv->Divide(3,3);

  //canv->cd(1);
  canv->cd(1)->SetLogy();
  TLegend* legpad1 = new TLegend(0.5,0.5,0.93,0.93,"","nbNDC");
  TH1* h1pmult = mH1F_PointMult->DrawNormalized("hist e");
  h1pmult->SetStats(0);
  h1pmult->SetLineColor(kBlack);
  AddHistStatsOneline(legpad1,h1pmult,"NoCut");
  TH1* h1phmult = mH1F_EpdPhPointMult->DrawNormalized("hist e same");
  h1phmult->SetLineColor(kBlue);
  AddHistStatsOneline(legpad1,h1phmult,"EpdNmip<0.7");
  TH1* h1chmult = mH1F_EpdChPointMult->DrawNormalized("hist e same");
  AddHistStatsOneline(legpad1,h1chmult,"EpdNmip>=0.7");
  h1chmult->SetLineColor(kGreen+2);
  legpad1->Draw();
  
  //canv->cd(2);
  canv->cd(2);//->SetLogy();
  TH1* h1pzgg = mH1F_BestPi0Zgg->DrawNormalized("hist e");
  h1pzgg->SetLineColor(kBlack);
  TH1* h1phzgg = mH1F_EpdPhZgg->DrawNormalized("hist e same");
  h1phzgg->SetLineColor(kBlue);
  TH1* h1chzgg  = mH1F_EpdChZgg->DrawNormalized("hist e same");
  h1chzgg->SetLineColor(kGreen+2);
  
  canv->cd(3);
  //canv->cd(3)->SetLogy();
  TH1* h1pphi = mH1F_BestPi0Phi->DrawNormalized("hist e");
  h1pphi->SetLineColor(kBlack);
  TH1* h1phphi = mH1F_EpdPhPhi->DrawNormalized("hist e same");
  h1phphi->SetLineColor(kBlue);
  TH1* h1chphi = mH1F_EpdChPhi->DrawNormalized("hist e same");
  h1chphi->SetLineColor(kGreen+2);
  
  canv->cd(4);
  //canv->cd(4)->SetLogy();
  TH1* h1peta = mH1F_BestPi0Eta->DrawNormalized("hist e");
  h1peta->SetLineColor(kBlack);
  TH1* h1pheta = mH1F_EpdPhEta->DrawNormalized("hist e same");
  h1pheta->SetLineColor(kBlue);
  TH1* h1cheta = mH1F_EpdChEta->DrawNormalized("hist e same");
  h1cheta->SetLineColor(kGreen+2);

  //canv->cd(5);
  canv->cd(5)->SetLogy();
  TH1* h1pen = mH1F_BestPi0En->DrawNormalized("hist e");
  h1pen->SetLineColor(kBlack);
  TH1* h1phen = mH1F_EpdPhEn->DrawNormalized("hist e same");
  h1phen->SetLineColor(kBlue);
  TH1* h1chen = mH1F_EpdChEn->DrawNormalized("hist e same");
  h1chen->SetLineColor(kGreen+2);
  
  canv->cd(6)->SetLogy();
  TH1* h1ppt = mH1F_BestPi0Pt->DrawNormalized("hist e");
  h1ppt->SetLineColor(kBlack);
  TH1* h1phpt = mH1F_EpdPhPt->DrawNormalized("hist e same");
  h1phpt->SetLineColor(kBlue);
  TH1* h1chpt = mH1F_EpdChPt->DrawNormalized("hist e same");
  h1chpt->SetLineColor(kGreen+2);

  canv->cd(7);
  //canv->cd(7)->SetLogy();
  mH1F_BestPi0Mass->SetLineColor(kBlack);
  mH1F_BestPi0Mass->Draw("hist e");
  mH1F_EpdPhInvMass->SetLineColor(kBlue);
  mH1F_EpdPhInvMass->Draw("hist e same");
  mH1F_EpdChInvMass->SetLineColor(kGreen+2);
  mH1F_EpdChInvMass->Draw("hist e same");
  
  canv->cd(8);
  //canv->cd(8)->SetLogy();
  TH1* h1pmass = mH1F_BestPi0Mass->DrawNormalized("hist e");
  h1pmass->GetXaxis()->SetRangeUser(0,0.3);
  h1pmass->SetLineColor(kBlack);
  //h1pmass->Rebin(2);
  TH1* h1phmass = mH1F_EpdPhInvMass->DrawNormalized("hist e same");
  h1phmass->SetLineColor(kBlue);
  //h1phmass->Rebin(2);
  TH1* h1chmass = mH1F_EpdChInvMass->DrawNormalized("hist e same");
  h1chmass->SetLineColor(kGreen+2);

  canv->cd(9);
  //canv->cd(9)->SetLogy();
  /*
  TH1* h1allmass = mH1F_AllPointPairMass->DrawNormalized("hist e");
  h1allmass->GetXaxis()->SetRangeUser(0,0.3);
  h1allmass->SetLineColor(kBlack);
  //h1allmass->Rebin(2);
  TH1* h1phallmass = mH1F_EpdPhAllPoints->DrawNormalized("hist e same");
  h1phallmass->SetLineColor(kBlue);
  //h1phallmass->Rebin(2);  
  TH1* h1challmass = mH1F_EpdChAllPoints->DrawNormalized("hist e same");
  h1challmass->SetLineColor(kGreen+2);
  */
  mH1F_AllPointPairMass->Draw("hist e");
  //mH1F_AllPointPairMass->GetXaxis()->SetRangeUser(0,0.3);
  mH1F_AllPointPairMass->SetLineColor(kBlack);
  mH1F_EpdPhAllPoints->Draw("hist e same");
  mH1F_EpdPhAllPoints->SetLineColor(kBlue);
  mH1F_EpdChAllPoints->Draw("hist e same");
  mH1F_EpdChAllPoints->SetLineColor(kGreen+2);


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
  canv->cd(1);
  mH1F_NFoundPhiBin->Draw("hist e");
  canv->cd(2);
  mH1F_AllCuts_xF->Draw("hist e");
  canv->cd(3);
  mH1F_AllCuts_Pi0En->Draw("hist e");
  canv->cd(4);
  //mH1F_AllCuts_Pi0Phi->Draw("hist e");
  mH2F_AllCuts_Pi0_etaVphi->Draw("colz");
  canv->cd(5);
  mH1F_InvMassAllCuts->Draw("hist e");
  canv->cd(6)->SetLogy();
  mH1F_NBadEpdProj->Draw("hist e");
  mH1F_NBadEpdProjVcut->Draw("hist e same");
  canv->cd(7);
  mH2F_AllCuts_Pi0_yVx->Draw("colz");
  canv->cd(8);
  mH1F_Pi0MultAllCuts->Draw("hist e");
  canv->cd(9);
  TH1* h1f_multcut_clone = (TH1*)mH1F_Pi0MultAllCuts->Clone("h1f_multcut_clone");
  h1f_multcut_clone->SetBinContent(1,0); //Artifically delete all entries in zero bin to zoom in on the higher regions
  h1f_multcut_clone->Draw("hist e");
  canv->Print(savename);
}

void StMuFcsPi0TreeMaker::PaintInvMassCuts(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  for( short ebin=0; ebin<NENERGYBIN; ++ebin ){
    canv->DivideSquare(NPHIBIN);
    for( short phibin=0; phibin<NPHIBIN; ++phibin ){
      canv->cd(phibin+1);
      std::stringstream histname;
      histname << "H1F_InvMass_en"<<ebin << "_phi"<<phibin;
      TH1D* hist_proj = ((TH3*)mH3F_AllCutsInvMass_enVphi)->ProjectionZ( histname.str().c_str(), phibin+1,phibin+1, ebin+1,ebin+1 );
      hist_proj->SetTitle( histname.str().c_str() );
      hist_proj->Draw("hist e");
      //mH1F_InvMassAllCutsByEnByPhi[ebin][phibin]->Draw("hist e");
    }
    canv->Print(savename);
    canv->Clear();  //Should hopefully delete the projection histograms 
  }
}

void StMuFcsPi0TreeMaker::PaintNpi0(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(2,2);
  int ipad = 1;
  for( int ibeam=0; ibeam<2; ++ibeam ){
    for( int ispin=0; ispin<2; ++ispin ){
      canv->cd(ipad++);
      mH2F_NPi0_enVphi[ibeam][ispin]->Draw("colz");
    }
  }
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
  mGE_AllCuts_InvMass->SetPoint(irun,irun,mH1F_InvMassAllCuts->GetMean());
  mGE_AllCuts_InvMass->SetPointError(irun,0,mH1F_InvMassAllCuts->GetRMS());
  //delete pi0gausfit;
  //pi0gausfit = 0;

  mGE_AllCuts_Pi0En->SetPoint(irun,irun,mH1F_AllCuts_Pi0En->GetMean());
  mGE_AllCuts_Pi0En->SetPointError(irun,0,mH1F_AllCuts_Pi0En->GetRMS());

  return;
}


void StMuFcsPi0TreeMaker::MergeForTssa( TH1* totalhist[][2], TH3* mergedinvmass )
{
  if( mH1F_RndmSpin->GetBinContent(1)>0.1 ){ return; } //Don't merge histograms from files with random spin patterns
  for( int ibeam=0; ibeam<2; ++ibeam ){
    for( int ispin=0; ispin<2; ++ispin ){
      totalhist[ibeam][ispin]->Add(mH2F_NPi0_enVphi[ibeam][ispin]);
    }
  }
  mergedinvmass->Add(mH3F_AllCutsInvMass_enVphi);
}

void StMuFcsPi0TreeMaker::DoTssaAna( TH1* npi0[][2], TH1* h1_rawasymVphi[][StMuFcsPi0TreeMaker::NENERGYBIN] )
{
  Double_t pi = TMath::Pi();
  //TH1* h1_rawasymVphi[2][NENERGYBIN];  //Histograms of computed raw asymmetries for blue and yellow beams for each energy bin
  //memset(h1_rawasymVphi,0,sizeof(h1_rawasymVphi));
  for( int ibeam = 0; ibeam<2; ++ibeam ){
    for( int iebin = 0; iebin<NENERGYBIN; ++iebin ){
      std::stringstream ss_hname;
      ss_hname << "H1_rawasymVphi_"<< (ibeam==0?"Blue":"Yellow") << "_En" << iebin;
      std::stringstream ss_htitle;
      ss_htitle << "Raw Asymmetry vs. Phi for "<<  (ibeam==0?"Blue":"Yellow") << " Beam Ebin=" << iebin << ";phi;Raw AN";
      h1_rawasymVphi[ibeam][iebin] = new TH1F(ss_hname.str().c_str(),ss_htitle.str().c_str(),NPHIBIN/2,-pi/2.0,pi/2.0);
      //AnAtPhi->SetNameTitle("GE_AnAtPhi","A_N for Pi0;cos(phi);A_N");
      for( int iphibin=0; iphibin<NPHIBIN/2 ; ++iphibin ){
	Double_t nupdown = sqrt(npi0[ibeam][0]->GetBinContent(npi0[ibeam][0]->GetBin(iphibin+1,iebin+1)) * npi0[ibeam][1]->GetBinContent(npi0[ibeam][1]->GetBin(iphibin+1+NPHIBIN/2,iebin+1)));
	Double_t ndownup = sqrt(npi0[ibeam][1]->GetBinContent(npi0[ibeam][1]->GetBin(iphibin+1,iebin+1)) * npi0[ibeam][0]->GetBinContent(npi0[ibeam][0]->GetBin(iphibin+1+NPHIBIN/2,iebin+1)));
	Double_t numerator = nupdown - ndownup;
	Double_t denominator = nupdown + ndownup;
	//std::cout << "|iphibin:"<<iphibin << "|iebin:"<<iebin << "|ibeam:"<<ibeam << "|numer:"<<numerator << "|denom:"<<denominator << std::endl;
	if( denominator != 0){
	  h1_rawasymVphi[ibeam][iebin]->SetBinContent(iphibin+1,numerator/denominator);
	}
	else{ h1_rawasymVphi[ibeam][iebin]->SetBinContent(iphibin+1,0); } //If denomnator is zero set bin to zero to avoid bad division
      }
    }
  }
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


