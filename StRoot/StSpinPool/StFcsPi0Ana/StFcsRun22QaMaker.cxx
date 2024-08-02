#include "StEvent/StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEventTypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StSpinPool/StFcsRawDaqReader/StFcsRawDaqReader.h"
#include "StThreeVectorF.hh"
#include "Stypes.h"

#include "StFcsRun22QaMaker.h"

ClassImp(StFcsRun22QaMaker)


StFcsRun22QaMaker::StFcsRun22QaMaker(const char* name):StMaker(name)
{
  mFileName = "";
  mSpinRndm.SetSeed(0);
  mAllHists = new TObjArray();

  memset(mH2F_Hit_adcVtb,0,sizeof(mH2F_Hit_adcVtb));
  memset(mH2F_Hit_enVid,0,sizeof(mH2F_Hit_enVid));
  memset(mH2F_Hit_fitpeakVid,0,sizeof(mH2F_Hit_fitpeakVid));
  memset(mH2F_Hit_chi2Vid,0,sizeof(mH2F_Hit_chi2Vid));
  memset(mH2F_Hit_npeaksVid,0,sizeof(mH2F_Hit_npeaksVid));
  memset(mH1F_Hit_NHits,0,sizeof(mH1F_Hit_NHits));
  memset(mH1F_Hit_ESum,0,sizeof(mH1F_Hit_ESum));

  memset(mH2F_HitPres_depVqt,0,sizeof(mH2F_HitPres_depVqt));
  memset(mH2F_HitPres_peakVtac,0,sizeof(mH2F_HitPres_peakVtac));
  memset(mH2F_HitEpd_tacVadcmip,0,sizeof(mH2F_HitEpd_tacVadcmip));
  memset(mH2F_HitEpd_nmipVchkey,0,sizeof(mH2F_HitEpd_nmipVchkey));

  memset(mH1F_NClusters,0,sizeof(mH1F_NClusters));
  memset(mH1F_Clu_NTowers,0,sizeof(mH1F_Clu_NTowers));
  memset(mH1F_Clu_NNei,0,sizeof(mH1F_Clu_NNei));
  memset(mH1F_Clu_NPoints,0,sizeof(mH1F_Clu_NPoints));
  memset(mH1F_Clu_En,0,sizeof(mH1F_Clu_En));
  memset(mH2F_Clu_yVx,0,sizeof(mH2F_Clu_yVx));
  memset(mH2F_Clu_sigmaxVsigmin,0,sizeof(mH2F_Clu_sigmaxVsigmin));
  memset(mH1F_Clu_theta,0,sizeof(mH1F_Clu_theta));
  memset(mH2F_Clu_Chi2NdfPhoton_2V1,0,sizeof(mH2F_Clu_Chi2NdfPhoton_2V1));

  memset(mH1F_NPoints,0,sizeof(mH1F_NPoints));
  memset(mH1F_Poi_En,0,sizeof(mH1F_Poi_En));
  memset(mH1F_Poi_NCluPhotons,0,sizeof(mH1F_Poi_NCluPhotons));
  memset(mH2F_Poi_yVx,0,sizeof(mH2F_Poi_yVx));
}

StFcsRun22QaMaker::~StFcsRun22QaMaker()
{
  delete mAllHists;
  delete mFileOutput;
}

Short_t StFcsRun22QaMaker::getRandomSpin()
{
  if( mSpinRndm.Rndm()<0.5 ){ return -1; }
  else{ return 1; }
}

void StFcsRun22QaMaker::spinFrom4BitSpin( int spin4bit, int& bpol, int& ypol )
{
  if( bpol!=0 || ypol!=0 ){ bpol=0; ypol=0; }
  if( spin4bit & 0x1 ) ypol = +1;
  if( spin4bit & 0x2 ) ypol = -1;
  if( spin4bit & 0x4 ) bpol = +1;
  if( spin4bit & 0x8 ) bpol = -1;
}

UInt_t StFcsRun22QaMaker::LoadHists(TFile* file)
{
  UInt_t loaded = 0;
  loaded += AddH1F(file,mH1F_Entries,"H1F_Entries","Entries",2,0,2);
  
  loaded += AddH1F(file,mH1F_Triggers,"H1F_Triggers","Triggers;TrigId",999,890000,890999);//@[June 3, 2024] > This is almost all of them as some ids are not in this range but good enough for now
  
  loaded += AddH1F(file,mH1F_VertexPrimZ,"H1F_VertexPrimZ","Primary Vertex (z);cm",201,-100.5,100.5);
  loaded += AddH1F(file,mH1F_VertexVpd,"H1F_VertexVpd","Vpd Vertex (z);cm",50,-200,200);
  loaded += AddH1F(file,mH1F_VertexBbc,"H1F_VertexBbc","Bbc Vertex (z);cm",50,-200,200);
  loaded += AddH1F(file,mH1F_BbcTimeDiff,"H1F_BbcTimeDiff","Bbc Time difference", 200,-5000,5000);
  loaded += AddH1F(file,mH1F_VertexZdc,"H1F_VertexZdc","Zdc Vertex (z);cm",50,-200,200);
  loaded += AddH1F(file,mH1F_VertexEpd,"H1F_VertexEpd","Epd Vertex (z);cm",50,-200,200);
  
  loaded += AddH2F(file,mH2F_BxId_7V48,"H2F_BxId7V48","Bunch Crossing Id;48 bit;7 bit", 121,-0.5,120.5, 121,-0.5,120.5);
  loaded += AddH2F(file,mH2F_Mult_tofVref,"H2F_Mult_tofVref","TOF multiplicity vs. Reference multiplicity;RefMult;TofMult",100,0,100,100,0,100);
  loaded += AddH2F(file,mH2F_Mult_tofVecal,"H2F_Mult_tofVecal","TOF multiplicity vs. Fcs Ecal multiplicity;EcalMult;TofMult",200,0,1000,100,0,100);

  loaded += AddH1F(file,mH1F_Spin,"H1F_Spin","Spin",3,-1.5,1.5);

  TString namesuffix[6] = {"EN","ES","HN","HS","PN","PS"};
  for( UShort_t i=0; i<kFcsNDet; ++i ){
    UInt_t nchs = 0;
    UInt_t ncol = 0;
    UInt_t nrow = 0;
    if( i<=kFcsEcalSouthDetId )                              { nchs = kFcsEcalMaxId; ncol = kFcsEcalNCol; nrow = kFcsEcalNRow; }
    else if( kFcsHcalNorthDetId<=i && i<=kFcsHcalSouthDetId ){ nchs = kFcsHcalMaxId; ncol = kFcsHcalNCol; nrow = kFcsHcalNRow; }
    else if( kFcsPresNorthDetId<=i && i<=kFcsPresSouthDetId ){ nchs = kFcsPresMaxId; ncol = kFcsPresNCol; nrow = kFcsPresNRow; }
    else{ LOG_ERROR << "StFcsFun22QaMaker::LoadHists() too many detector ids for some reason" << endm; return kStErr; }
    TString histname = "H2F_Hit_adcVtb_" + namesuffix[i];
    TString histtitle = "Adc vs. Tb for " + namesuffix[i] + ";tb;adc";
    if( mFcsAdcTbOn ){
      mH2F_Hit_adcVtb[i] = new TObjArray(); //Create new array for each detector
      loaded += AddH2FArr(file,(mH2F_Hit_adcVtb[i]),nchs,histname.Data(),histtitle.Data(),100,0,100,400,0,4000);
    }
    histname = "H2F_Hit_enVid_" + namesuffix[i];
    histtitle = "Energy vs. Id for " + namesuffix[i] + ";id;energy (GeV)";    
    loaded += AddH2F(file,(mH2F_Hit_enVid[i]),histname.Data(),histtitle.Data(),nchs,0,nchs,100,0,100);
    histname = "H2F_Hit_fitpeakVid_" + namesuffix[i];
    histtitle = "Fitted Peak Location vs. Id for " + namesuffix[i] + ";id;tb";
    loaded += AddH2F(file,(mH2F_Hit_fitpeakVid[i]),histname.Data(),histtitle.Data(),nchs,0,nchs,100,0,100);
    histname = "H2F_Hit_chi2Vid_" + namesuffix[i];
    histtitle = "Chi^2/NDF for fitted peaks (npeaks>1) vs. Id for " + namesuffix[i] + ";id;Chi^2/NDF";
    loaded += AddH2F(file,(mH2F_Hit_chi2Vid[i]),histname.Data(),histtitle.Data(),nchs,0,nchs,110,-10,100);
    histname = "H2F_Hit_npeaksVid_" + namesuffix[i];
    histtitle = "Number of fitted peaks vs. Id for " + namesuffix[i] + ";id;NPeaks";
    loaded += AddH2F(file,(mH2F_Hit_npeaksVid[i]),histname.Data(),histtitle.Data(),nchs,0,nchs,15,0,15);
    histname = "H1F_Hit_Nhits_" + namesuffix[i];
    histtitle = "Hit Multiplicity for " + namesuffix[i] + ";NChs";
    loaded += AddH1F(file,(mH1F_Hit_NHits[i]),histname.Data(),histtitle.Data(),nchs,0,nchs);

    histname = "H1F_NClusters_" + namesuffix[i];
    histtitle = "Number of Clusters for " + namesuffix[i] + ";NClusters";
    loaded += AddH1F(file,(mH1F_NClusters[i]),histname.Data(),histtitle.Data(),50,0,50);
    histname = "H1F_Clu_NTowers_" + namesuffix[i];
    histtitle = "Number of towers in a cluster for " + namesuffix[i] + ";NTowers/cluster";
    loaded += AddH1F(file,(mH1F_Clu_NTowers[i]),histname.Data(),histtitle.Data(),100,0,100);
    histname = "H1F_Clu_NNei_" + namesuffix[i];
    histtitle = "Number of neighbor clusters for a given cluster for " + namesuffix[i] + ";NNeighbors/cluster";
    loaded += AddH1F(file,(mH1F_Clu_NNei[i]),histname.Data(),histtitle.Data(),10,0,10);
    histname = "H1F_Clu_NPoints_" + namesuffix[i];
    histtitle = "Number of points in a cluster for " + namesuffix[i] + ";NPoints/cluster";
    loaded += AddH1F(file,(mH1F_Clu_NPoints[i]),histname.Data(),histtitle.Data(),4,0,4);
    histname = "H1F_Clu_En_" + namesuffix[i];
    histtitle = "Energy in a cluster for " + namesuffix[i] + ";energy (GeV)";
    loaded += AddH1F(file,(mH1F_Clu_En[i]),histname.Data(),histtitle.Data(),100,0,100);
    histname = "H2F_Clu_yVx_" + namesuffix[i];
    histtitle = "Cluster position in row vs. column space for " + namesuffix[i] + ";col;row";
    loaded += AddH2F(file,(mH2F_Clu_yVx[i]),histname.Data(),histtitle.Data(), ncol+1,1,ncol+1, nrow+1,1,nrow+1);
    histname = "H2F_Clu_sigmaxVsigmin_" + namesuffix[i];
    histtitle = "Cluster Sigma Max vs. Sigma Min for " + namesuffix[i] + ";sigma min;sigma max";
    loaded += AddH2F(file,(mH2F_Clu_sigmaxVsigmin[i]),histname.Data(),histtitle.Data(), 100,0,1, 100,0,1);
    histname = "H1F_Clu_theta_" + namesuffix[i];
    histtitle = "Cluster theta for " + namesuffix[i];
    loaded += AddH1F(file,(mH1F_Clu_theta[i]),histname.Data(),histtitle.Data(), 20,-1.0*TMath::Pi(),TMath::Pi());
    histname = "H2F_Clu_Chi2NdfPhoton_2V1_" + namesuffix[i];
    histtitle = "Cluster Chi^2/NDF for 2 photon fit vs. 1 photon fit for " + namesuffix[i] + ";1 photon Chi^2/NDF;2 photon Chi^2/NDF";
    loaded += AddH2F(file,(mH2F_Clu_Chi2NdfPhoton_2V1[i]),histname.Data(),histtitle.Data(), 100,0,100, 100,0,100);

    histname = "H1F_NPoints_" + namesuffix[i];
    histtitle = "Number of points for " + namesuffix[i] + ";NPoints";
    loaded += AddH1F(file,(mH1F_NPoints[i]),histname.Data(),histtitle.Data(), 50,0,50);
    histname = "H1F_Poi_En_" + namesuffix[i];
    histtitle = "Point Energy for " + namesuffix[i] + ";Energy (GeV)";
    loaded += AddH1F(file,(mH1F_Poi_En[i]),histname.Data(),histtitle.Data(), 100,0,100);
    histname = "H1F_Poi_NCluPhotons_" + namesuffix[i];
    histtitle = "Number of photons from parent cluster for " + namesuffix[i] + ";NPoints in Parent Cluster";
    loaded += AddH1F(file,(mH1F_Poi_NCluPhotons[i]),histname.Data(),histtitle.Data(), 4,0,4);
    histname = "H2F_Poi_yVx_" + namesuffix[i];
    histtitle = "Point positions for " + namesuffix[i] + ";x (cm);y (cm)";
    loaded += AddH2F(file,(mH2F_Poi_yVx[i]),histname.Data(),histtitle.Data(), ncol+1,1,ncol+1, nrow+1,1,nrow+1 );
  }
  loaded += AddH1F(file,mH1F_Hit_ESum[0],"H1F_Hit_ESum_Ecal","Total energy sum in Ecal;Energy (GeV)", 200,0,200);
  loaded += AddH1F(file,mH1F_Hit_ESum[1],"H1F_Hit_ESum_Hcal","Total energy sum in Hcal;Energy (GeV)", 100,0,100);
  loaded += AddH1F(file,mH1F_Hit_ESum[2],"H1F_Hit_ESum_Pres","Total energy sum in Pres;Energy (GeV)", 100,0,100);

  if( mEpdAdcQaOn || mEpdTacQaOn || mEpdTacAdcOn ){
    loaded += AddH1F(file,mH1F_Epd_NHits,"H1F_Epd_NHits","Number of Hits from EPD collection;NHits",500,0,500);
    loaded += AddH1F(file,mH1F_Epd_NHitsWest,"H1F_Epd_NHitsWest","Number of Hits form EPD collection (West);NHits",300,0,300);
    if( mEpdAdcQaOn ){
      mH2F_HitPres_depVqt[0]   = new TObjArray();
      mH2F_HitPres_depVqt[1]   = new TObjArray();
      loaded += AddH2FArr(file,mH2F_HitPres_depVqt[0],192,"H2F_HitPres_depVqt_PN","QT sum vs. DEP ADC sum for Fcs Preshower North hits;QtSum;DepSum", 64,0,4096, 64,0,4096);
      loaded += AddH2FArr(file,mH2F_HitPres_depVqt[1],192,"H2F_HitPres_depVqt_PS","QT sum vs. DEP ADC sum for Fcs Preshower South hits;QtSum;DepSum", 64,0,4096, 64,0,4096);
    }
    if( mEpdTacQaOn ){
      mH2F_HitPres_peakVtac[0] = new TObjArray();
      mH2F_HitPres_peakVtac[1] = new TObjArray();
      loaded += AddH2FArr(file,mH2F_HitPres_peakVtac[0],192,"H2F_HitPres_peakVtac_PN","Qt TAC vs. Found peak tb for Fcs Preshower North hits;peak (tb);Qt TAC", 100,0,100, 100,0,100);
      loaded += AddH2FArr(file,mH2F_HitPres_peakVtac[1],192,"H2F_HitPres_peakVtac_PS","Qt TAC vs. Found peak tb for Fcs Preshower South hits;peak (tb);Qt TAC", 100,0,100, 100,0,100); 
    }
    if( mEpdTacAdcOn ){
      mH2F_HitEpd_tacVadcmip[0] = new TObjArray(); //Create new array for east side
      mH2F_HitEpd_tacVadcmip[1] = new TObjArray(); //Create new array for west side
      //EPD has 372 tiles on one side
      loaded += AddH2FArr(file,mH2F_HitEpd_tacVadcmip[0],372,"H2F_HitEpd_tacVadcmip_E","Qt TAC vs. ADC/ADC_1mip;ADC/ADC_1mip;TAC", 50,0,25, 200,0,4000);
      loaded += AddH2FArr(file,mH2F_HitEpd_tacVadcmip[1],372,"H2F_HitEpd_tacVadcmip_W","Qt TAC vs. ADC/ADC_1mip;ADC/ADC_1mip;TAC", 50,0,25, 200,0,4000 );
    }

    loaded += AddH2F(file,mH2F_HitEpd_nmipVchkey[0],"H2F_HitEpd_nmipVchkey_E","NMIP values for East EPD channels by key;key ((pp-1)*31+(tt-1));nmip", 372,0,372, 100,0,50);
    loaded += AddH2F(file,mH2F_HitEpd_nmipVchkey[1],"H2F_HitEpd_nmipVchkey_W","NMIP values for West EPD channels by key;key ((pp-1)*31+(tt-1));nmip", 372,0,372, 100,0,50);
    
    loaded += AddH2F(file,mH2F_Epd_earlywVearlye,"H2F_Epd_earlywVearlye","EPD Earliest TAC West vs. East;Earliest East TAC;Earliest West TAC", 300,0,4200, 300,0,4200);
    loaded += AddH2F(file,mH2F_Epd_avgwVavge,"H2F_Epd_avgwVavge","EPD Average TAC West vs. East;Average East TAC;Average West TAC", 200,0,1000, 200,0,1000);
    loaded += AddH1F(file,mH1F_EpdTacDiff_Early,"H1F_EpdTacDiff_Early","Epd TAC difference from Earliest TAC;TacDiff Early",200,-3000,3000);
    loaded += AddH1F(file,mH1F_EpdTacDiff_Avg,"H1F_EpdTacDiff_Avg","Epd TAC difference from Average TAC;TacDiff Avg",200,-3000,3000);
    
    loaded += AddH2F(file,mH2F_EpdCut_earlywVearlye,"H2F_EpdCut_earlywVearlye","EPD Earliest TAC West vs. East with 1<nMIP<15;Earliest East TAC;Earliest West TAC", 300,0,4200, 300,0,4200);
    loaded += AddH2F(file,mH2F_EpdCut_avgwVavge,"H2F_EpdCut_avgwVavge","EPD Average TAC West vs. East with 1<nMIP<15;Average East TAC;Average West TAC", 200,0,2000, 200,0,2000);
    loaded += AddH1F(file,mH1F_EpdCutTacDiff_Early,"H1F_EpdCutTacDiff_Early","Epd TAC difference from Earliest TAC and 1<nMIP<15;TacDiff Early",200,-3000,3000);
    loaded += AddH1F(file,mH1F_EpdCutTacDiff_Avg,"H1F_EpdCutTacDiff_Avg","Epd TAC difference from Average TAC and 1<nMIP<15;TacDiff Avg",200,-3000,3000);
  }
  
  loaded += AddH2F(file,mH2F_CluHigh_angleVesum,"H2F_CluHigh_angleVesum", "Highest two energy clusters opening angle vs. total cluster energy;esum (GeV);opening angle",100,0,100, 60,0,TMath::Pi());
  loaded += AddH2F(file,mH2F_CluHighEn_lowVhigh,"H2F_CluHighEn_lowVhigh","Highest two energy clusters energies energy 1 vs. energy 2;E1 (GeV); E2(GeV)", 100,0,100, 100,0,100);
  loaded += AddH2F(file,mH2F_CluHigh_dggVesum,"H2F_CluHigh_dggVesum","Highest two energy clusters Dgg vs. energy sum;esum (GeV);Dgg (cm)", 100,0,100, 100,0,100);
  loaded += AddH2F(file,mH2F_CluHigh_invmassVesum,"H2F_CluHigh_invmassVesum","Highest two energy clusters invariant mass vs. energy sum;esum (GeV);invariant mass (GeV/c^2)", 100,0,100, 100,0,1.0);
  loaded += AddH2F(file,mH2F_CluHigh_invmassVdgg,"H2F_CluHigh_invmassVdgg","Highest two energy clusters invariant mass vs. Dgg;Dgg (cm);invariant mass (GeV/c^2)", 100,0,100, 100,0,1.0);
  loaded += AddH2F(file,mH2F_CluHigh_invmassVzgg,"H2F_CluHigh_invmassVzgg","Highest two energy clusters invariant mass vs. Zgg;Zgg |E1-E2|/(E1+E2);invariant mass (GeV/c^2)", 100,0,1.0, 100,0,1.0);

  loaded += AddH2F(file,mH2F_PoiHigh_angleVesum,"H2F_PoiHigh_angleVesum", "Highest two energy points opening angle vs. total point energy;esum (GeV);opening angle",100,0,100, 60,0,TMath::Pi());
  loaded += AddH2F(file,mH2F_PoiHighEn_lowVhigh,"H2F_PoiHighEn_lowVhigh","Highest two energy points energies energy 1 vs. energy 2;E1 (GeV); E2(GeV)", 100,0,100, 100,0,100);
  loaded += AddH2F(file,mH2F_PoiHigh_dggVesum,"H2F_PoiHigh_dggVesum","Highest two energy points Dgg vs. energy sum;esum (GeV);Dgg (cm)", 100,0,100, 100,0,100);
  loaded += AddH2F(file,mH2F_PoiHigh_invmassVesum,"H2F_PoiHigh_invmassVesum","Highest two energy points invariant mass vs. energy sum;esum (GeV);invariant mass (GeV/c^2)", 100,0,100, 100,0,1.0);
  loaded += AddH2F(file,mH2F_PoiHigh_invmassVdgg,"H2F_PoiHigh_invmassVdgg","Highest two energy points invariant mass vs. Dgg;Dgg (cm);invariant mass (GeV/c^2)", 100,0,100, 100,0,1.0);
  loaded += AddH2F(file,mH2F_PoiHigh_invmassVzgg,"H2F_PoiHigh_invmassVzgg","Highest two energy points invariant mass vs. Zgg;Zgg |E1-E2|/(E1+E2);invariant mass (GeV/c^2)", 100,0,1.0, 100,0,1.0);

  return loaded;
}

Int_t StFcsRun22QaMaker::Init()
{
  if( mFileName.Length() == 0){ mFileName="test.root"; } //Ensure a TFile is always created
  mFileOutput = new TFile(mFileName.Data(), "RECREATE");
  UInt_t totalhists = this->LoadHists(0); //This is total of histograms loaded from a file not created. Don't use mFileOutput as you are not trying to load from #mFileOutput
  SetOwner(kTRUE);
  LOG_INFO << "StFcsRun22QaMaker::Init() - Loaded " << totalhists << " histograms" << endm;
  return kStOk;
}

Int_t StFcsRun22QaMaker::InitRun(int runnumber)
{
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  //mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
    LOG_ERROR << "StFcsRun22QaMaker::InitRun Failed to get StFcsDbMaker" << endm;
    return kStFatal;
  }
  // if( !mEpdGeo ){ mEpdGeo = new StEpdGeom(); }
  // else{
  //   LOG_ERROR << "StFcsRun22QaMaker::InitRun - StEpdGeom Exists!" << endm;
  //   return kStFatal;
  // }
  return kStOk;
}

Int_t StFcsRun22QaMaker::Make()
{
  mMuDstMkr = (StMuDstMaker*)GetInputDS("MuDst");
  if( mMuDstMkr==0 ){ LOG_ERROR <<"StFcsRun22QaMaker::Make - !MuDstMkr" <<endm; return kStErr; }
  mMuDst = mMuDstMkr->muDst();
  if( mMuDst==0 ){ LOG_ERROR << "StFcsRun22QaMaker::Make - !MuDst" << endm; return kStErr; }
  mMuEvent = mMuDst->event();
  if( mMuEvent==0 ){ LOG_ERROR <<"StFcsRun22QaMaker::Make - !MuEvent" <<endm; return kStErr; }
  mTrigData = mMuEvent->triggerData();
  if( mTrigData==0 ){ LOG_ERROR <<"StFcsRun22QaMaker::Make - !TrigData" <<endm; return kStErr; }
  mRunInfo = &(mMuEvent->runInfo());
  if( mRunInfo==0 ){ LOG_ERROR <<"StFcsRun22QaMaker::Make - !RunInfo" <<endm; return kStErr; }

  mH1F_Entries->Fill(1);
  
  Int_t infostatus = this->FillEventInfo();
  switch( infostatus ){
  case kStEOF: return kStEOF;
  case kStErr: return kStErr;
  case kStFatal: return kStFatal;
  case kStSkip: return kStSkip;
  case kStStop: return kStStop;
  }
  //std::cout << "Filled Event" << std::endl;
  
  if( mEpdAdcQaOn || mEpdTacQaOn || mEpdTacAdcOn ){
    mMuEpdHits = 0;
    mEpdColl = 0;
    mMuEpdHits = mMuDst->epdHits();
    if( mMuEpdHits!=0 ){ if( mMuEpdHits->GetEntriesFast()==0 ){mMuEpdHits=0;} }//If mMuEpdHits is not zero but has no hits set it to zero so rest of code processes from StEpdHitMaker
    if( mMuEpdHits==0 ){ LOG_INFO << "StFcsRun22QaMaker::Make - No MuEPD hits" << endm;
      mEpdHitMkr = (StEpdHitMaker*)GetMaker("epdHit");
      if( mEpdHitMkr==0 ){ LOG_WARN << "StFcsRun22QaMaker::Make - No StEpdHitMaker(\"epdHit\")" << endm; }
      else{ mEpdColl = mEpdHitMkr->GetEpdCollection(); }
      if( mEpdColl==0 ){ LOG_WARN << "StFcsRun22QaMaker::FillFcsInfo - No Epd hit information found" << endm; mEpdHitMkr=0; }//Set the hit maker back to zero so it can be used as a check that the epd collection doesn't exist
    }
    Int_t epdstatus = this->FillEpdInfo();
    switch( epdstatus ){
    case kStEOF: return kStEOF;
    case kStErr: return kStErr;
    case kStFatal: return kStFatal;
    case kStSkip: return kStSkip;
    case kStStop: return kStStop;
    }
  }
  
  Int_t fcsstatus = this->FillFcsInfo();
  switch( fcsstatus ){
  case kStEOF: return kStEOF;
  case kStErr: return kStErr;
  case kStFatal: return kStFatal;
  case kStSkip: return kStSkip;
  case kStStop: return kStStop;
  }
  //std::cout << "Filled Fcs" << std::endl;

  unsigned int totalecalhits = mMuFcsColl->numberOfHits(kFcsEcalNorthDetId) + mMuFcsColl->numberOfHits(kFcsEcalSouthDetId);
  mH2F_Mult_tofVecal->Fill(totalecalhits,mTrigData->tofMultiplicity());

  //std::cout << "Finished Make" << std::endl;
  if( infostatus==kStWarn || fcsstatus==kStWarn ){ return kStWarn; } //Now check if either returned a warning and if so returning warning
  else{ return kStOk; } //Both were ok
}

Int_t StFcsRun22QaMaker::FillEventInfo()
{
  mH2F_BxId_7V48->Fill(mTrigData->bunchId48Bit(),mTrigData->bunchId7Bit());
  mH2F_Mult_tofVref->Fill(mMuEvent->refMult(),mTrigData->tofMultiplicity());

  //Spin information
  if( mSpinDbMkr==0 ){
    mH1F_Spin->Fill(getRandomSpin()); //Random Blue beam spin
    mH1F_Spin->Fill(getRandomSpin()); //Random Yellow beam spin
  }
  else{
    int spin4bit = mSpinDbMkr->spin4usingBX48( mTrigData->bunchId48Bit() );
    int bluespin, yellowspin;
    spinFrom4BitSpin( spin4bit, bluespin, yellowspin );
    mH1F_Spin->Fill(bluespin);
    mH1F_Spin->Fill(yellowspin);
  }
  
  //Trigger Information
  StMuTriggerIdCollection* TrigMuColl = &(mMuEvent->triggerIdCollection());
  if( !TrigMuColl ){ LOG_ERROR <<"StFcsRun22QaMaker::FillEventInfo - !TrigMuColl" <<endl; return kStErr; }
  const StTriggerId& trgIDs = TrigMuColl->nominal();
  Int_t ntrig = trgIDs.triggerIds().size();
  for( Int_t i=0; i<ntrig; ++i ){ mH1F_Triggers->Fill(trgIDs.triggerIds().at(i));}

  //Vertex Information
  //For vertex one possible function is MuEvent->primaryVertexPosition()
  //To get all possible vertex StMuPrimaryVertex *muprimv = MuDst->primaryVertex(index);
  mH1F_VertexPrimZ->Fill( mMuEvent->primaryVertexPosition().z() );
  if( mMuDst->btofHeader() ){ mH1F_VertexVpd->Fill( mMuDst->btofHeader()->vpdVz() ); }
  //@[April 7, 2021] > No Slewing correction for BBC yet, see StFmsJetMaker2015 in BrightSTAR??
  const float bbcTdiff = mTrigData->bbcTimeDifference() - 4096; //subtract 4096 since 0 means bad event and distribution is Gaussian around 4096
  mH1F_BbcTimeDiff->Fill(bbcTdiff);
  if( fabs(bbcTdiff)>1.e-6 ){ mH1F_VertexBbc->Fill( bbcTdiff * -0.2475 ); } //0.2475 = 0.0165*30/2

  //StZdcTriggerDetector& zdc = mMuEvent->zdcTriggerDetector();
  //std::cout <<"|ZdcV:"<<zdc.vertexZ() << std::endl;
  mH1F_VertexZdc->Fill( mTrigData->zdcVertexZ() );
  
  return kStOk;
}

Int_t StFcsRun22QaMaker::FillFcsInfo()
{
  mMuFcsColl = mMuDst->muFcsCollection();
  if (!mMuFcsColl) { LOG_ERROR << "StFcsRun22QaMaker::Make did not find MuFcsCollection" << endm; return kStErr; }
  
  bool check_fillclu = false;
  bool check_fillpoi = false;

  float bestclu_invmass         = -1;
  float bestclu_totale          = -1;
  float bestclu_dgg             = -1;
  float bestclu_zgg             = -1;
  float bestclu_opening_angle   = -1;
  float bestclu_hightowerenergy = -1;
  float bestclu_lowtowerenergy  = -1;
  
  float bestpoi_invmass         = -1;
  float bestpoi_totale          = -1;
  float bestpoi_dgg             = -1;
  float bestpoi_zgg             = -1;
  float bestpoi_opening_angle   = -1;
  float bestpoi_hightowerenergy = -1;
  float bestpoi_lowtowerenergy  = -1;

  double esum[3] = {0,0,0}; //total energy deposited from all hits in the Ecal, Hcal, and Pres respectively

  TClonesArray* hits = mMuFcsColl->getHitArray();
  if( hits==0 ){ LOG_INFO << "StFcsRun22QaMaker::FillFcsInfo - No FCS hits" << endm; }
  TClonesArray* clusters = mMuFcsColl->getClusterArray();
  if( clusters==0 ){ LOG_INFO << "StFcsRun22QaMaker::FillFcsInfo - No FCS clusters" << endm; }
  TClonesArray* points = mMuFcsColl->getPointArray();
  if( points==0 ){ LOG_INFO << "StFcsRun22QaMaker::FillFcsInfo - No FCS points" << endm; }
  
  //std::cout << "|hits:"<<hits << "|clusters:"<<clusters << "|points:"<<points << std::endl;
    
  for( UInt_t idet=0; idet<kFcsNDet; ++idet ){
    //std::cout << "+ |idet:"<<idet << "|maxdet:"<<kFcsNDet;
    if( hits!=0 ){
      unsigned int nh = mMuFcsColl->numberOfHits(idet);
      mH1F_Hit_NHits[idet]->Fill(nh);
      unsigned int ihit = mMuFcsColl->indexOfFirstHit(idet);
      nh += ihit; //Need to correct for the fact that number of hits is just that and doesn't correspond to max index to loop to
      //std::cout << "|nh:" << mMuFcsColl->numberOfHits(idet) << "|index:"<<mMuFcsColl->indexOfFirstHit(idet) << "|actualmax:"<<nh << std::endl;
      for( ; ihit<nh; ++ihit ){
	StMuFcsHit* hit = (StMuFcsHit*)hits->At(ihit);
	//unsigned short hit_det = hit->detectorId();
	unsigned short hit_ehp = hit->ehp();
	unsigned short hit_id = hit->id();
	float hit_energy = hit->energy();
	float hit_adcsum = hit->adcSum();
	float hit_tbpeak = hit->fitPeak();
	int   hit_npeak  = hit->nPeak();
	float hit_chi2   = hit->fitChi2();
	unsigned int ntb = hit->nTimeBin();
	//std::cout << "   + |ihit:"<<ihit << "|idet:"<<hit_det << "|hit_id:"<<hit_id << "|hit_ehp:"<<hit_ehp << "|hit_adcsum:"<<hit_adcsum << "|hit_energy:"<<hit_energy << "|hit_ntb:"<<ntb << "|nh:"<<nh << "|th:"<<mMuFcsColl->numberOfHits() << std::endl;
	esum[hit_ehp] += hit_energy;
	mH2F_Hit_enVid[idet]     ->Fill(hit_id,hit_energy);
	mH2F_Hit_fitpeakVid[idet]->Fill(hit_id,hit_tbpeak);
	mH2F_Hit_chi2Vid[idet]   ->Fill(hit_id,hit_chi2);
	mH2F_Hit_npeaksVid[idet] ->Fill(hit_id,hit_npeak);
	if( mFcsAdcTbOn ){
	  for( unsigned int i=0; i<ntb; ++i ){ ((TH1*) mH2F_Hit_adcVtb[idet]->UncheckedAt(hit_id))->Fill(hit->timebin(i),hit->adc(i)); }
	}
	if( kFcsPresNorthDetId<=idet && idet<=kFcsPresSouthDetId ){
	  if( mEpdAdcQaOn || mEpdTacQaOn ){
	    unsigned int nepdhits = 0;

	    StSPtrVecEpdHit* epdhits = 0;
	    if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
	    else if( mEpdColl!=0 ){
	      epdhits = &(mEpdColl->epdHits());
	      nepdhits = epdhits->size();
	    }
	    else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A SERIOUS BUG IN THE CODE" << endm; return kStErr; }
	    
	    int fcspp; int fcstt;
	    mFcsDb->getEPDfromId(idet,hit_id,fcspp,fcstt);
	    //For processing epd hits
	    int adc = 0;
	    int tac = 0;
	    StMuEpdHit* muepdhit = 0;
	    StEpdHit* epdhit = 0;
	    for(unsigned int i=0; i<nepdhits; ++i ){
	      if( mMuEpdHits!=0 ){ muepdhit = (StMuEpdHit*)mMuEpdHits->UncheckedAt(i); } //To match similar in StMuDstMaker->epdHit(int i)
	      else if( epdhits!=0 ){ epdhit = (StEpdHit*)(*epdhits)[i]; }
	      else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A VERY SERIOUS BUG IN THE CODE" << endm; return kStErr; }
	      int ew    = muepdhit!=0 ? muepdhit->side()    : epdhit->side();      //east=-1, west=1
	      int epdpp = muepdhit!=0 ? muepdhit->position(): epdhit->position();  //Supersector runs [1,12]
	      int epdtt = muepdhit!=0 ? muepdhit->tile()    : epdhit->tile();      //Tile number [1,31]
	      if( 1==ew ){
		adc = muepdhit!=0 ? muepdhit->adc() : epdhit->adc();
		tac = muepdhit!=0 ? muepdhit->tac() : epdhit->tac();
		if( fcspp==epdpp && fcstt==epdtt ){ break; }
	      }
	    }
	    //Here adc and tac is equal to ADC and TAC from matched epd hit
	    if( mEpdAdcQaOn ){ ((TH1*) mH2F_HitPres_depVqt[idet-kFcsPresNorthDetId]->UncheckedAt(hit_id))->Fill(adc,hit_adcsum); }
	    if( mEpdTacQaOn ){ ((TH1*) mH2F_HitPres_peakVtac[idet-kFcsPresNorthDetId]->UncheckedAt(hit_id))->Fill(tac,hit_tbpeak); }
	  }
	}
      }//fcs hits
    }
    else{ LOG_INFO <<"|hits is empty:"<<hits << endm; }

    if( clusters!=0 ){
      unsigned int nc = mMuFcsColl->numberOfClusters(idet);
      mH1F_NClusters[idet]->Fill(nc);
      unsigned int iclus=mMuFcsColl->indexOfFirstCluster(idet);
      nc += iclus;
      for( ; iclus<nc; ++iclus){
	StMuFcsCluster* clu = (StMuFcsCluster*)clusters->At(iclus);
	float iclu_x = clu->x();
	float iclu_y = clu->y();
	float iclu_energy = clu->energy();
	mH1F_Clu_NTowers[idet]->Fill(clu->nTowers());
	mH1F_Clu_NNei[idet]->Fill(clu->nNeighbor());
	mH1F_Clu_NPoints[idet]->Fill(clu->nPoints());
	mH1F_Clu_En[idet]->Fill(iclu_energy);
	mH2F_Clu_yVx[idet]->Fill(iclu_x,iclu_y);
	mH2F_Clu_sigmaxVsigmin[idet]->Fill(clu->sigmaMin(),clu->sigmaMax());
	mH1F_Clu_theta[idet]->Fill(clu->theta());
	mH2F_Clu_Chi2NdfPhoton_2V1[idet]->Fill(clu->chi2Ndf1Photon(),clu->chi2Ndf2Photon());

	StThreeVectorD iclu_pos = mFcsDb->getStarXYZfromColumnRow( idet, iclu_x, iclu_y );
	StLorentzVectorD iclu_p = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, 0 );

	if( idet<=kFcsEcalSouthDetId ){
	  if( iclus==(nc-1) ){ continue; }
	  for( unsigned int j=iclus+1; j<nc; j++ ){
	    StMuFcsCluster* cluj = (StMuFcsCluster*)clusters->At(j);
	    float jclu_energy = cluj->energy();
	    float jclu_x = cluj->x();
	    float jclu_y = cluj->y();
	    StThreeVectorD jclu_pos = mFcsDb->getStarXYZfromColumnRow( idet, jclu_x, jclu_y );
	    double ensum = iclu_energy + jclu_energy;
	    float zgg = (fabs(iclu_energy - jclu_energy)) / (ensum);
	    StLorentzVectorD jclu_p = mFcsDb->getLorentzVector( jclu_pos, jclu_energy, 0 );
	  
	    if( jclu_energy<mEnCut ){ continue; }
	    if( zgg>0.7 ){ continue; }
	    if( ensum>bestclu_totale ){
	      check_fillclu = true;
	      bestclu_invmass = ((iclu_p + jclu_p).m());
	      bestclu_totale = ensum;
	      bestclu_dgg = sqrt( (iclu_pos[0]-jclu_pos[0])*(iclu_pos[0]-jclu_pos[0]) + (iclu_pos[1]-jclu_pos[1])*(iclu_pos[1]-jclu_pos[1]) + (iclu_pos[2]-jclu_pos[2])*(iclu_pos[2]-jclu_pos[2]) );
	      bestclu_zgg = zgg;
	      double cluidotj = iclu_pos[0]*jclu_pos[0] + iclu_pos[1]*jclu_pos[1] + iclu_pos[2]*jclu_pos[2];          //dot product of vectors for the current cluster position and cluster j position
	      double iclumag = sqrt( iclu_pos[0]*iclu_pos[0] + iclu_pos[1]*iclu_pos[1] + iclu_pos[2]*iclu_pos[2] );  //magnitude of position vector for current cluster
	      double jclumag = sqrt( jclu_pos[0]*jclu_pos[0] + jclu_pos[1]*jclu_pos[1] + jclu_pos[2]*jclu_pos[2] );//magnitude of position vector for cluster j
	      bestclu_opening_angle = acos( cluidotj / (iclumag*jclumag) );
	      if( iclu_energy>jclu_energy ){ bestclu_hightowerenergy = iclu_energy; bestclu_lowtowerenergy = jclu_energy; }
	      else                         { bestclu_hightowerenergy = iclu_energy; bestclu_lowtowerenergy = iclu_energy; }
	    }
	  }//jclu
	}
      }//iclu
    }

    if( points!=0 ){
      unsigned int np = mMuFcsColl->numberOfPoints(idet);
      mH1F_NPoints[idet]->Fill(np);
      unsigned int ipoint=mMuFcsColl->indexOfFirstPoint(idet);
      np += ipoint;
      for( ; ipoint<np; ++ipoint ){
	StMuFcsPoint* point = (StMuFcsPoint*)points->At(ipoint);
	float ipoi_x = point->x();
	float ipoi_y = point->y();
	float ipoi_energy = point->energy();
	mH1F_Poi_En[idet]->Fill(ipoi_energy);
	mH1F_Poi_NCluPhotons[idet]->Fill(point->nParentClusterPhotons());
	mH2F_Poi_yVx[idet]->Fill(ipoi_x,ipoi_y);
	StThreeVectorD ipoi_pos = mFcsDb->getStarXYZfromColumnRow( idet, ipoi_x, ipoi_y );
	StLorentzVectorD ipoi_p = mFcsDb->getLorentzVector(ipoi_pos, ipoi_energy, 0);
      
	if( idet<=kFcsEcalSouthDetId ){
	  if( ipoint==(np-1) ){ continue; }
	  for( unsigned int j=(ipoint+1); j<np; ++j ){
	    StMuFcsPoint* poij = (StMuFcsPoint*)points->At(j);
	    float jpoi_energy = poij->energy();
	    float poiesum = ipoi_energy+jpoi_energy;
	    float zgg = (fabs(ipoi_energy - jpoi_energy)) / (poiesum);
	    StThreeVectorD jpoi_pos = mFcsDb->getStarXYZfromColumnRow(idet, poij->x(), poij->y());
	    StLorentzVectorD jpoi_p = mFcsDb->getLorentzVector(jpoi_pos, jpoi_energy, 0);

	    if( idet<=kFcsEcalSouthDetId ){
	      if( jpoi_energy<mEnCut ){ continue; }
	      if( zgg>=0.7 ){ continue; }
	      if( poiesum>bestpoi_totale ){
		check_fillpoi = true;
		bestpoi_invmass = (ipoi_p + jpoi_p).m();
		bestpoi_totale = poiesum;
		bestpoi_dgg = sqrt( (ipoi_pos[0]-jpoi_pos[0])*(ipoi_pos[0]-jpoi_pos[0]) + (ipoi_pos[1]-jpoi_pos[1])*(ipoi_pos[1]-jpoi_pos[1]) + (ipoi_pos[2]-jpoi_pos[2])*(ipoi_pos[2]-jpoi_pos[2]) );
		bestpoi_zgg = zgg;
		double poiidotj = ipoi_pos[0]*jpoi_pos[0] + ipoi_pos[1]*jpoi_pos[1] + ipoi_pos[2]*jpoi_pos[2];          //dot product of vectors for the current point position and point j position
		double ipoimag = sqrt( ipoi_pos[0]*ipoi_pos[0] + ipoi_pos[1]*ipoi_pos[1] + ipoi_pos[2]*ipoi_pos[2] );  //magnitude of position vector for current point
		double jpoimag = sqrt( jpoi_pos[0]*jpoi_pos[0] + jpoi_pos[1]*jpoi_pos[1] + jpoi_pos[2]*jpoi_pos[2] );//magnitude of position vector for point j
		bestpoi_opening_angle = acos( poiidotj / (ipoimag*jpoimag) );
		if( ipoi_energy>jpoi_energy ){ bestpoi_hightowerenergy = ipoi_energy; bestpoi_lowtowerenergy = jpoi_energy; }
		else                         { bestpoi_hightowerenergy = jpoi_energy; bestpoi_lowtowerenergy = ipoi_energy; }
	      }
	    }
	  }//j point
	}
      }//i point
    }
    
  }//fcs dets

  mH1F_Hit_ESum[0]->Fill(esum[0]);
  mH1F_Hit_ESum[1]->Fill(esum[1]);
  mH1F_Hit_ESum[2]->Fill(esum[2]);

  if( check_fillclu ){
    mH2F_CluHigh_angleVesum->Fill(bestclu_totale,bestclu_opening_angle);
    mH2F_CluHighEn_lowVhigh->Fill(bestclu_lowtowerenergy,bestclu_hightowerenergy);
    mH2F_CluHigh_dggVesum->Fill(bestclu_totale,bestclu_dgg);
    mH2F_CluHigh_invmassVesum->Fill(bestclu_totale,bestclu_invmass);
    mH2F_CluHigh_invmassVdgg->Fill(bestclu_dgg,bestclu_invmass);
    mH2F_CluHigh_invmassVzgg->Fill(bestclu_zgg,bestclu_invmass);
  }
  if( check_fillpoi ){
    mH2F_PoiHigh_angleVesum->Fill(bestpoi_totale,bestpoi_opening_angle);
    mH2F_PoiHighEn_lowVhigh->Fill(bestpoi_lowtowerenergy,bestpoi_hightowerenergy);
    mH2F_PoiHigh_dggVesum->Fill(bestpoi_totale,bestpoi_dgg);
    mH2F_PoiHigh_invmassVesum->Fill(bestpoi_totale,bestpoi_invmass);
    mH2F_PoiHigh_invmassVdgg->Fill(bestpoi_dgg,bestpoi_invmass);
    mH2F_PoiHigh_invmassVzgg->Fill(bestpoi_zgg,bestpoi_invmass);
  }
  
  return kStOk;
}

Int_t StFcsRun22QaMaker::Finish()
{
  if( mFileOutput!=0 ){
    mFileOutput->cd();
    mAllHists->Write();
    return kStOk;
  }
  else{
    LOG_WARN << "StFcsRun22QaMaker::Finish() - No file created because pointer is null" << endm;
    return kStWarn;
  }
}

Int_t StFcsRun22QaMaker::FillEpdInfo()
{
  if( mMuEpdHits==0 && mEpdColl==0 ){
    LOG_WARN << "StFcsRun22QaMaker::FillEpdInfo() has no epd hits" << endm;
    return kStWarn;
  }
  unsigned int nepdhits = 0;
  StSPtrVecEpdHit* epdhits = 0;
  if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
  else if( mEpdColl!=0 ){
    epdhits = &(mEpdColl->epdHits());
    nepdhits = epdhits->size();
  }
  else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A SERIOUS BUG IN THE CODE" << endm; return kStErr; }
  if( mH1F_Epd_NHits!=0 ){ mH1F_Epd_NHits->Fill(nepdhits); }
  
  //For processing epd hits
  //Larger TAC values mean earlier times so a tac of 0 means latest possible time
  int earliesttace  = 0;  //Stores the largest TAC value in EPD East
  int earliesttacw  = 0;  //Stores the largest TAC value in EPD West
  int navgtacw   = 0;     //Number of TAC values summed for EPD West
  int navgtace   = 0;     //Number of TAC values summed for EPD East
  double sumtacw = 0;     //Sum of TAC values for EPD West
  double sumtace = 0;     //Sum of TAC values for EPD East

  int cut_earliesttace  = 0;  //Stores the largest TAC value in EPD East
  int cut_earliesttacw  = 0;  //Stores the largest TAC value in EPD West
  int cut_navgtacw   = 0;     //Number of TAC values summed for EPD West
  int cut_navgtace   = 0;     //Number of TAC values summed for EPD East
  double cut_sumtacw = 0;     //Sum of TAC values for EPD West
  double cut_sumtace = 0;     //Sum of TAC values for EPD East

  int nhitwest   = 0;
  StMuEpdHit* muepdhit = 0;
  StEpdHit* epdhit = 0;
  //std::cout << "|mMuEpdHits:"<<mMuEpdHits << "|mEpdColl:"<<mEpdColl << "|epdhits:"<<epdhits << "|nepdhits:"<<nepdhits << std::endl;
  for(unsigned int i=0; i<nepdhits; ++i ){
    if( mMuEpdHits!=0 ){ muepdhit = (StMuEpdHit*)mMuEpdHits->UncheckedAt(i); } //To match similar in StMuDstMaker->epdHit(int i)
    else if( epdhits!=0 ){ epdhit = (StEpdHit*)((*epdhits)[i]); }
    else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A VERY SERIOUS BUG IN THE CODE" << endm; return kStErr; } 
    //std::cout << "|i:"<<i << "|muepdhit:"<<muepdhit << "|epdhit:"<<epdhit << std::endl;
    int ew    = muepdhit!=0 ? muepdhit->side()    : epdhit->side();      //east=-1, west=1
    int epdpp = muepdhit!=0 ? muepdhit->position(): epdhit->position();  //Supersector runs [1,12]
    int epdtt = muepdhit!=0 ? muepdhit->tile()    : epdhit->tile();      //Tile number [1,31]
    //int adc = muepdhit!=0 ? muepdhit->adc() : epdhit->adc();
    float nmip = muepdhit!=0 ? muepdhit->nMIP(): epdhit->nMIP();
    //float adcnmip = 0;
    //if( nmip>0.0001 ){ adcnmip = static_cast<float>(adc)/nmip; }
    int tac = muepdhit!=0 ? muepdhit->tac() : epdhit->tac();
    //std::cout << "|ew:"<<ew << "|pp:"<<epdpp << "|tt:"<<epdtt << "|adc:"<<adc << "|nmip:"<<nmip <<"tac:"<<tac << std::endl
    if( ew==-1 ){ //This is east side
      if( mEpdTacAdcOn ){ ((TH1*)mH2F_HitEpd_tacVadcmip[0]->UncheckedAt( (epdpp-1)*31+(epdtt-1) ))->Fill(nmip,tac); }
      mH2F_HitEpd_nmipVchkey[0]->Fill((epdpp-1)*31+(epdtt-1),nmip);
      sumtace += tac;
      ++navgtace;
      if( tac > earliesttace ){ earliesttace = tac; }
      if( 1<nmip && nmip<15 ){
	cut_sumtace += tac;
	++cut_navgtace;
	if( tac > cut_earliesttace ){ cut_earliesttace = tac; }
      }
    }
    if( ew==1 ){ //This is west side
      if( mEpdTacAdcOn ){ ((TH1*)mH2F_HitEpd_tacVadcmip[1]->UncheckedAt( (epdpp-1)*31+(epdtt-1) ))->Fill(nmip,tac); }
      mH2F_HitEpd_nmipVchkey[1]->Fill((epdpp-1)*31+(epdtt-1),nmip);
      sumtacw += tac;
      ++navgtacw;
      ++nhitwest;
      if( tac > earliesttacw ){ earliesttacw = tac; }
      if( 1<nmip && nmip<15 ){
	cut_sumtacw += tac;
	++cut_navgtacw;
	if( tac > cut_earliesttacw ){ cut_earliesttacw = tac; }
      }
    }
  }// EPD hit loop
  
  double avgtace = 0;
  double avgtacw = 0;
  double cut_avgtace = 0;
  double cut_avgtacw = 0;
  if( navgtace>0 )    { avgtace = sumtace/static_cast<double>(navgtace); }
  if( navgtacw>0 )    { avgtacw = sumtacw/static_cast<double>(navgtacw); }
  if( cut_navgtace>0 ){ cut_avgtace = cut_sumtace/static_cast<double>(cut_navgtace); }
  if( cut_navgtacw>0 ){ cut_avgtacw = cut_sumtacw/static_cast<double>(cut_navgtacw); }
  
  if( mH1F_Epd_NHits!=0 ){ mH1F_Epd_NHitsWest->Fill(nhitwest); }

  //std::cout << "|earlye:"<<earliesttace << "|earlyw:"<<earliesttacw << "|avge:"<<avgtace << "|avgw:"<<avgtacw << std::endl;
  if( mH2F_Epd_avgwVavge!=0 )    { mH2F_Epd_avgwVavge->Fill( avgtace,avgtacw  ); }
  if( mH2F_Epd_earlywVearlye!=0 ){ mH2F_Epd_earlywVearlye->Fill(earliesttace,earliesttacw); }
  if( mH1F_EpdTacDiff_Early!=0 ) { mH1F_EpdTacDiff_Early->Fill( earliesttacw - earliesttace ); }
  if( mH1F_EpdTacDiff_Avg!=0 )   { mH1F_EpdTacDiff_Avg->Fill( avgtacw - avgtace ); }

  //std::cout << "|cut_earlye:"<<cut_earliesttace << "|cut_earlyw:"<<cut_earliesttacw << "|cut_avge:"<<cut_avgtace << "|cut_avgw:"<<cut_avgtacw << std::endl;
  if( mH2F_EpdCut_avgwVavge!=0 )    { mH2F_EpdCut_avgwVavge->Fill( cut_avgtace,cut_avgtacw  ); }
  if( mH2F_EpdCut_earlywVearlye!=0 ){ mH2F_EpdCut_earlywVearlye->Fill(cut_earliesttace,cut_earliesttacw); }
  if( mH1F_EpdCutTacDiff_Early!=0 ) { mH1F_EpdCutTacDiff_Early->Fill( cut_earliesttacw - cut_earliesttace ); }
  if( mH1F_EpdCutTacDiff_Avg!=0 )   { mH1F_EpdCutTacDiff_Avg->Fill( cut_avgtacw - cut_avgtace ); }

  mH1F_VertexEpd->Fill( (avgtacw-avgtace) * -0.2475 );
  //@[Julye 17, 2024] > Using the same logic as BBC since haven't looked at EpdTacDiff histograms, also because EPD has same 15.6ps/TAC as BBC.
  //@[July 21, 2024]>After looking at TAC differences it seems average with no cuts gives actual results and ADC_Nmip cut tac values do not.
  //@[July 23, 2024]>Looking at TAC diff don't need to subtract 4096 like we do for BBC
  return kStOk;
}

UInt_t StFcsRun22QaMaker::AddH1F(TFile* file, TH1*& h1, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh)
{
  UInt_t status = 0;
  if( h1!=0 ){
    //Here I am using bit 22 since that is unused by TH1
    if( h1->TestBit(22) ){ h1=0; } //It is true that the file was loaded so safe to change pointer without delete
    else{ mAllHists->Add(h1); return status; } //Object was new so stop and return 0 and re-add to array
  }
  if( file!=0 ){ h1 = (TH1F*)file->Get(name); }
  if( h1==0 ){
    h1 = new TH1F(name,title,nbins,xlow,xhigh);
    h1->Sumw2();
  }
  else{
    h1->SetBit(22);
    ++status;
  }
  h1->SetTitle(title);
  mAllHists->Add(h1);
  return status;//1 if histogram loaded or exists, 0 otherwise
}

UInt_t StFcsRun22QaMaker::AddH1FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh)
{
  UInt_t status = 0;
  for( UInt_t iobj = 0; iobj<nobjs; ++iobj ){
    std::stringstream ss_name;
    ss_name << name << "_" << iobj;
    TH1F* h1 = 0;
    if( file!=0 ){ h1 = (TH1F*)file->Get(ss_name.str().c_str()); }
    if( h1==0 ){
      h1 = new TH1F(ss_name.str().c_str(),title,nbins,xlow,xhigh);
      h1->Sumw2();
    }
    else{
      h1->SetBit(22);
      ++status;
    }
    h1->SetTitle(title);
    arr->Add(h1);
    mAllHists->Add(h1);
  }
  return status;
}

UInt_t StFcsRun22QaMaker::AddH2F(TFile* file, TH1*& h2, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh)
{
  UInt_t status = 0;
  if( h2!=0 ){
    //Here I am using bit 22 since that is unused by TH1
    if( h2->TestBit(22) ){ h2=0; } //It is true that the file was loaded so safe to change pointer without delete
    else{ mAllHists->Add(h2); return status; } //Object was new so stop and return 0 and re-add to array
  }
  if( file!=0 ){ h2 = (TH2F*)file->Get(name); }
  if( h2==0 ){
    h2 = new TH2F(name,title, nbinsx,xlow,xhigh, nbinsy,ylow,yhigh);
    h2->Sumw2();
  }
  else{
    h2->SetBit(22);
    ++status;
  }
  h2->SetTitle(title);
  mAllHists->Add(h2);
  return status;//1 if histogram loaded, 0 if new
}

UInt_t StFcsRun22QaMaker::AddH2FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh)
{
  UInt_t status = 0;
  for( UInt_t iobj = 0; iobj<nobjs; ++iobj ){
    std::stringstream ss_name;
    ss_name << name << "_" << iobj;
    TH2F* h2 = 0;
    if( file!=0 ){ h2 = (TH2F*)file->Get(ss_name.str().c_str()); }
    if( h2==0 ){
      h2 = new TH2F(ss_name.str().c_str(),title,nbinsx,xlow,xhigh, nbinsy,ylow,yhigh);
      h2->Sumw2();
    }
    else{
      h2->SetBit(22);
      ++status;
    }
    h2->SetTitle(title);
    arr->Add(h2);
    mAllHists->Add(h2);
  }
  return status; 
}

/*
virtual void Paint(Option_t opt="")
{
  TString option(opt);
  option.ToUpper();
  if( option.Contains("E") ){ PaintEventInfo(); }
  if( option.Contains("A") ){ PaintAdcVTb(); }
  if( option.Contains("H") ){ PaintFcsHitQa(); }
  if( option.Contains("C") ){ PaintFcsClusterQa(); }
  if( option.Contains("M") ){ PaintFcsClusterPi0(); }
  if( option.Contains("P") ){ PaintFcsPointQa(); }
  if( option.Contains("R") ){ PaintFcsPointPi0(); }
}
*/
void StFcsRun22QaMaker::DrawEventInfo(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(4,3);
  canv->cd(1);
  mH1F_Entries->Draw("hist e");
  canv->cd(2);
  mH1F_Triggers->Draw("hist e p");
  canv->cd(3);
  mH1F_VertexPrimZ->Draw("hist e");
  canv->cd(4);
  mH1F_VertexVpd->Draw("hist e");
  canv->cd(5);
  mH1F_VertexBbc->Draw("hist e");
  canv->cd(6);
  mH1F_BbcTimeDiff->Draw("hist e");
  canv->cd(7)->SetLogy(1);
  mH1F_VertexZdc->Draw("hist e");
  canv->cd(8);
  mH1F_VertexEpd->Draw("hist e");
  canv->cd(9);
  mH2F_BxId_7V48->Draw("colz");
  canv->cd(10);
  //mH2F_Mult_tofVref->Draw("colz");
  TH1* tofmult = ((TH2*)mH2F_Mult_tofVecal)->ProjectionY("tofmult");
  tofmult->SetTitle("TOF Multiplicty");
  tofmult->Draw("hist e");
  canv->cd(11);
  mH2F_Mult_tofVecal->Draw("colz");
  canv->cd(12);
  mH1F_Spin->Draw("hist e");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawVertex(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1);
  mH1F_VertexVpd->Draw("hist e");
  canv->cd(2)->SetLogy();
  mH1F_VertexZdc->Draw("hist e");
  if( mH1F_VertexEpd!=0 ){
    canv->cd(3);
    mH1F_VertexEpd->Draw("hist e");
  }
  canv->cd(4);
  mH1F_VertexBbc->Draw("hist e");
  canv->cd(5);
  mH1F_BbcTimeDiff->Draw("hist e");

  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawBxId(TCanvas* canv, const char* savename)
{
  canv->Clear();
  mH2F_BxId_7V48->Draw("colz");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawFcsHitSingle(TCanvas* canv, unsigned int det, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1)->SetLogz(1);
  mH2F_Hit_enVid[det]->Draw("colz");
  canv->cd(2)->SetLogz(1);
  mH2F_Hit_fitpeakVid[det]->Draw("colz");
  canv->cd(3)->SetLogz(1);
  mH2F_Hit_chi2Vid[det]->Draw("colz");
  canv->cd(4)->SetLogz(1);
  mH2F_Hit_npeaksVid[det]->Draw("colz");
  canv->cd(5);
  mH1F_Hit_NHits[det]->Draw("hist e");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawFcsClusterSingle(TCanvas* canv, unsigned int det, const char* savename)
{
  canv->Clear();
  canv->Divide(3,3);
  canv->cd(1);
  mH1F_NClusters[det]->Draw("hist e");
  canv->cd(2);
  mH1F_Clu_NTowers[det]->Draw("hist e");
  canv->cd(3);
  mH1F_Clu_NNei[det]->Draw("hist e");
  canv->cd(4);
  mH1F_Clu_NPoints[det]->Draw("hist e");
  canv->cd(5);
  mH1F_Clu_En[det]->Draw("hist e");
  canv->cd(6);
  mH2F_Clu_yVx[det]->Draw("colz");
  canv->cd(7);
  mH2F_Clu_sigmaxVsigmin[det]->Draw("colz");
  canv->cd(8);
  mH1F_Clu_theta[det]->Draw("hist e");
  canv->cd(9);
  mH2F_Clu_Chi2NdfPhoton_2V1[det]->Draw("colz");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawFcsPointSingle(TCanvas* canv, unsigned int det, const char* savename)
{
  canv->Clear();
  canv->Divide(2,2);
  canv->cd(1);
  mH1F_NPoints[det]->Draw("hist e");
  canv->cd(2);
  mH1F_Poi_En[det]->Draw("hist e");
  canv->cd(3);
  mH1F_Poi_NCluPhotons[det]->Draw("hist e");
  canv->cd(4);
  mH2F_Poi_yVx[det]->Draw("colz");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawFcsTotalE(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(2,2);
  canv->cd(1);
  mH1F_Hit_ESum[0]->Draw("hist e");
  canv->cd(2);
  mH1F_Hit_ESum[1]->Draw("hist e");
  canv->cd(3);
  mH1F_Hit_ESum[2]->Draw("hist e");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawAdcVTb(TCanvas* canv, const char* savename)
{
  if( mFcsAdcTbOn ){
    canv->Clear();
    canv->Divide(5,5);
    for( UInt_t i=0; i<kFcsNDet; ++i ){
      for( Int_t ich=0, ipad=1; ich<mH2F_Hit_adcVtb[i]->GetEntriesFast(); ++ich,++ipad ){
	if( ipad>25 ){ ipad=1; canv->Print(savename); canv->Clear(); canv->Divide(5,5); }
	canv->cd(ipad);
	((TH1*)mH2F_Hit_adcVtb[i]->UncheckedAt(ich))->Draw("colz");
      }
    }
    canv->Print(savename);
  }
}

void StFcsRun22QaMaker::DrawFcsHitQa(TCanvas* canv, const char* savename)
{
  for( UShort_t i=0; i<kFcsNDet; ++i ){
    canv->Clear();
    canv->Divide(3,2);
    canv->cd(1)->SetLogz(1);
    mH2F_Hit_enVid[i]->Draw("colz");
    canv->cd(2)->SetLogz(1);
    mH2F_Hit_fitpeakVid[i]->Draw("colz");
    canv->cd(3)->SetLogz(1);
    mH2F_Hit_chi2Vid[i]->Draw("colz");
    canv->cd(4)->SetLogz(1);
    mH2F_Hit_npeaksVid[i]->Draw("colz");
    canv->cd(5);
    mH1F_Hit_NHits[i]->Draw("hist e");
    canv->Print(savename);
  }
  canv->Clear();
  canv->Divide(2,2);
  canv->cd(1);
  mH1F_Hit_ESum[0]->Draw("hist e");
  canv->cd(2);
  mH1F_Hit_ESum[1]->Draw("hist e");
  canv->cd(3);
  mH1F_Hit_ESum[2]->Draw("hist e");
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawEpdHitQa(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(2,2);
  canv->cd(1);
  if( mH1F_Epd_NHits!=0 ){ mH1F_Epd_NHits->Draw("hist e"); }
  canv->cd(2);
  if( mH1F_Epd_NHitsWest!=0 ){ mH1F_Epd_NHitsWest->Draw("hist e"); }
  canv->cd(3)->SetLogz(true);
  if( mH2F_HitEpd_nmipVchkey[0]!=0 ){ mH2F_HitEpd_nmipVchkey[0]->Draw("colz"); }
  canv->cd(4)->SetLogz(true);
  if( mH2F_HitEpd_nmipVchkey[1]!=0 ){ mH2F_HitEpd_nmipVchkey[1]->Draw("colz"); }
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawEpdTacQa(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(2,2);
  canv->cd(1)->SetLogz(true);
  if( mH2F_Epd_earlywVearlye!=0 ){ mH2F_Epd_earlywVearlye->Draw("colz"); }
  canv->cd(2)->SetLogz(true);
  if( mH2F_Epd_avgwVavge!=0 ){ mH2F_Epd_avgwVavge->Draw("colz"); }
  canv->cd(3);
  if( mH1F_EpdTacDiff_Early!=0 ){ mH1F_EpdTacDiff_Early->Draw("hist e"); }
  canv->cd(4);
  if( mH1F_EpdTacDiff_Avg!=0 ){ mH1F_EpdTacDiff_Avg->Draw("hist e"); }
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawEpdTacCutQa(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(2,2);
  canv->cd(1)->SetLogz(true);
  if( mH2F_EpdCut_earlywVearlye!=0 ){ mH2F_EpdCut_earlywVearlye->Draw("colz"); }
  canv->cd(2)->SetLogz(true);
  if( mH2F_EpdCut_avgwVavge!=0 ){ mH2F_EpdCut_avgwVavge->Draw("colz"); }
  canv->cd(3);
  if( mH1F_EpdCutTacDiff_Early!=0 ){ mH1F_EpdCutTacDiff_Early->Draw("hist e"); }
  canv->cd(4);
  if( mH1F_EpdCutTacDiff_Avg!=0 ){ mH1F_EpdCutTacDiff_Avg->Draw("hist e"); }
  canv->Print(savename);
}

void StFcsRun22QaMaker::DrawEpdDepAdcQa(TCanvas* canv, const char* savename)
{
  if( mEpdAdcQaOn || (mH2F_HitPres_depVqt[0]!=0 && mH2F_HitPres_depVqt[1]!=0) ){
    canv->Clear();
    canv->Divide(5,5);
    for( UInt_t i=0; i<2; ++i ){
      for( Int_t ich=0, ipad=1; ich<mH2F_HitPres_depVqt[i]->GetEntriesFast(); ++ich,++ipad ){
	if( ipad>25 ){ ipad=1; canv->Print(savename); canv->Clear(); canv->Divide(5,5); }
	canv->cd(ipad)->SetLogz(true);
	((TH1*)mH2F_HitPres_depVqt[i]->UncheckedAt(ich))->Draw("colz");
      }
    }
    canv->Print(savename);
  }  
}

void StFcsRun22QaMaker::DrawEpdDepTacQa(TCanvas* canv, const char* savename)
{
  if( mEpdTacQaOn || (mH2F_HitPres_peakVtac[0]!=0 && mH2F_HitPres_peakVtac[1]!=0) ){
    canv->Clear();
    canv->Divide(5,5);
    for( UInt_t i=0; i<2; ++i ){
      for( Int_t ich=0, ipad=1; ich<mH2F_HitPres_peakVtac[i]->GetEntriesFast(); ++ich,++ipad ){
	if( ipad>25 ){ ipad=1; canv->Print(savename); canv->Clear(); canv->Divide(5,5); }
	canv->cd(ipad);
	((TH1*)mH2F_HitPres_peakVtac[i]->UncheckedAt(ich))->Draw("colz");
      }
    }
  }
}

void StFcsRun22QaMaker::DrawEpdTacAdcQa(TCanvas* canv, const char* savename)
{
  if( mEpdTacAdcOn || (mH2F_HitEpd_tacVadcmip[0]!=0 && mH2F_HitEpd_tacVadcmip[1]!=0) ){
    canv->Clear();
    canv->Divide(5,5);
    for( UInt_t i=0; i<2; ++i ){
      for( Int_t ich=0, ipad=1; ich<mH2F_HitEpd_tacVadcmip[i]->GetEntriesFast(); ++ich,++ipad ){
	if( ipad>25 ){ ipad=1; canv->Print(savename); canv->Clear(); canv->Divide(5,5); }
	canv->cd(ipad)->SetLogz(1);
	((TH1*)mH2F_HitEpd_tacVadcmip[i]->UncheckedAt(ich))->Draw("colz");
      }
    }
  }
}

void StFcsRun22QaMaker::DrawEpdAllQa(TCanvas* canv, const char* savename)
{
  DrawEpdHitQa(canv, savename);
  DrawEpdTacQa( canv, savename);
  DrawEpdTacCutQa( canv, savename);
  DrawEpdDepAdcQa( canv, savename);
  DrawEpdDepTacQa( canv, savename);
  DrawEpdTacAdcQa( canv, savename);
}

void StFcsRun22QaMaker::DrawFcsClusterQa(TCanvas* canv, const char* savename)
{
  for( UShort_t i=0; i<kFcsNDet; ++i ){
    canv->Clear();
    canv->Divide(3,3);
    canv->cd(1);
    mH1F_NClusters[i]->Draw("hist e");
    canv->cd(2);
    mH1F_Clu_NTowers[i]->Draw("hist e");
    canv->cd(3);
    mH1F_Clu_NNei[i]->Draw("hist e");
    canv->cd(4);
    mH1F_Clu_NPoints[i]->Draw("hist e");
    canv->cd(5);
    mH1F_Clu_En[i]->Draw("hist e");
    canv->cd(6);
    mH2F_Clu_yVx[i]->Draw("colz");
    canv->cd(7);
    mH2F_Clu_sigmaxVsigmin[i]->Draw("colz");
    canv->cd(8);
    mH1F_Clu_theta[i]->Draw("hist e");
    canv->cd(9);
    mH2F_Clu_Chi2NdfPhoton_2V1[i]->Draw("colz");
    canv->Print(savename);
  }
}

void StFcsRun22QaMaker::DrawFcsClusterPi0(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1);
  mH2F_CluHigh_angleVesum->Draw("colz");
  canv->cd(2);
  mH2F_CluHighEn_lowVhigh->Draw("colz");
  canv->cd(3);
  mH2F_CluHigh_dggVesum->Draw("colz");
  canv->cd(4);
  mH2F_CluHigh_invmassVesum->Draw("colz");
  canv->cd(5);
  mH2F_CluHigh_invmassVdgg->Draw("colz");
  canv->cd(6);
  mH2F_CluHigh_invmassVzgg->Draw("colz");
  canv->Print(savename);
}
void StFcsRun22QaMaker::DrawFcsPointQa(TCanvas* canv, const char* savename)
{
  for( UShort_t i=0; i<kFcsNDet; ++i ){
    canv->Clear();
    canv->Divide(2,2);
    canv->cd(1);
    mH1F_NPoints[i]->Draw("hist e");
    canv->cd(2);
    mH1F_Poi_En[i]->Draw("hist e");
    canv->cd(3);
    mH1F_Poi_NCluPhotons[i]->Draw("hist e");
    canv->cd(4);
    mH2F_Poi_yVx[i]->Draw("colz");
    canv->Print(savename);
  }
}

void StFcsRun22QaMaker::DrawFcsPointPi0(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1);
  mH2F_PoiHigh_angleVesum->Draw("colz");
  canv->cd(2);
  mH2F_PoiHighEn_lowVhigh->Draw("colz");
  canv->cd(3);
  mH2F_PoiHigh_dggVesum->Draw("colz");
  canv->cd(4);
  mH2F_PoiHigh_invmassVesum->Draw("colz");
  canv->cd(5);
  mH2F_PoiHigh_invmassVdgg->Draw("colz");
  canv->cd(6);
  mH2F_PoiHigh_invmassVzgg->Draw("colz");
  canv->Print(savename);
}

