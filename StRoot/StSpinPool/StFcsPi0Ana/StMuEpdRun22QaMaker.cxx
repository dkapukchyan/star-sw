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

#include "StMuEpdRun22QaMaker.h"

ClassImp(StMuEpdRun22QaMaker)


StMuEpdRun22QaMaker::StMuEpdRun22QaMaker(const char* name):StMaker(name)
{
  //mFileName = "";
  //mAllHists = new TObjArray();
  memset(mH2F_HitEpd_nmipVchkey,0,sizeof(mH2F_HitEpd_nmipVchkey));
  memset(mH2F_HitEpd_tacVadcmip,0,sizeof(mH2F_HitEpd_tacVadcmip));

  
}

StMuEpdRun22QaMaker::~StMuEpdRun22QaMaker()
{
  if( mInternalHists ){ delete mHists; }
  //delete mAllHists;
  //delete mFileOutput;
}

void StMuEpdRun22QaMaker::setHistManager( HistManager* hm )
{
  if( mInternalHists ){ delete mHists; mHists = 0; }
  mInternalHists = false;
  if( mHists!=0 ){ LOG_WARN << "StMuEpdRun22QaMaker::setHistManager() - HistManager exists and is external - no changes made" << endm; return; }
  else{ mHists = hm; }
}

UInt_t StMuEpdRun22QaMaker::LoadHists(TFile* file)
{
  if( mHists==0 ){ return 0; }
  UInt_t loaded = 0;
  loaded += mHists->AddH1F(file,mH1F_VertexEpd,"H1F_VertexEpd","Epd Vertex (z);cm",50,-200,200);
  loaded += mHists->AddH1F(file,mH1F_Epd_NHits,"H1F_Epd_NHits","Number of Hits from EPD collection;NHits",500,0,500);
  loaded += mHists->AddH1F(file,mH1F_Epd_NHits_Cut,"H1F_Epd_NHits_Cut","Number of Hits from EPD collection with nMIP>0.7;NHits",500,0,500);
  loaded += mHists->AddH1F(file,mH1F_Epd_NHitsWest,"H1F_Epd_NHitsWest","Number of Hits from EPD collection (West);NHits",300,0,300);
  loaded += mHists->AddH1F(file,mH1F_Epd_NHitsWest_Cut,"H1F_Epd_NHitsWest_Cut","Number of Hits from EPD collection (West) with nMIP>0.7;NHits",300,0,300);

  loaded += mHists->AddH2F(file,mH2F_VertexZ_vpdVepd,"H2F_VertexZ_vpdVepd","VPD vs. Epd Vertex (z);EPD (z) cm;VPD (z) cm", 50,-200,200, 50,-200,200);
  loaded += mHists->AddH2F(file,mH2F_VertexZ_zdcVepd,"H2F_VertexZ_zdcVepd","ZDC vs. Epd Vertex (z);EPD (z) cm;ZDC (z) cm", 50,-200,200, 50,-200,200);
  loaded += mHists->AddH2F(file,mH2F_VertexZ_bbcVepd,"H2F_VertexZ_bbcVepd","BBC vs. Epd Vertex (z);EPD (z) cm;BBC (z) cm", 50,-200,200, 50,-200,200);
  loaded += mHists->AddH2F(file,mH2F_VertexZ_vpdVbbc,"H2F_VertexZ_vpdVbbc","VPD vs. BBC Vertex (z);BBC (z) cm;VPD (z) cm", 50,-200,200, 50,-200,200);
  loaded += mHists->AddH2F(file,mH2F_VertexZ_vpdVzdc,"H2F_VertexZ_vpdVzdc","VPD vs. ZDC Vertex (z);ZDC (z) cm;VPD (z) cm", 50,-200,200, 50,-200,200);
  loaded += mHists->AddH2F(file,mH2F_VertexZ_zdcVbbc,"H2F_VertexZ_zdcVbbc","ZPD vs. BBC Vertex (z);BBC (z) cm;ZDC (z) cm", 50,-200,200, 50,-200,200);

  if( mEpdTacAdcOn ){
    mH2F_HitEpd_tacVadcmip[0] = new TObjArray(); //Create new array for east side
    mH2F_HitEpd_tacVadcmip[1] = new TObjArray(); //Create new array for west side
    //EPD has 372 tiles on one side
    loaded += mHists->AddH2FArr(file,mH2F_HitEpd_tacVadcmip[0],372,"H2F_HitEpd_tacVadcmip_E","Qt TAC vs. ADC/ADC_1mip;ADC/ADC_1mip;TAC", 50,0,25, 200,0,4000);
    loaded += mHists->AddH2FArr(file,mH2F_HitEpd_tacVadcmip[1],372,"H2F_HitEpd_tacVadcmip_W","Qt TAC vs. ADC/ADC_1mip;ADC/ADC_1mip;TAC", 50,0,25, 200,0,4000 );
  }

  loaded += mHists->AddH2F(file,mH2F_HitEpd_nmipVchkey[0],"H2F_HitEpd_nmipVchkey_E","NMIP values for East EPD channels by key;key ((pp-1)*31+(tt-1));nmip", 372,0,372, 100,0,25);
  loaded += mHists->AddH2F(file,mH2F_HitEpd_nmipVchkey[1],"H2F_HitEpd_nmipVchkey_W","NMIP values for West EPD channels by key;key ((pp-1)*31+(tt-1));nmip", 372,0,372, 100,0,25);
    
  loaded += mHists->AddH2F(file,mH2F_Epd_earlywVearlye,"H2F_Epd_earlywVearlye","EPD Earliest TAC West vs. East;Earliest East TAC;Earliest West TAC", 421,-10,4200, 421,-10,4200);
  loaded += mHists->AddH2F(file,mH2F_Epd_avgwVavge,"H2F_Epd_avgwVavge","EPD Average TAC West vs. East;Average East TAC;Average West TAC", 200,0,1000, 200,0,1000);
  loaded += mHists->AddH2F(file,mH2F_EpdTacDiff_avgVearly,"H2F_EpdTacDiff_avgVearly","Epd TAC difference from Average TAC vs. Early TAC;TacDiff Early;TacDiff Avg",200,-3000,3000, 200,-3000,3000);
    
  loaded += mHists->AddH2F(file,mH2F_EpdCut_earlywVearlye,"H2F_EpdCut_earlywVearlye","EPD Earliest TAC West vs. East with 1<nMIP<15;Earliest East TAC;Earliest West TAC", 300,0,4200, 300,0,4200);
  loaded += mHists->AddH2F(file,mH2F_EpdCut_avgwVavge,"H2F_EpdCut_avgwVavge","EPD Average TAC West vs. East with 1<nMIP<15;Average East TAC;Average West TAC", 300,0,4200, 300,0,4200);
  loaded += mHists->AddH2F(file,mH2F_EpdCutTacDiff_avgVearly,"H2F_EpdCutTacDiff_avgVearly","#splitline{Epd TAC difference from Average TAC vs. Early TAC}{with cuts 1<nMIP<15 & TAC>50};TacDiff Early;TacDiff Avg",200,-3000,3000, 200,-3000,3000);

  return loaded;
}

Int_t StMuEpdRun22QaMaker::Init()
{
  //LOG_INFO << "StMuEpdRun22QaMaker::Init()" << endm;
  //if( mFileName.Length() == 0){ mFileName="test.root"; } //Ensure a TFile is always created
  if( mHists==0 ){
    LOG_INFO << "StMuEpdRun22QaMaker::Init() - No HistManager specified. Creating a new one with file name StMuEpdRun22Qa.root. Potential conflicts exist if a HistManager exists with same file name" << endm;
    mHists = new HistManager();
    mInternalHists = true;
    mHists->InitFile("StMuEpdRun22Qa.root");
  }

  //if( mFileName.Length() == 0){ mFileName="test.root"; } //Ensure a TFile is always created
  //mFileOutput = new TFile(mFileName.Data(), "RECREATE");
  UInt_t totalhists = this->LoadHists(0); //This is total of histograms loaded from a file not created. Don't use mFileOutput as you are not trying to load from #mFileOutput
  mHists->SetOwner(kTRUE);
  LOG_INFO << "StMuEpdRun22QaMaker::Init() - Loaded " << totalhists << " histograms" << endm;
  return kStOk;
}

Int_t StMuEpdRun22QaMaker::InitRun(int runnumber)
{
  //LOG_INFO << "StMuEpdRun22QaMaker::InitRun()" << endm;
  /*
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  //mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
    LOG_ERROR << "StMuEpdRun22QaMaker::InitRun Failed to get StFcsDbMaker" << endm;
    return kStFatal;
    }*/
  // if( !mEpdGeo ){ mEpdGeo = new StEpdGeom(); }
  // else{
  //   LOG_ERROR << "StMuEpdRun22QaMaker::InitRun - StEpdGeom Exists!" << endm;
  //   return kStFatal;
  // }
  return kStOk;
}

Int_t StMuEpdRun22QaMaker::Make()
{
  //LOG_INFO << "StMuEpdRun22QaMaker::Make()" << endm;
  mMuDstMkr = (StMuDstMaker*)GetInputDS("MuDst");
  if( mMuDstMkr==0 ){ LOG_ERROR <<"StMuEpdRun22QaMaker::Make - !MuDstMkr" <<endm; return kStErr; }
  mMuDst = mMuDstMkr->muDst();
  if( mMuDst==0 ){ LOG_ERROR << "StMuEpdRun22QaMaker::Make - !MuDst" << endm; return kStErr; }
  mMuEvent = mMuDst->event();
  if( mMuEvent==0 ){ LOG_ERROR <<"StMuEpdRun22QaMaker::Make - !MuEvent" <<endm; return kStErr; }
  mTrigData = mMuEvent->triggerData();
  if( mTrigData==0 ){ LOG_ERROR <<"StMuEpdRun22QaMaker::Make - !TrigData" <<endm; return kStErr; }
  mRunInfo = &(mMuEvent->runInfo());
  if( mRunInfo==0 ){ LOG_ERROR <<"StMuEpdRun22QaMaker::Make - !RunInfo" <<endm; return kStErr; }
  
  mMuEpdHits = 0;
  mEpdColl = 0;
  mMuEpdHits = mMuDst->epdHits();
  if( mMuEpdHits!=0 ){ if( mMuEpdHits->GetEntriesFast()==0 ){mMuEpdHits=0;} }//If mMuEpdHits is not zero but has no hits set it to zero so rest of code processes from StEpdHitMaker
  if( mMuEpdHits==0 ){ LOG_INFO << "StMuEpdRun22QaMaker::Make - No MuEPD hits" << endm;
    mEpdHitMkr = (StEpdHitMaker*)GetMaker("epdHit");
    if( mEpdHitMkr==0 ){ LOG_WARN << "StMuEpdRun22QaMaker::Make - No StEpdHitMaker(\"epdHit\")" << endm; }
    else{ mEpdColl = mEpdHitMkr->GetEpdCollection(); }
    if( mEpdColl==0 ){ LOG_WARN << "StMuEpdRun22QaMaker::FillFcsInfo - No Epd hit information found" << endm; mEpdHitMkr=0; }//Set the hit maker back to zero so it can be used as a check that the epd collection doesn't exist
  }
  Int_t epdstatus = this->FillEpdInfo();
  switch( epdstatus ){
  case kStEOF: return kStEOF;
  case kStErr: return kStErr;
  case kStFatal: return kStFatal;
  case kStSkip: return kStSkip;
  case kStStop: return kStStop;
  }

  //Do this after EPD vertex was filled
  Double_t vpdz = 0;
  if( mMuDst->btofHeader() ){ vpdz = mMuDst->btofHeader()->vpdVz(); }
  Double_t zdcz = mTrigData->zdcVertexZ();
  const float bbcTdiff = mTrigData->bbcTimeDifference() - 4096; //subtract 4096 since 0 means bad event and distribution is Gaussian around 4096
  Double_t bbcz = 0;
  if( fabs(bbcTdiff)>1.e-6 ){ bbcz = bbcTdiff * -0.2475; }

  mH2F_VertexZ_vpdVepd->Fill(mEpdVertex,vpdz);
  mH2F_VertexZ_zdcVepd->Fill(mEpdVertex,zdcz);
  mH2F_VertexZ_bbcVepd->Fill(mEpdVertex,bbcz);
  mH2F_VertexZ_vpdVbbc->Fill(bbcz,vpdz);
  mH2F_VertexZ_vpdVzdc->Fill(zdcz,vpdz);
  mH2F_VertexZ_zdcVbbc->Fill(bbcz,zdcz);
  
  //std::cout << "Finished Make" << std::endl;
  if( epdstatus==kStWarn ){ return kStWarn; } //Now check if either returned a warning and if so returning warning
  else{ return kStOk; } //Both were ok
}

Int_t StMuEpdRun22QaMaker::FillEpdInfo()
{
  if( mMuEpdHits==0 && mEpdColl==0 ){
    LOG_WARN << "StMuEpdRun22QaMaker::FillEpdInfo() has no epd hits" << endm;
    return kStWarn;
  }
  unsigned int nepdhits = 0;
  StSPtrVecEpdHit* epdhits = 0;
  if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
  else if( mEpdColl!=0 ){
    epdhits = &(mEpdColl->epdHits());
    nepdhits = epdhits->size();
  }
  else{ LOG_ERROR << "StMuEpdRun22QaMaker::FillEpdinfo() - If you see this error then there is a bug that is setting EPD hits improperly" << endm; return kStErr; }
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
  int nhitscut = 0;
  int nhitswestcut = 0;
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
    if( nmip>0.7 ){ ++nhitscut; }
    if( ew==-1 ){ //This is east side
      if( mEpdTacAdcOn ){ ((TH1*)mH2F_HitEpd_tacVadcmip[0]->UncheckedAt( (epdpp-1)*31+(epdtt-1) ))->Fill(nmip,tac); }
      mH2F_HitEpd_nmipVchkey[0]->Fill((epdpp-1)*31+(epdtt-1),nmip);
      sumtace += tac;
      ++navgtace;
      if( tac > earliesttace ){ earliesttace = tac; }
      if( 1<nmip && nmip<15 && tac>50){
	cut_sumtace += tac;
	++cut_navgtace;
	if( tac > cut_earliesttace ){ cut_earliesttace = tac; }
      }
    }
    if( ew==1 ){ //This is west side
      if( nmip>0.7 ){ ++nhitswestcut; }
      if( mEpdTacAdcOn ){ ((TH1*)mH2F_HitEpd_tacVadcmip[1]->UncheckedAt( (epdpp-1)*31+(epdtt-1) ))->Fill(nmip,tac); }
      mH2F_HitEpd_nmipVchkey[1]->Fill((epdpp-1)*31+(epdtt-1),nmip);
      sumtacw += tac;
      ++navgtacw;
      ++nhitwest;
      if( tac > earliesttacw ){ earliesttacw = tac; }
      if( 1<nmip && nmip<15 && tac>50 ){
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

  if( mH1F_Epd_NHits_Cut!=0 ){ mH1F_Epd_NHits_Cut->Fill(nhitscut); }
  if( mH1F_Epd_NHitsWest!=0 ){ mH1F_Epd_NHitsWest->Fill(nhitwest); }
  if( mH1F_Epd_NHitsWest_Cut!=0 ){ mH1F_Epd_NHitsWest_Cut->Fill(nhitswestcut); }

  //std::cout << "|earlye:"<<earliesttace << "|earlyw:"<<earliesttacw << "|avge:"<<avgtace << "|avgw:"<<avgtacw << std::endl;
  if( mH2F_Epd_avgwVavge!=0 )    { mH2F_Epd_avgwVavge->Fill( avgtace,avgtacw  ); }
  if( mH2F_Epd_earlywVearlye!=0 ){ mH2F_Epd_earlywVearlye->Fill(earliesttace,earliesttacw); }
  //if( mH1F_EpdTacDiff_Early!=0 ) { mH1F_EpdTacDiff_Early->Fill( earliesttacw - earliesttace ); }
  //if( mH1F_EpdTacDiff_Avg!=0 )   { mH1F_EpdTacDiff_Avg->Fill( avgtacw - avgtace ); }
  if( mH2F_EpdTacDiff_avgVearly!=0 ){ mH2F_EpdTacDiff_avgVearly->Fill( earliesttacw-earliesttace, avgtacw - avgtace ); }

  //std::cout << "|cut_earlye:"<<cut_earliesttace << "|cut_earlyw:"<<cut_earliesttacw << "|cut_avge:"<<cut_avgtace << "|cut_avgw:"<<cut_avgtacw << std::endl;
  if( mH2F_EpdCut_avgwVavge!=0 )    { mH2F_EpdCut_avgwVavge->Fill( cut_avgtace,cut_avgtacw  ); }
  if( mH2F_EpdCut_earlywVearlye!=0 ){ mH2F_EpdCut_earlywVearlye->Fill(cut_earliesttace,cut_earliesttacw); }
  //if( mH1F_EpdCutTacDiff_Early!=0 ) { mH1F_EpdCutTacDiff_Early->Fill( cut_earliesttacw - cut_earliesttace ); }
  //if( mH1F_EpdCutTacDiff_Avg!=0 )   { mH1F_EpdCutTacDiff_Avg->Fill( cut_avgtacw - cut_avgtace ); }
  if( mH2F_EpdCutTacDiff_avgVearly!=0 ){ mH2F_EpdCutTacDiff_avgVearly->Fill( cut_earliesttacw-cut_earliesttace, cut_avgtacw-cut_avgtace ); }

  //For vertex only fill if the average TAC is >50 which is determined by eye for Run 22
  if( cut_avgtace>50 && cut_avgtace>50 ){
    mEpdVertex = (cut_avgtacw-cut_avgtace) * -0.2475;
    mH1F_VertexEpd->Fill( mEpdVertex );
  }
  else{ mEpdVertex = 0; }
  //@[Julye 17, 2024] > Using the same logic as BBC since haven't looked at EpdTacDiff histograms, also because EPD has same 15.6ps/TAC as BBC.
  //@[July 21, 2024]>After looking at TAC differences it seems average with no cuts gives actual results and ADC_Nmip cut tac values do not.
  //@[July 23, 2024]>Looking at TAC diff don't need to subtract 4096 like we do for BBC
  return kStOk;
}

Int_t StMuEpdRun22QaMaker::Finish()
{
  //LOG_INFO << "StMuEpdRun22QaMaker::Finish()" << endm;
  TFile* file = mHists->InitFile(); //No arguments just returns the file pointer in #HistManager
  if( file!=0 ){
    LOG_INFO << "StMuFcsRun22QaMaker::Finish() - Writing to file:" << file->GetName() << endm;
    file->cd();
    mHists->Write();
    return kStOk;
  }
  /*
  if( mFileOutput!=0 ){
    mFileOutput->cd();
    mAllHists->Write();
    return kStOk;
  }
  */
  else{
    LOG_WARN << "StMuEpdRun22QaMaker::Finish() - No file created because pointer is null" << endm;
    return kStWarn;
  }
}

void StMuEpdRun22QaMaker::DrawEpdAllQa(TCanvas* canv, const char* savename)
{
  DrawVertex(canv, savename);
  DrawEpdHitQa(canv, savename);
  DrawEpdTacQa( canv, savename);
  DrawEpdTacCutQa( canv, savename);
  DrawEpdTacAdcQa( canv, savename);
}

void StMuEpdRun22QaMaker::DrawVertex(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1);
  if( mH2F_VertexZ_vpdVepd!=0 ){ mH2F_VertexZ_vpdVepd->Draw("colz"); }
  canv->cd(2);
  if( mH2F_VertexZ_zdcVepd!=0 ){ mH2F_VertexZ_zdcVepd->Draw("colz"); }
  canv->cd(3);
  if( mH2F_VertexZ_bbcVepd!=0 ){ mH2F_VertexZ_bbcVepd->Draw("colz"); }
  canv->cd(4);
  if( mH2F_VertexZ_vpdVbbc!=0 ){ mH2F_VertexZ_vpdVbbc->Draw("colz"); }
  canv->cd(5);
  if( mH2F_VertexZ_vpdVzdc!=0 ){ mH2F_VertexZ_vpdVzdc->Draw("colz"); }
  canv->cd(6);
  if( mH2F_VertexZ_zdcVbbc!=0 ){ mH2F_VertexZ_zdcVbbc->Draw("colz"); }
  canv->Print(savename);
}

void StMuEpdRun22QaMaker::DrawEpdHitQa(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1);
  if( mH1F_Epd_NHits!=0 ){ mH1F_Epd_NHits->Draw("hist e"); }
  canv->cd(2);
  if( mH1F_Epd_NHits_Cut!=0 ){ mH1F_Epd_NHits_Cut->Draw("hist e"); }
  canv->cd(3)->SetLogz(true);
  if( mH2F_HitEpd_nmipVchkey[0]!=0 ){ mH2F_HitEpd_nmipVchkey[0]->Draw("colz"); }
  canv->cd(4);
  if( mH1F_Epd_NHitsWest!=0 ){ mH1F_Epd_NHitsWest->Draw("hist e"); }
  canv->cd(5);
  if( mH1F_Epd_NHitsWest_Cut!=0 ){ mH1F_Epd_NHitsWest_Cut->Draw("hist e"); }  
  canv->cd(6)->SetLogz(true);
  if( mH2F_HitEpd_nmipVchkey[1]!=0 ){ mH2F_HitEpd_nmipVchkey[1]->Draw("colz"); }
  canv->Print(savename);
}

void StMuEpdRun22QaMaker::DrawEpdTacQa(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1)->SetLogz(true);
  if( mH2F_Epd_earlywVearlye!=0 ){ mH2F_Epd_earlywVearlye->Draw("colz"); }
  canv->cd(2)->SetLogz(true);
  if( mH2F_Epd_avgwVavge!=0 ){ mH2F_Epd_avgwVavge->Draw("colz"); }
  if( mH2F_EpdTacDiff_avgVearly!=0 ){
    canv->cd(3);
    mH2F_EpdTacDiff_avgVearly->Draw("colz");
    TH1D* earlytac = ((TH2*)mH2F_EpdTacDiff_avgVearly)->ProjectionX("H1F_EpdTacDiff_early");
    earlytac->SetTitle("Epd TAC difference from Earliest TAC;TacDiff Early");
    TH1D* avgtac   = ((TH2*)mH2F_EpdTacDiff_avgVearly)->ProjectionY("H1F_EpdTacDiff_avg");
    avgtac->SetTitle("Epd TAC difference from Average TAC;TacDiff Avg");
    canv->cd(4);
    earlytac->Draw("hist e");
    canv->cd(5);
    avgtac->Draw("hist e");    
  }
  canv->Print(savename);
}

void StMuEpdRun22QaMaker::DrawEpdTacCutQa(TCanvas* canv, const char* savename)
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1)->SetLogz(true);
  if( mH2F_EpdCut_earlywVearlye!=0 ){ mH2F_EpdCut_earlywVearlye->Draw("colz"); }
  canv->cd(2)->SetLogz(true);
  if( mH2F_EpdCut_avgwVavge!=0 ){ mH2F_EpdCut_avgwVavge->Draw("colz"); }
  if( mH2F_EpdCutTacDiff_avgVearly!=0 ){
    canv->cd(3);
    mH2F_EpdCutTacDiff_avgVearly->Draw("colz");
    TH1D* earlytac = ((TH2*)mH2F_EpdCutTacDiff_avgVearly)->ProjectionX("H1F_EpdCutTacDiff_early");
    earlytac->SetTitle("Epd TAC difference from Earliest TAC and 1<nMIP<15 and TAC>50;TacDiff Early");
    TH1D* avgtac   = ((TH2*)mH2F_EpdCutTacDiff_avgVearly)->ProjectionY("H1F_EpdCutTacDiff_avg");
    avgtac->SetTitle("Epd TAC difference from Average TAC and 1<nMIP<15 and TAC>50;TacDiff Avg");
    canv->cd(4);
    earlytac->Draw("hist e");
    canv->cd(5);
    avgtac->Draw("hist e");
  }
  canv->cd(6);
  mH1F_VertexEpd->Draw("hist e");
  canv->Print(savename);
}

/*
void StMuEpdRun22QaMaker::DrawEpdDepAdcQa(TCanvas* canv, const char* savename)
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

void StMuEpdRun22QaMaker::DrawEpdDepTacQa(TCanvas* canv, const char* savename)
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
*/
void StMuEpdRun22QaMaker::DrawEpdTacAdcQa(TCanvas* canv, const char* savename)
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



