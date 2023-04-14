#include "StFcsShowerAnaMaker.h"

#include "StEnumerations.h"
#include "StMessMgr.h"
#include "Stypes.h"
#include "StEventTypes.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StThreeVectorF.hh"

#include "StFcsDbMaker/StFcsDb.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsCollection.h"
#include "StFcsCluster.h"
#include "StFcsPoint.h"


StFcsShowerAnaMaker::StFcsShowerAnaMaker(const char* name):StMaker(name)
{
}

StFcsShowerAnaMaker::~StFcsShowerAnaMaker()
{
  CleanHists();
}

Int_t StFcsShowerAnaMaker::Init()
{
  if( !mFcsDb ){ mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb")); }
  //mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
    LOG_ERROR << "StFcsShowerAnaMaker::InitRun Failed to get StFcsDb" << endm;
    return kStFatal;
  }
  InitHists();
  return kStOk;
}

Int_t StFcsShowerAnaMaker::Finish()
{
  if (mFileName.Length() == 0) return kStOk;
  TFile* outfile = new TFile(mFileName.Data(), "RECREATE");

  outfile->cd();
  WriteHists();
  return kStOk;
}

void StFcsShowerAnaMaker::InitHists()
{
  mH1F_PointXLocal = new TH1F("H1F_PointXLocal","",100,0,1);
  mH1F_PointYLocal = new TH1F("H1F_PointYLocal","",100,0,1);
  
  mH1F_ClusSigMax = new TH1F("H1F_ClusSigMax","",200,0,2);
  mH1F_ClusSigMin = new TH1F("H1F_ClusSigMin","",200,0,2);
  
  mH2F_ClusSigMaxEn = new TH2F("H2F_ClusSigMaxEn","",100,0,100,200,0,2);
  mH2F_ClusSigMinEn = new TH2F("H2F_ClusSigMinEn","",100,0,100,200,0,2);
  /*
  else if (mFile!=0 ){
    mH1F_PointXLocal = (TH1F*)mFile->Get("H1F_PointXLocal");
    mH1F_PointYLocal = (TH1F*)mFile->Get("H1F_PointYLocal");

    mH1F_mClusSigMax = (TH1F*)mFile->Get("H1F_ClusSigMax");
    mH1F_mClusSigMin = (TH1F*)mFile->Get("H1F_ClusSigMin");

    mH2F_mClusSigMaxEn = (TH1F*)mFile->Get("H2F_ClusSigMaxEn");
    mH2F_mClusSigMinEn = (TH1F*)mFile->Get("H2F_ClusSigMinEn");
  }
  else{ std::cout << "ERROR - StFcsShowerAnaMaker::InitHists():Unable to create or load histograms" << std::endl; }
  */
}

void StFcsShowerAnaMaker::CleanHists()
{
  delete mH1F_PointXLocal;
  delete mH1F_PointYLocal;
  
  delete mH1F_ClusSigMax;
  delete mH1F_ClusSigMin;
  
  delete mH2F_ClusSigMaxEn;
  delete mH2F_ClusSigMinEn;
}

void StFcsShowerAnaMaker::WriteHists()
{
  mH1F_PointXLocal->Write();
  mH1F_PointYLocal->Write();
  
  mH1F_ClusSigMax->Write();
  mH1F_ClusSigMin->Write();
  
  mH2F_ClusSigMaxEn->Write();
  mH2F_ClusSigMinEn->Write();
}

Int_t StFcsShowerAnaMaker::Make()
{
  StEvent* event = (StEvent*)GetInputDS("StEvent");
  if (!event) {
    LOG_ERROR << "StFcsShowerAnaMaker::Make did not find StEvent" << endm;
    return kStErr;
  }
  mFcsColl = event->fcsCollection();
  if (!mFcsColl) {
    LOG_ERROR << "StFcsShowerAnaMaker::Make did not find StEvent->StFcsCollection" << endm;
    return kStErr;
  }
  
  for( int det=0; det<kFcsNDet; det++ ){
    
    StSPtrVecFcsCluster& clusters = mFcsColl->clusters(det);
    int nc = mFcsColl->numberOfClusters(det);
    for( int iclus=0; iclus<nc; ++iclus ){
      StFcsCluster* cluster = clusters[iclus];
      if( cluster->energy() > mEnCut ){
	mH1F_ClusSigMax->Fill(cluster->sigmaMax());
	mH1F_ClusSigMin->Fill(cluster->sigmaMin());
	mH2F_ClusSigMaxEn->Fill(cluster->energy(),cluster->sigmaMax());
	mH2F_ClusSigMinEn->Fill(cluster->energy(),cluster->sigmaMin());
      }
    }

    StSPtrVecFcsPoint& points = mFcsColl->points(det);
    int np = mFcsColl->numberOfPoints(det);
    for( int ipoint=0; ipoint<np; ++ipoint ){
      StFcsPoint* point=points[ipoint];
      if( point->energy() > mEnCut ){
	float xfull = point->x();
	float xwhole = floor(xfull);
	mH1F_PointXLocal->Fill( xfull-xwhole ); //Only store fractional part
	float yfull = point->y();
	float ywhole = floor(yfull);
	mH1F_PointYLocal->Fill( yfull-ywhole ); //Only store fractional part
      }
    }
  }
  return kStOk;
}
