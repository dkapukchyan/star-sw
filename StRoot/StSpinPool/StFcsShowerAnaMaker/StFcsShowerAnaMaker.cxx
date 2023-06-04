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
    StThreeVectorD xyzoff = mFcsDb->getDetectorOffset(det);
    double alpha = mFcsDb->getDetectorAngle(det);
    std::cout << "|det:"<<det << "|nclus:"<<nc << "|npoints:"<<np << "|ShowerZmax:"<<mFcsDb->getShowerMaxZ(det) << "|det:"<<det << "|alpha:"<<alpha << "|xoff:"<<xyzoff.x() << "|yoff:"<<xyzoff.y() << "|zoff:"<<xyzoff.z() << std::endl;
    for( int ipoint=0; ipoint<np; ++ipoint ){
      StFcsPoint* point=points[ipoint];
      if( point->energy() > mEnCut ){
	float xfull = point->x();
	float xwhole = floor(xfull);
	mH1F_PointXLocal->Fill( xfull-xwhole ); //Only store fractional part
	float yfull = point->y();
	float ywhole = floor(yfull);
	mH1F_PointYLocal->Fill( yfull-ywhole ); //Only store fractional part

	StThreeVectorF pointxyz = point->xyz();
	std::cout << "|ipoint:"<<ipoint << "|E:"<<point->energy() << "|x:"<<pointxyz.x() << "|y:"<<pointxyz.y() << "|z:"<<pointxyz.z() << std::endl;

	StFcsCluster* pointclus = point->cluster();
	g2t_track_st* g2ttrk = 0;
	St_g2t_track* trackTable = static_cast<St_g2t_track*>(GetDataSet("g2t_track"));
	if( !trackTable ){ std::cout<< "g2t_track Table not found" << std::endl; continue; }
	else{
	  const int nTrk = trackTable->GetNRows();
	  std::cout << "g2t_track table has "<< nTrk << " tracks" << std::endl;
	  if( nTrk>0 ){
	    g2ttrk = trackTable->GetTable();
	    if( !g2ttrk) { std::cout << " g2t_track GetTable failed" << std::endl; continue; }
	  }
	}
	if( g2ttrk ){
	  float frac=0;
	  int ntrk=0;
	  //const g2t_track_st* parenttrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  //std::cout << "|parenttrk|Id:"<<parenttrk->id << "|Pid:"<<parenttrk->ge_pid << "|E:"<<parenttrk->e << "|eta:"<<parenttrk->eta << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	  const g2t_track_st* primtrk = mFcsDb->getPrimaryG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  double pt2 = sqrt( primtrk->p[0]*primtrk->p[0] + primtrk->p[1]*primtrk->p[1] );
	  double phix = acos(primtrk->p[0]/primtrk->pt);
	  double phiy = asin(primtrk->p[1]/primtrk->pt);
	  double eta2 = asinh(primtrk->p[2]/primtrk->pt);
	  double phi = 0;
	  if( primtrk->p[0]>0 && primtrk->p[1]>0 ){ phi=phix; }            //first quadrant angles match
	  if( primtrk->p[0]<0 && primtrk->p[1]>0 ){ phi=phix; }            //second quadrant take arccos
	  if( primtrk->p[0]<0 && primtrk->p[1]<0 ){ phi=phix+3.1415/2.0; } //third quadrant take arccos+90
	  if( primtrk->p[0]>0 && primtrk->p[1]<0 ){ phi=phiy; }            //fourth quadrant take arcsin
	  //double theta = 2.0*atan2(exp(-1.0*primtrk->eta),1);
	  //double radius = mFcsDb->getShowerMaxZ(det)/cos(theta);
	  //double xprim = radius*cos(phi)*sin(theta);
	  //double yprim = radius*sin(phi)*sin(theta);
	  //double Doffset = sqrt(xyzoff.x()*xyzoff.x()+xyzoff.z()*xyzoff.z());
	  //double xprim = Doffset*cos(theta)/cos(theta-alpha);
	  //double zprim = Doffset*sin(theta)/cos(theta-alpha);
	  //std::cout << "|primtrk|Id:"<<primtrk->id << "|Pid:"<<primtrk->ge_pid << "|E:"<<primtrk->e << "|px:"<<primtrk->p[0] << "|py:"<<primtrk->p[1] << "|pz:"<<primtrk->p[2] << "|pt:"<<primtrk->pt << "|ptot:"<<primtrk->ptot << "|eta:"<<primtrk->eta << "|phix:"<<phix*(180.0/3.1415) << "|phiy:"<<phiy*(180.0/3.1415) << "|phi:"<<phi*180.0/3.1415 << "|theta:"<<theta*(180.0/3.1415) << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	  //std::cout << "|primtrk|Id:"<<primtrk->id << "|Pid:"<<primtrk->ge_pid << "|E:"<<primtrk->e << "|px:"<<primtrk->p[0] << "|py:"<<primtrk->p[1] << "|pz:"<<primtrk->p[2] << "|pt:"<<primtrk->pt << "|ptot:"<<primtrk->ptot << "|eta:"<<primtrk->eta << "|phi:"<<phi*180.0/3.1415 << "|theta:"<<theta*180.0/3.1415 << "|r:"<<radius << "|x:"<<xprim << "|y:"<<yprim << "|z:"<<pointxyz.z() << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	  std::cout << "|primtrk|Id:"<<primtrk->id << "|Pid:"<<primtrk->ge_pid << "|E:"<<primtrk->e << "|px:"<<primtrk->p[0] << "|py:"<<primtrk->p[1] << "|pz:"<<primtrk->p[2] << "|pt:"<<primtrk->pt << "|pt2:"<<pt2 << "|ptot:"<<primtrk->ptot << "|eta:"<<primtrk->eta << "|eta2:"<<eta2 << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	}
      }
    }
  }
  return kStOk;
}

