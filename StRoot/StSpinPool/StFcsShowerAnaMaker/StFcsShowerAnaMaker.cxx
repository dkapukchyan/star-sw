#include "StFcsShowerAnaMaker.h"

#include "StEnumerations.h"
#include "StMessMgr.h"
#include "Stypes.h"
#include "StEventTypes.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StThreeVectorF.hh"
#include "tables/St_g2t_vertex_Table.h"

#include "StFcsDbMaker/StFcsDb.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsCollection.h"
#include "StFcsCluster.h"
#include "StFcsPoint.h"

//#include "/star/u/dkap7827/Tools2/Tools/MyTools/inc/Ctools.h"
//#include "/star/u/dkap7827/Tools2/Tools/MyTools/inc/Rtools.h"
ClassImp(StFcsShowerAnaMaker)

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

  if( mDataTree==0 ){
    mDataTree = new Rtools::ClonesArrTree("ShowerAnaTree","");
    mDataTree->AddArr("FcsHit",    "StFcsPicoHit");
    mDataTree->AddArr("FcsCluster","StFcsPicoCluster");
    mDataTree->AddArr("FcsPoint",  "StFcsPicoPoint");
    mDataTree->AddArr("G2tPrim",   "StPicoG2tTrack");
  }
  else{ std::cout << "WARNING StFcsShowerAnaMaker::Init - DataTree Exists" << std::endl; }
  
  return kStOk;
}

Int_t StFcsShowerAnaMaker::Finish()
{
  if (mFileName.Length() == 0) return kStOk;
  TFile* outfile = new TFile(mFileName.Data(), "RECREATE");

  outfile->cd();
  WriteHists();
  if( mDataTree ){ mDataTree->Write(); }
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

  int totalhits = 0;   //Running counter for all hits for all detectors
  int totalclus = 0;   //Running counter for all clusters for all detectors
  int totalpoints = 0; //Running counter for all points for all detectors
  
  TClonesArray* mHitsArr = mDataTree->GetClassArr("StFcsPicoHit");
  TClonesArray* mClusArr = mDataTree->GetClassArr("StFcsPicoCluster");
  TClonesArray* mPointArr = mDataTree->GetClassArr("StFcsPicoPoint");
  TClonesArray* mG2tPrimArr = mDataTree->GetClassArr("StPicoG2tTrack");
  //std::cout << "|HitsArr:"<<mHitsArr << "|ClusArr:"<<mClusArr << "|PointArr:"<<mPointArr << "|G2tArr:"<<mG2tPrimArr << std::endl;
  for( int det=0; det<kFcsNDet; det++ ){
    StSPtrVecFcsHit& hits = mFcsColl->hits(det);
    int nh = mFcsColl->numberOfHits(det);
    //std::cout << "|det:"<<det << "|nh:"<<nh << std::endl;
    for( int ihit=0; ihit<nh; ++ihit ){
      StFcsHit* hit=hits[ihit];
      StFcsPicoHit* picohit = (StFcsPicoHit*)mHitsArr->ConstructedAt(totalhits++);
      //std::cout << "|det:"<<det<<"|ihit:"<<ihit << "|hit:"<<hit << "|picohit:"<<picohit << std::endl;
      picohit->mDetId  = hit->detectorId();
      picohit->mChId   = hit->id();
      picohit->mAdcSum = hit->adcSum();
      picohit->mEnergy = hit->energy();
      
      if( det<=kFcsHcalSouthDetId ){ //ECAL and HCAL
	StThreeVectorF xyz = mFcsDb->getStarXYZ(hit);
	picohit->mXstar = xyz.x();
	picohit->mYstar = xyz.y();
	picohit->mZstar = xyz.z();
      }
      else if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){//EPD as Pres
	picohit->mZstar = 375.0; //Taken from StFcsEventDisplay for zepd this is in cm
      }
      //std::cout << "|i:"<<ihit<<"|th:"<<totalhits<<"("<<picohit->mXstar<<","<<picohit->mYstar<<","<<picohit->mEnergy<<")"<<std::endl;
    }

    StSPtrVecFcsCluster& clusters = mFcsColl->clusters(det);
    int nc = mFcsColl->numberOfClusters(det);
    //std::cout << "|det:"<<det << "|nc:"<<nc << std::endl;
    for( int iclus=0; iclus<nc; ++iclus ){
      StFcsCluster* cluster = clusters[iclus];
      StFcsPicoCluster* picoclus = (StFcsPicoCluster*)mClusArr->ConstructedAt(totalclus++);
      picoclus->mId             = cluster->id();
      picoclus->mDetId          = cluster->detectorId();
      picoclus->mCategory       = cluster->category();
      picoclus->mNTowers        = cluster->nTowers();
      picoclus->mNNeighbor      = cluster->nNeighbor();
      picoclus->mNPoints        = cluster->nPoints();
      picoclus->mEnergy         = cluster->energy();
      picoclus->mX              = cluster->x();
      picoclus->mY              = cluster->y();
      picoclus->mSigmaMin       = cluster->sigmaMin();
      picoclus->mSigmaMax       = cluster->sigmaMax();
      picoclus->mTheta          = cluster->theta();
      picoclus->mChi2Ndf1Photon = cluster->chi2Ndf1Photon();
      picoclus->mChi2Ndf2Phoron = cluster->chi2Ndf2Photon();
      //Lorentz 4 momentum of cluster
      StLorentzVectorD clusp = cluster->fourMomentum();
      picoclus->mLorentzX = clusp.x();
      picoclus->mLorentzY = clusp.y();
      picoclus->mLorentzZ = clusp.z();
      picoclus->mLorentzE = clusp.e();
      
      if( cluster->energy() > mEnCut ){
	mH1F_ClusSigMax->Fill(cluster->sigmaMax());
	mH1F_ClusSigMin->Fill(cluster->sigmaMin());
	mH2F_ClusSigMaxEn->Fill(cluster->energy(),cluster->sigmaMax());
	mH2F_ClusSigMinEn->Fill(cluster->energy(),cluster->sigmaMin());
      }
    }
    
    StSPtrVecFcsPoint& points = mFcsColl->points(det);
    int np = mFcsColl->numberOfPoints(det);
    //std::cout << "|det:"<<det << "|np:"<<np << std::endl;
    for( int ipoint=0; ipoint<np; ++ipoint ){
      StFcsPoint* point=points[ipoint];
      //std::cout << "|det:"<<det<<"|ipoint:"<<ipoint << "|point:"<<point << std::endl;
      if( point->energy() > mEnCut ){
	float xfull = point->x();
	float xwhole = floor(xfull);
	mH1F_PointXLocal->Fill( xfull-xwhole ); //Only store fractional part
	float yfull = point->y();
	float ywhole = floor(yfull);
	mH1F_PointYLocal->Fill( yfull-ywhole ); //Only store fractional part

	StFcsCluster* pointclus = point->cluster();
	g2t_track_st* g2ttrk = 0;
	g2t_vertex_st* g2tvert = 0;
	St_g2t_track* trackTable = static_cast<St_g2t_track*>(GetDataSet("g2t_track"));
	St_g2t_vertex* vertexTable = static_cast<St_g2t_vertex*>(GetDataSet("g2t_vertex"));
	double speedperc = 0;
	if( !trackTable ){ std::cout<< "g2t_track Table not found" << std::endl; continue; }
	else{
	  const int nTrk = trackTable->GetNRows();
	  std::cout << "g2t_track table has "<< nTrk << " tracks" << std::endl;
	  if( nTrk>0 ){
	    g2ttrk = trackTable->GetTable();
	    if( !g2ttrk) { std::cout << " g2t_track GetTable failed" << std::endl; continue; }
	  }
	}
	if( !vertexTable ){ std::cout<< "g2t_vertex Table not found" << std::endl; continue; }
	else{
	  const int nVertex = vertexTable->GetNRows();
	  std::cout << "g2t_vertex table has "<< nVertex << " vertices" << std::endl;
	  if( nVertex>0 ){
	    g2tvert = vertexTable->GetTable();
	    if( !g2tvert) { std::cout << " g2t_vertex GetTable failed" << std::endl; continue; }
	  }
	}
	if( g2ttrk && g2tvert ){
	  StFcsPicoPoint* picopoint = (StFcsPicoPoint*)mPointArr->ConstructedAt(totalpoints);
	  picopoint->mNS                    = point->detectorId();
	  picopoint->mEnergy                = point->energy();
	  picopoint->mXlocal                = point->x();
	  picopoint->mYlocal                = point->y();
	  picopoint->mNParentClusterPhotons = point->nParentClusterPhotons();
	  //STAR xyx
	  StThreeVectorF pointxyz = point->xyz();
	  //std::cout << "|ipoint:"<<ipoint << "|E:"<<point->energy() << "|x:"<<pointxyz.x() << "|y:"<<pointxyz.y() << "|z:"<<pointxyz.z() << std::endl;
	  picopoint->mXstar = pointxyz.x();
	  picopoint->mYstar = pointxyz.y();
	  picopoint->mZstar = pointxyz.z();
	  //Lorentz 4 momentum of point
	  StLorentzVectorD pointp = point->fourMomentum();
	  picopoint->mLorentzX = pointp.x();
	  picopoint->mLorentzY = pointp.y();
	  picopoint->mLorentzZ = pointp.z();
	  picopoint->mLorentzE = pointp.e();

	  float frac=0;
	  int ntrk=0;
	  //const g2t_track_st* parenttrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  //std::cout << "|parenttrk|Id:"<<parenttrk->id << "|Pid:"<<parenttrk->ge_pid << "|E:"<<parenttrk->e << "|eta:"<<parenttrk->eta << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	  //const g2t_track_st* primtrk = mFcsDb->getPrimaryG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  const g2t_track_st* primtrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  StPicoG2tTrack* picotrk = (StPicoG2tTrack*)mG2tPrimArr->ConstructedAt(totalpoints++);
	  StThreeVectorD projshowerxyz1 = mFcsDb->projectTrackToEcal(primtrk);
	  StThreeVectorD projshowerxyz2 = mFcsDb->projectTrackToEcalSMax(primtrk);
	  //std::cout << "|picotrk:"<<picotrk << "|primtrk:"<<primtrk << std::endl;
	  picotrk->id = primtrk->id;
	  picotrk->ge_pid = primtrk->ge_pid;
	  picotrk->mPx = primtrk->p[0];
	  picotrk->mPy = primtrk->p[1];
	  picotrk->mPz = primtrk->p[2];
	  picotrk->mE  = primtrk->e;
	  picotrk->mEta = primtrk->eta;
	  picotrk->mXProj = projshowerxyz2.x();
	  picotrk->mYProj = projshowerxyz2.y();
	  picotrk->mZProj = projshowerxyz2.z();
	  
	  //double phi = atan2(primtrk->p[1],primtrk->p[0]);
	  //double theta = 2.0*atan(exp(-1.0*primtrk->eta));
	  //double mass = sqrt(primtrk->e*primtrk->e - primtrk->ptot*primtrk->ptot);
	  std::cout << "|primtrk|Id:"<<primtrk->id << "|Pid:"<<primtrk->ge_pid << "|E:"<<primtrk->e
		    << "|px:"<<primtrk->p[0] << "|py:"<<primtrk->p[1] << "|pz:"<<primtrk->p[2] << "|pt:"<<primtrk->pt << "|ptot:"<<primtrk->ptot
		    << "|eta:"<<primtrk->eta //<< "|theta:"<<theta << "|phi:"<<phi << "|mass:"<< mass
	    //<< "|point:("<<pointxyz.x() <<","<<pointxyz.y()<<","<<pointxyz.z() <<")"
	    //<< "|projE:("<<projshowerxyz1.x() <<","<<projshowerxyz1.y()<<","<<projshowerxyz1.z() <<")"
	    //<< "|projS:("<<projshowerxyz2.x() <<","<<projshowerxyz2.y()<<","<<projshowerxyz2.z() <<")"
		    << "|start_vertex:"<<primtrk->start_vertex_p << "|stop_vertex:"<<primtrk->stop_vertex_p
		    << "|frac:"<<frac << "|ntrk:"<<ntrk
		    << std::endl;
	  /*
	  std::cout << "|primtrk|Id:"<<picotrk->id << "|Pid:"<<picotrk->ge_pid << "|E:"<<picotrk->mE
		    << "|px:"<<picotrk->mPx << "|py:"<<picotrk->mPy << "|pz:"<<picotrk->mPz
		    << "|pt:"<<picotrk->pt()<<"|g2t:"<<primtrk->pt << "|ptot:"<<picotrk->ptot() << "|g2t:"<<primtrk->ptot
		    << "|eta:"<<picotrk->mEta <<"|g2t:"<<primtrk->eta << "|theta:"<<picotrk->theta() << "|phi:"<<picotrk->phi()
		    << "|point:("<<pointxyz.x() <<","<<pointxyz.y()<<","<<pointxyz.z() <<")"
		    << "|projE:("<<projshowerxyz1.x() <<","<<projshowerxyz1.y()<<","<<projshowerxyz1.z() <<")"
		    << "|projS:("<<projshowerxyz2.x() <<","<<projshowerxyz2.y()<<","<<projshowerxyz2.z() <<")"
		    << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	  */
	  int totalntrk = ntrk;
	  for( int itrk=0; itrk<totalntrk; ++itrk ){
	    const g2t_track_st* sectrk = mFcsDb->getPrimaryG2tTrack(pointclus,g2ttrk,frac,ntrk,itrk);
	    if( sectrk==0 ){ continue; }
	    double phi = atan2(sectrk->p[1],sectrk->p[0]);
	    double theta = 2.0*atan(exp(-1.0*sectrk->eta));
	    double mass = sqrt(sectrk->e*sectrk->e - sectrk->ptot*sectrk->ptot);
	    double diff = speedperc  - sectrk->ptot/sectrk->e;
	    std::cout << " + |sectrk:"<<itrk<<"|Id:"<<sectrk->id << "|Pid:"<<sectrk->ge_pid << "|E:"<<sectrk->e
		      << "|px:"<<sectrk->p[0] << "|py:"<<sectrk->p[1] << "|pz:"<<sectrk->p[2] << "|pt:"<<sectrk->pt << "|ptot:"<<sectrk->ptot
		      << "|eta:" << sectrk->eta << "|theta:"<<theta << "|phi:"<<phi << "|mass:"<< mass << "|beta:"<<sectrk->ptot/sectrk->e
		      << "|diff:"<< diff
		      << "|start_vertex:"<<sectrk->start_vertex_p << "|stop_vertex:"<<sectrk->stop_vertex_p
		      << "|frac:"<<frac << "|ntrk:"<<ntrk
		      << std::endl;
	  }
	  for( int i=0; i<vertexTable->GetNRows(); ++i ){
	    double xdiff = g2tvert[i].ge_x[0] - g2tvert[0].ge_x[0];
	    double ydiff = g2tvert[i].ge_x[1] - g2tvert[0].ge_x[1];
	    double zdiff = g2tvert[i].ge_x[2] - g2tvert[0].ge_x[2];
	    //double dist = sqrt( g2tvert[i].ge_x[0]*g2tvert[i].ge_x[0] + g2tvert[i].ge_x[1]*g2tvert[i].ge_x[1] + g2tvert[i].ge_x[2]*g2tvert[i].ge_x[2]  );
	    double dist = sqrt( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff );
	    double speed = (dist/g2tvert[i].ge_tof) / 100.0; //convert to meters/sec
	    speedperc = speed/299792458.0;
	    double masspi0 = 0.1349768;
	    //double masspi0 = 0.135031;
	    double MyE = masspi0 * sqrt( (1.0)/(1.0-speedperc*speedperc) );
	    double MyP = sqrt( MyE*MyE - masspi0*masspi0 );
	    double MyM = sqrt( MyE*MyE - MyP*MyP );
	    std::cout << "|g2tvert|i:"<<i << "|dp:"<<(g2tvert[i].daughter_p)+1 << "|eg_label:"<<g2tvert[i].eg_label << "|eg_proc:" << g2tvert[i].eg_proc
		      << "|event_p:"<<g2tvert[i].event_p << "|id:"<< g2tvert[i].id << "|is_itrmd:"<<g2tvert[i].is_itrmd << "|parent_p:"<<g2tvert[i].parent_p
		      << "|eg_tof:"<<g2tvert[i].eg_tof << "|eg_x:"<<g2tvert[i].eg_x[0] <<","<< g2tvert[i].eg_x[1]<<","<<g2tvert[i].eg_x[2]
		      << "|ge_tof:"<<g2tvert[i].ge_tof << "|ge_x:"<<g2tvert[i].ge_x[0]<<","<<g2tvert[i].ge_x[1]<<","<<g2tvert[i].ge_x[2]
		      << "|ge_vz:"<<speedperc << "|ge_MyE:"<< MyE << "|ge_MyP:"<< MyP << "|ge_MyM:"<< MyM
		      << std::endl;
	  }
	}
      }
    }
  }

  mDataTree->Fill();
  mHitsArr->Clear();
  mClusArr->Clear();
  mPointArr->Clear();
  mG2tPrimArr->Clear();
  
  return kStOk;
}

