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
  delete mHistsArr;
  delete mDataTree;
  
}

Int_t StFcsShowerAnaMaker::Init()
{
  if( !mFcsDb ){ mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb")); }
  //mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
    LOG_ERROR << "StFcsShowerAnaMaker::InitRun Failed to get StFcsDb" << endm;
    return kStFatal;
  }

  mHistsArr = new TObjArray();
  if( LoadHistograms(mHistsArr) ){ std::cout << "StFcsShowerAnaMaker::InitRun Failed to make new histgorams" << std::endl; }
  
  if( mDataTree==0 ){
    mDataTree = new Rtools::ClonesArrTree("ShowerAnaTree","");
    mDataTree->AddArr("FcsHit",    "StFcsPicoHit");
    mDataTree->AddArr("FcsCluster","StFcsPicoCluster");
    mDataTree->AddArr("FcsPoint",  "StFcsPicoPoint");
    mDataTree->AddArr("G2tPrim",   "StPicoG2tTrack");
    mDataTree->AddArr("G2tParent", "StPicoG2tTrack");
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
  outfile->Close();
  delete outfile;
  return kStOk;
}

bool StFcsShowerAnaMaker::LoadHistograms(TObjArray* arr, TFile* file)
{
  if( arr->GetEntriesFast()!=0 ){ return true; }
  UInt_t nloaded = 0;

  nloaded += Rtools::LoadH1(arr,file,mH1F_PointXLocal,"H1F_PointXLocal","",100,0,1);
  nloaded += Rtools::LoadH1(arr,file,mH1F_PointYLocal,"H1F_PointYLocal","",100,0,1);
  nloaded += Rtools::LoadH1(arr,file,mH1F_ClusSigMax,"H1F_ClusSigMax","",200,0,2);
  nloaded += Rtools::LoadH1(arr,file,mH1F_ClusSigMin,"H1F_ClusSigMin","",200,0,2);
  nloaded += Rtools::LoadH2(arr,file,mH2F_PointXProjX,"H2F_PointXProjX","",200,-100,100,200,-100,100);
  nloaded += Rtools::LoadH2(arr,file,mH2F_PointYProjY,"H2F_PointYProjY","",200,-100,100,200,-100,100);  
  nloaded += Rtools::LoadH2(arr,file,mH2F_ClusSigMaxEn,"H2F_ClusSigMaxEn","",100,0,100,200,0,2);
  nloaded += Rtools::LoadH2(arr,file,mH2F_ClusSigMinEn,"H2F_ClusSigMinEn","",100,0,100,200,0,2);
  
  nloaded += Rtools::LoadH1(arr,file,mH1F_primid,"H1F_primid",";GEANT ID", 11,-0.5,10.5);
  nloaded += Rtools::LoadH1(arr,file,mH1F_parentid,"H1F_parentid",";GEANT ID", 11,-0.5,10.5);
  nloaded += Rtools::LoadH2(arr,file,mH2F_npoiVnclus,"H2F_npoiVnclus",";NClusters;NPoints", 7,-0.5,6.5, 7,-0.5,6.5);
  nloaded += Rtools::LoadH2(arr,file,mH2F_cluseVlore,"H2F_cluseVlore",";Cluster Lorentz E;Cluster E", 100,0,50, 100,0,50);
  nloaded += Rtools::LoadH2(arr,file,mH2F_trkeVpoie,"H2F_trkeVpoie",";point E;trk E", 300,0,30, 300,0,30);
  nloaded += Rtools::LoadH1(arr,file,mH1F_invmasspoi,"H1F_invmasspoi",";inv mass point (GeV)", 100,0,1);
  nloaded += Rtools::LoadH1(arr,file,mH1F_invmasstrk,"H1F_invmasstrk",";inv mass g2trk (GeV)", 100,0,1);
  nloaded += Rtools::LoadH1(arr,file,mH1F_dpoitrk,"H1F_dpoitrk",";r (cm)", 100,0,10);
  nloaded += Rtools::LoadH2(arr,file,mH2F_massVdgg,"H2F_massVdgg",";d_{gg} point (cm);inv masss point (GeV)", 100,0,50, 100,0,1);
  nloaded += Rtools::LoadH2(arr,file,mH2F_trkmassVdgg,"H2F_trkmassVdgg",";d_{gg} point (cm);inv masss track (GeV)", 100,0,50, 100,0,1);

  std::cout << "|nloaded:"<<nloaded << std::endl;
  //If number of new histograms == array size then all histograms are new so return true i.e. it is true that all histograms were "made"
  if( nloaded == arr->GetEntriesFast() ){ return true; }
  else{ return false; }  //i.e. It is not true that all histograms were made
}

void StFcsShowerAnaMaker::CleanHists()
{
  for( UInt_t i=0; i<mHistsArr->GetEntriesFast(); ++i ){ delete mHistsArr->At(i); }
  mHistsArr->Clear();
}

void StFcsShowerAnaMaker::WriteHists()
{
  for( UInt_t i=0; i<mHistsArr->GetEntriesFast(); ++i ){ mHistsArr->At(i)->Write(); }
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
  
  TClonesArray* mHitsArr = mDataTree->GetBranchArr("FcsHit");
  TClonesArray* mClusArr = mDataTree->GetBranchArr("FcsCluster");
  TClonesArray* mPointArr = mDataTree->GetBranchArr("FcsPoint");
  TClonesArray* mG2tPrimArr = mDataTree->GetBranchArr("G2tPrim");
  TClonesArray* mG2tParArr = mDataTree->GetBranchArr("G2tParent");
  //std::cout << "|HitsArr:"<<mHitsArr << "|ClusArr:"<<mClusArr << "|PointArr:"<<mPointArr << "|G2tPrimArr:"<<mG2tPrimArr << "|G2tParArr:"<<mG2tParArr << std::endl;
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
      picoclus->mPx = clusp.px();
      picoclus->mPy = clusp.py();
      picoclus->mPz = clusp.pz();
      picoclus->mE  = clusp.e();

      mH2F_cluseVlore->Fill(clusp.e(),cluster->energy());
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
	  if( GetDebug()>0 ){ std::cout << "g2t_track table has "<< nTrk << " tracks" << std::endl; }
	  if( nTrk>0 ){
	    g2ttrk = trackTable->GetTable();
	    if( !g2ttrk ) { std::cout << " g2t_track GetTable failed" << std::endl; continue; }
	  }
	}
	if( !vertexTable ){ std::cout<< "g2t_vertex Table not found" << std::endl; continue; }
	else{
	  const int nVertex = vertexTable->GetNRows();
	  if( GetDebug()>0 ){ std::cout << "g2t_vertex table has "<< nVertex << " vertices" << std::endl; }
	  if( nVertex>0 ){
	    g2tvert = vertexTable->GetTable();
	    if( !g2tvert) { std::cout << " g2t_vertex GetTable failed" << std::endl; continue; }
	  }
	}
	if( g2ttrk && g2tvert ){
	  StFcsPicoPoint* picopoint = (StFcsPicoPoint*)mPointArr->ConstructedAt(totalpoints);
	  picopoint->mDetId                 = point->detectorId();
	  picopoint->mEnergy                = point->energy();
	  picopoint->mXlocal                = point->x();
	  picopoint->mYlocal                = point->y();
	  picopoint->mParentClusterId       = point->parentClusterId();
	  picopoint->mNParentClusterPhotons = point->nParentClusterPhotons();
	  //STAR xyx
	  StThreeVectorF pointxyz = point->xyz();
	  //std::cout << "|ipoint:"<<ipoint << "|E:"<<point->energy() << "|x:"<<pointxyz.x() << "|y:"<<pointxyz.y() << "|z:"<<pointxyz.z() << std::endl;
	  picopoint->mXstar = pointxyz.x();
	  picopoint->mYstar = pointxyz.y();
	  picopoint->mZstar = pointxyz.z();
	  //Lorentz 4 momentum of point
	  StLorentzVectorD pointp = point->fourMomentum();
	  picopoint->mPx = pointp.px();
	  picopoint->mPy = pointp.py();
	  picopoint->mPz = pointp.pz();
	  picopoint->mE  = pointp.e();

	  float frac=0;
	  int ntrk=0;
	  //const g2t_track_st* parenttrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  //std::cout << "|parenttrk|Id:"<<parenttrk->id << "|Pid:"<<parenttrk->ge_pid << "|E:"<<parenttrk->e << "|eta:"<<parenttrk->eta << "|frac:"<<frac << "|ntrk:"<<ntrk << std::endl;
	  const g2t_track_st* primtrk = mFcsDb->getPrimaryG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  //const g2t_track_st* primtrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  StPicoG2tTrack* picotrk = (StPicoG2tTrack*)mG2tPrimArr->ConstructedAt(totalpoints);
	  StThreeVectorD projshowerxyz1 = mFcsDb->projectTrackToEcal(primtrk,g2tvert);
	  StThreeVectorD projshowerxyz2 = mFcsDb->projectTrackToEcalSMax(primtrk,g2tvert);
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

	  mH1F_primid->Fill(primtrk->ge_pid);
	  mH2F_PointXProjX->Fill(pointxyz.x(), projshowerxyz2.x());
	  mH2F_PointYProjY->Fill(pointxyz.y(), projshowerxyz2.y());
	  mH2F_trkeVpoie->Fill( point->energy(), primtrk->e );

	  mH2F_trkeVpoie->Fill( picopoint->mEnergy, picotrk->mE );
	  double xdiff = picopoint->mXstar - picotrk->mXProj;
	  double ydiff = picopoint->mYstar - picotrk->mYProj;
	  mH1F_dpoitrk->Fill( sqrt(xdiff*xdiff + ydiff*ydiff) );

	  const g2t_track_st* parenttrk = mFcsDb->getParentG2tTrack(pointclus,g2ttrk,frac,ntrk);
	  StPicoG2tTrack* picoparent = (StPicoG2tTrack*)mG2tParArr->ConstructedAt(totalpoints++);
	  StThreeVectorD projparentxyz = mFcsDb->projectTrackToEcalSMax(parenttrk,g2tvert);

	  picoparent->id = parenttrk->id;
	  picoparent->ge_pid = parenttrk->ge_pid;
	  picoparent->mPx = parenttrk->p[0];
	  picoparent->mPy = parenttrk->p[1];
	  picoparent->mPz = parenttrk->p[2];
	  picoparent->mE  = parenttrk->e;
	  picoparent->mEta = parenttrk->eta;
	  picoparent->mXProj = projparentxyz.x();
	  picoparent->mYProj = projparentxyz.y();
	  picoparent->mZProj = projparentxyz.z();

	  mH1F_parentid->Fill(parenttrk->ge_pid);

	  //double phi = atan2(primtrk->p[1],primtrk->p[0]);
	  //double theta = 2.0*atan(exp(-1.0*primtrk->eta));
	  //double mass = sqrt(primtrk->e*primtrk->e - primtrk->ptot*primtrk->ptot);
	  if( GetDebug()>1 ){
	    std::cout << "|primtrk|Id:"<<primtrk->id << "|Pid:"<<primtrk->ge_pid << "|E:"<<primtrk->e
		      << "|px:"<<primtrk->p[0] << "|py:"<<primtrk->p[1] << "|pz:"<<primtrk->p[2] << "|pt:"<<primtrk->pt << "|ptot:"<<primtrk->ptot
		      << "|eta:"<<primtrk->eta //<< "|theta:"<<theta << "|phi:"<<phi << "|mass:"<< mass
		      << "|point:("<<pointxyz.x() <<","<<pointxyz.y()<<","<<pointxyz.z() <<")"
		      << "|projE:("<<projshowerxyz1.x() <<","<<projshowerxyz1.y()<<","<<projshowerxyz1.z() <<")"
		      << "|projS:("<<projshowerxyz2.x() <<","<<projshowerxyz2.y()<<","<<projshowerxyz2.z() <<")"
		      << "|start_vertex:"<<primtrk->start_vertex_p << "|stop_vertex:"<<primtrk->stop_vertex_p
		      << "|frac:"<<frac << "|ntrk:"<<ntrk
		      << std::endl;
	  }
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
	  if( GetDebug()>1 ){
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
    if( mPointArr->GetEntriesFast()>=2 ){
      TLorentzVector p1; //Lorentz vector of 1 point
      TLorentzVector p2; //Lorentz vector of other point
      StFcsPicoPoint* poi0 = (StFcsPicoPoint*)mPointArr->At(0);
      p1.SetPxPyPzE( poi0->mPx, poi0->mPy, poi0->mPz, poi0->mE );
      StFcsPicoPoint* poi1 = (StFcsPicoPoint*)mPointArr->At(1);
      p2.SetPxPyPzE( poi1->mPx, poi1->mPy, poi1->mPz, poi1->mE );
      TLorentzVector inv = p1+p2;
      double xdiff = poi0->mXstar - poi1->mXstar;
      double ydiff = poi0->mYstar - poi1->mYstar;
      double zdiff = poi0->mZstar - poi1->mZstar;
      double dgg = sqrt( xdiff*xdiff + ydiff*ydiff + zdiff*zdiff );
      mH1F_invmasspoi->Fill(inv.M());
      mH2F_massVdgg->Fill(dgg,inv.M());
    }
    if( mG2tPrimArr->GetEntriesFast()>=2 ){
      TLorentzVector p1; //Lorentz vector of 1 g2t track
      TLorentzVector p2; //Lorentz vector of other g2t track
      StPicoG2tTrack* trk0 = (StPicoG2tTrack*)mG2tPrimArr->At(0);
      p1.SetPxPyPzE( trk0->mPx, trk0->mPy, trk0->mPz, trk0->mE );
      StPicoG2tTrack* trk1 = (StPicoG2tTrack*)mG2tPrimArr->At(1);
      p2.SetPxPyPzE( trk1->mPx, trk1->mPy, trk1->mPz, trk1->mE );
      TLorentzVector inv = p1+p2;
      double xdiff = trk0->mXProj - trk1->mXProj;
      double ydiff = trk0->mYProj - trk1->mYProj;
      //double zdiff = trk0->mZProj - trk1->mZProj;
      double dgg = sqrt( xdiff*xdiff + ydiff*ydiff );
      mH1F_invmasstrk->Fill(inv.M());
      mH2F_trkmassVdgg->Fill(dgg,inv.M());
    }
  }

  mH2F_npoiVnclus->Fill(totalclus,totalpoints);

  mDataTree->Fill();
  mHitsArr->Clear();
  mClusArr->Clear();
  mPointArr->Clear();
  mG2tPrimArr->Clear();
  mG2tParArr->Clear();
  
  return kStOk;
}

