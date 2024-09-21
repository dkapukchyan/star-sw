#include "StEvent/StEnumerations.h"
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
}

StMuFcsPi0TreeMaker::~StMuFcsPi0TreeMaker()
{
  delete mPi0Tree;
  delete mH1F_Entries;
  delete mEvtInfo;
  delete mPhArr;
  delete mPi0Arr;
  if( mFile_Output!=0 ){ mFile_Output->Close(); }
  delete mFile_Output;
}

Int_t StMuFcsPi0TreeMaker::Init()
{
  if( mFilename.Length() == 0){ mFilename="test.root"; } //Ensure a TFile is always created
  mFile_Output = new TFile(mFilename.Data(), "RECREATE");
  mPi0Tree     = new TTree("Pi0Tree","Tree with FcsPi0Candidate");
  mEvtInfo     = new TClonesArray("FcsEventInfo");
  mPhArr       = new TClonesArray("FcsPhotonCandidate");
  mPi0Arr      = new TClonesArray("FcsPi0Candidate");

  mPi0Tree->Branch("EventInfo",&mEvtInfo);
  mPi0Tree->Branch("TriggerInfo",0,"NTrig/I:Trig[NTrig]/I");
  ((TLeaf*)mPi0Tree->GetBranch("TriggerInfo")->GetListOfLeaves()->At(0))->SetAddress(&mNTrig);
  ((TLeaf*)mPi0Tree->GetBranch("TriggerInfo")->GetListOfLeaves()->At(1))->SetAddress(&mTriggers);
  mPi0Tree->Branch("Photon",&mPhArr);
  mPi0Tree->Branch("Pi0",&mPi0Arr);

  //std::cout << "|mEvtInfo:"<<mEvtInfo << "|bEvtInfo:"<< mPi0Tree->GetBranch("EventInfo")->GetAddress() << std::endl;
  
  mH1F_Entries = new TH1F("H1_Entries", "Number of entries;", 3,-0.5, 2.5);
  mH1F_Entries->Sumw2();

  mFcsTrigMap = (StFcsRun22TriggerMap*)GetMaker("fcsRun22TrigMap");
  if( mFcsTrigMap==0 ){ LOG_WARN << "StMuFcsRun22QaMaker::Init() - No Trigger Map found" << endm; }
  else{ if( mFcsTrigMap->sizeOfTriggers()<=0 ){ LOG_WARN << "StMuFcsRun22QaMaker::Init() - Trigger Map is empty" << endm; } }

  return kStOk;
}

void StMuFcsPi0TreeMaker::LoadDataFromFile(TFile* file, TTree* tree, TClonesArray* arr, TH1* hist)
{
  if( file==0 || file->IsZombie() ){ std::cout << "LoadDataFromFile - ERROR:Unable to load from null file or zombie file" << std::endl; return; }
  if( tree!=0 ){ std::cout << "LoadDataFromFile - WARNING:Overwriting TTree pointer" << std::endl; }
  tree = (TTree*)file->Get("Pi0Tree");
  if( tree==0 ){ std::cout << "LoadDataFromFile - ERROR:Pi0Tree not found in file" << std::endl; return; }
  if( arr!=0 ){ std::cout << "LoadDataFromFile - WARNING:Overwriting TClonesArray pointer" << std::endl; }
  arr = new TClonesArray("FcsPi0Candidate");
  tree->SetBranchAddress("Pi0",&arr);
  if( hist!=0 ){ std::cout << "LoadDataFromFile - WARNING:Overwriting TH1 pointer" << std::endl; }
  hist = (TH1*)file->Get("H1_Entries");
  if( hist==0 ){ std::cout << "LoadDataFromFile - WARNING:Entries histogram not found in file" << std::endl; }
  return;
}

Int_t StMuFcsPi0TreeMaker::InitRun(int runnumber) {
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  //mFcsDb->setDbAccess(0);
  if(!mFcsDb){
    LOG_ERROR << "StMuFcsPi0TreeMaker::InitRun - Failed to get StFcsDbMaker" << endm;
    return kStFatal;
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
    mPi0Tree->Write();
    mH1F_Entries->Write();
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
  if( mIgnoreTrig ){ validtrigfound = true; } //If ignoring triggers then set validtrigfound to true so event is not skipped
  StMuTriggerIdCollection* TrigMuColl = &(mMuEvent->triggerIdCollection());
  if( !TrigMuColl ){ LOG_ERROR <<"StMuFcsPi0TreeMaker::FillEventInfo - !TrigMuColl" <<endl; return kStErr; }
  const StTriggerId& trgIDs = TrigMuColl->nominal();
  mNTrig = trgIDs.triggerIds().size();
  for( Int_t i=0; i<mNTrig; ++i ){
    unsigned int trig = trgIDs.triggerId(i);
    mTriggers[i] = trig;
    if( !mIgnoreTrig ){
      if( mFcsTrigMap!=0 ){
	std::string thistrig = mFcsTrigMap->nameFromId(trig,mMuEvent->runNumber());
	for( unsigned int j=0; j<mTargetTrig.size(); ++j ){
	  if( mTargetTrig.at(j)==thistrig ){ validtrigfound=true; }
	}
      }
    }
  }
  if( !validtrigfound ){ return kStSkip; }

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

  FcsEventInfo* evtinfo = (FcsEventInfo*) mEvtInfo->ConstructedAt(0);  
  evtinfo->mRunTime         = mMuEvent->eventInfo().time();
  evtinfo->mRunNum          = mMuEvent->runNumber();
  evtinfo->mFill            = mRunInfo->beamFillNumber(StBeamDirection::east);//Yellow beam
  evtinfo->mEvent           = mMuEvent->eventId();
  evtinfo->mBx48Id          = mTrigData->bunchId48Bit();
  evtinfo->mBx7Id           = mTrigData->bunchId7Bit();
  evtinfo->mTofMultiplicity = mTrigData->tofMultiplicity();

  //std::cout << "|runtime:"<<evtinfo->mRunTime << "|runnum:"<<evtinfo->mRunNum << "|event:"<<evtinfo->mEvent << std::endl;
  
  //Spin information
  if( mSpinDbMkr==0 ){
    Double_t rndm = mSpinRndm.Rndm(); //random number between 0 and 1
    //The numbers below are chosen for their bit representation as described and correspond to the source polarization. Need to flip to convert it to STAR polarization direction because of the Siberian Snakes; i.e. '+' -> '-' and '-' -> '+'
    if( rndm<=0.25 ){ evtinfo->mSpin = 5; }                    //Bits 0101 is B+ and Y+
    else if( 0.25<rndm && rndm<=0.5 ){ evtinfo->mSpin = 6; }   //Bits 0110 is B+ and Y-
    else if( 0.5<rndm && rndm<=0.75 ){ evtinfo->mSpin = 9; }   //Bits 1001 is B- and Y+
    else{ evtinfo->mSpin = 10; }                               //Bits 1010 is B- and Y-
  }
  else{
    evtinfo->mSpin = mSpinDbMkr->spin4usingBX48( mTrigData->bunchId48Bit() ); //This is also source polarization
  }
  

  //Vertex Information
  evtinfo->mVpdVz = -999;
  if( mMuDst->btofHeader() ){ evtinfo->mVpdVz = mMuDst->btofHeader()->vpdVz(); }
  //@[April 7, 2021] > No Slewing correction for BBC yet, see StFmsJetMaker2015 in BrightSTAR??
  evtinfo->mBbcTacDiff = mTrigData->bbcTimeDifference() - 4096; //subtract 4096 since 0 means bad event and distribution is Gaussian around 4096
  evtinfo->mBbcVz = -999;
  if( fabs(evtinfo->mBbcTacDiff)>1.e-6 ){ evtinfo->mBbcVz = evtinfo->mBbcTacDiff * -0.2475; } //0.2475 = 0.0165*30/2x
  //StZdcTriggerDetector& zdc = mMuEvent->zdcTriggerDetector();
  //std::cout <<"|ZdcV:"<<zdc.vertexZ() << std::endl;
  evtinfo->mZdcVz =  mTrigData->zdcVertexZ();

  if( mEpdQaMkr!=0 ){
    evtinfo->mEpdTacEarlyE = mEpdQaMkr->epdTacEarlyE();
    evtinfo->mEpdTacEarlyW = mEpdQaMkr->epdTacEarlyW();
    evtinfo->mEpdAvgE = mEpdQaMkr->epdTacAvgE();
    evtinfo->mEpdAvgW = mEpdQaMkr->epdTacAvgW();
    evtinfo->mEpdVz = mEpdQaMkr->epdVertex();
  }
  else{
    evtinfo->mEpdVz = -999;
  }

  short foundvertex = 0;   //Bit vectorfor knowing which vertex was used 0 means no vertex, 1=Vpd,2=Epd,4=Bbc
  if( evtinfo->mVpdVz > -998){ foundvertex |= 0x001; }
  else if( evtinfo->mEpdVz > -998 ){ foundvertex |= 0x010; }
  else if( evtinfo->mBbcVz > -998 ){ foundvertex |= 0x100; }
  else{ foundvertex = 0; } //Redundant but don't want a dangling "else if" statement
  evtinfo->mFoundVertex = foundvertex;

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

	//std::cout << "Cluster|detid:"<<ph->mDetId << "|mX:"<<ph->mX << "|mY:"<<ph->mY << "|mZ:"<<ph->mZ << std::endl;

	StLorentzVectorD iclu_p_withz;
	if( foundvertex == 0x001 ){ iclu_p_withz = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, evtinfo->mVpdVz ); }
	if( foundvertex == 0x010 ){ iclu_p_withz = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, evtinfo->mEpdVz ); }
	if( foundvertex == 0x100 ){ iclu_p_withz = mFcsDb->getLorentzVector( iclu_pos, iclu_energy, evtinfo->mBbcVz ); }
	ph->mPxVert = foundvertex!=0 ? iclu_p_withz.px() : 0;
	ph->mPyVert = foundvertex!=0 ? iclu_p_withz.py() : 0;
	ph->mPzVert = foundvertex!=0 ? iclu_p_withz.pz() : 0;
      }
    }
  }
  
  evtinfo->mClusterSize = ncandidates;
  Int_t clustersize = ncandidates; //local copy of evtinfo->mClusterSize
  
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

	ph->mEn = ipoi_energy;
	ph->mPxRaw = ipoi_p.px();
	ph->mPyRaw = ipoi_p.py();
	ph->mPzRaw = ipoi_p.pz();

	StLorentzVectorD ipoi_p_withz;
	if( foundvertex == 0x001 ){ ipoi_p_withz = mFcsDb->getLorentzVector( ipoi_pos, ipoi_energy, evtinfo->mVpdVz ); }
	if( foundvertex == 0x010 ){ ipoi_p_withz = mFcsDb->getLorentzVector( ipoi_pos, ipoi_energy, evtinfo->mEpdVz ); }
	if( foundvertex == 0x100 ){ ipoi_p_withz = mFcsDb->getLorentzVector( ipoi_pos, ipoi_energy, evtinfo->mBbcVz ); }
	ph->mPxVert = foundvertex!=0 ? ipoi_p_withz.px() : 0;
	ph->mPyVert = foundvertex!=0 ? ipoi_p_withz.py() : 0;
	ph->mPzVert = foundvertex!=0 ? ipoi_p_withz.pz() : 0;

	//std::cout << "Point|detid:"<<ph->mDetId << "|mX:"<<ph->mX << "|mY:"<<ph->mY << "|mZ:"<<ph->mZ << std::endl;
	
      }//i point
    }//fcs dets
  }

  //std::cout << "|ncandidates:"<<ncandidates <<"|clustersize:"<<clustersize <<"|Size:"<<mPhArr->GetEntriesFast() << std::endl;
  Int_t npi0candidate = 0;
  for( Int_t ic = 0; ic<clustersize; ++ic ){
    FcsPhotonCandidate* iclus = (FcsPhotonCandidate*) mPhArr->ConstructedAt(ic);
    if( ! iclus->mFromCluster ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
    if( ic==(clustersize-1) ){ continue; }
    for( Int_t jc=ic+1; jc<evtinfo->mClusterSize; jc++ ){
      FcsPhotonCandidate* jclus = (FcsPhotonCandidate*) mPhArr->ConstructedAt(jc);
      if( ! jclus->mFromCluster ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
      FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
      pi0c->mFromCluster = true;
      pi0c->mPhoton1Idx = ic;
      pi0c->mPhoton2Idx = jc;

      TLorentzVector pi0Vert_LV = iclus->lvVert() + jclus->lvVert();
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
    }
  }
  for( Int_t ip = clustersize; ip<mPhArr->GetEntriesFast(); ++ip ){
    FcsPhotonCandidate* ipoi = (FcsPhotonCandidate*) mPhArr->ConstructedAt(ip);
    if( ipoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
    if( ip==(mPhArr->GetEntriesFast()-1) ){ continue; }
    for( Int_t jp=ip+1; jp<mPhArr->GetEntriesFast(); jp++ ){
      FcsPhotonCandidate* jpoi = (FcsPhotonCandidate*) mPhArr->ConstructedAt(jp);
      if( jpoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
      FcsPi0Candidate* pi0c = (FcsPi0Candidate*) mPi0Arr->ConstructedAt(npi0candidate++);
      pi0c->mFromCluster = false;
      pi0c->mPhoton1Idx = ip;
      pi0c->mPhoton2Idx = jp;

      TLorentzVector pi0Vert_LV = ipoi->lvVert() + jpoi->lvVert();
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
    }
  }

  mPi0Tree->Fill();

  mEvtInfo->Clear();
  mNTrig = 0; //Since ROOT only writes up to the size of mNTrig then only need to reset this back to zero and next loop will overwrite array as neccessary
  mPhArr->Clear();
  mPi0Arr->Clear();
    
  return kStOk;
}

