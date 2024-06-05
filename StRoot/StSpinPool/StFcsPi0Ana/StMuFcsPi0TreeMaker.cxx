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

ClassImp(FcsPi0Info)

FcsPi0Info::FcsPi0Info()
{}

FcsPi0Info::~FcsPi0Info()
{}

ClassImp(StMuFcsPi0TreeMaker)

StMuFcsPi0TreeMaker::StMuFcsPi0TreeMaker(const Char_t* name) : StMaker(name)
{}

StMuFcsPi0TreeMaker::~StMuFcsPi0TreeMaker()
{
  delete mH1F_Entries;
  delete mPi0Tree;
  delete mPi0Arr;
  delete mFile_Output;
}

Int_t StMuFcsPi0TreeMaker::Init()
{
  if( mFilename.Length() == 0){ mFilename="test.root"; } //Ensure a TFile is always created
  mFile_Output = new TFile(mFilename.Data(), "RECRATE");
  mPi0Tree = new TTree("Pi0Tree","Tree with FcsPi0Info");
  mPi0Arr = new TClonesArray("FcsPi0Info");
  mPi0Tree->Branch("Pi0",&mPi0Arr);
  
  mH1F_Entries = new TH1F("H1_Entries", "Number of entries;", 3,-0.5, 2.5);
  mH1F_Entries->Sumw2();

  return kStOk;
}

void StMuFcsPi0TreeMaker::LoadDataFromFile(TFile* file, TTree* tree, TClonesArray* arr, TH1* hist)
{
  if( file==0 || file->IsZombie() ){ std::cout << "LoadDataFromFile - ERROR:Unable to load from null file or zombie file" << std::endl; return; }
  if( tree!=0 ){ std::cout << "LoadDataFromFile - WARNING:Overwriting TTree pointer" << std::endl; }
  tree = (TTree*)file->Get("Pi0Tree");
  if( tree==0 ){ std::cout << "LoadDataFromFile - ERROR:Pi0Tree not found in file" << std::endl; return; }
  if( arr!=0 ){ std::cout << "LoadDataFromFile - WARNING:Overwriting TClonesArray pointer" << std::endl; }
  arr = new TClonesArray("FcsPi0Info");
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
  return kStOK;
}

//-----------------------
Int_t StMuFcsPi0TreeMaker::Finish() {

  mFile_Output->cd();
  mPi0Tree->Write();
  mH1F_Entries->Write();
  mFile_Output->Close();
  
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

  mMuFcsColl = mMuDst->muFcsCollection();
  if (!mMuFcsColl) { LOG_ERROR << "StMuFcsPi0TreeMaker::Make did not find MuFcsCollection" << endm; return kStErr; }

  mH1F_Entries->Fill(1);
  /*
  //TOF mult cut
  int tofMult = 0;
  //const StTriggerData* trgdata = event->triggerData();
  //if(!trgdata && StMuDst::event()) trgdata = StMuDst::event()->triggerData();
  if(mTrigData){
    tofMult = mTrigData->tofMultiplicity();
    LOG_DEBUG<<"TOF mult="<<tofMult<<endm;
    if (tofMult > 100) return kStOK;
  }else{
    LOG_WARN << "No TriggerData found in Mudst. No TOFMult cut"<<endm;
  }

  //ZVERTEX
  Double_t zTPC=-999.0;
  StMuPrimaryVertex* tpcvtx = mMuDst->primaryVertex();
  if(tpcvtx){ zTPC=tpcvtx->position().z(); }
  Double_t zVPD = -999.0;     
  if(mMuDst->btofHeader()){ zVPD=mMuDst->btofHeader()->vpdVz(); }
  Double_t zBBC=-999.0;
  if(mTrigData){ zBBC = (4096 - mTrigData->bbcTimeDifference())*0.016*30.0/2.0; }  
  
  for (int det = 0; det < 2; det++) {
    TClonesArray* clusters = mMuFcsColl->getClusterArray();
    if( clusters!=0 ){
      int nc = mMuFcsColl->numberOfClusters(det);
    }
        
    for( int iclus=mMuFcsColl->indexOfFirstCluster(det); iclus<nc; iclus++) {
      if( i == nc-1 ){ break; }
      StMuFcsCluster* clu = (StMuFcsCluster*)clusters->At(iclus);
      float clu_energy = clu->energy();
      float clu_x = clu->x();
      float clu_y = clu->y();
      StThreeVectorD cluPos = mFcsDb->getStarXYZfromColumnRow(det, clu_x, clu_y);
      float cluPos_x = cluPos.x();
      float cluPos_y = cluPos.y();
      StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(det, clu_x, clu_y);
      StLorentzVectorD p = mFcsDb->getLorentzVector((xyz), clu_energy, 0);
      
      for( int jclus=iclus+1; jclus<nc; jclus++) {
	StMuFcsCluster* cluj = (StMuFcsCluster*)clusters->At(jclus);
	float cluj_energy = cluj->energy();
	float cluj_x = cluj->x();
	float cluj_y = cluj->y();
	StThreeVectorD clujPos = mFcsDb->getStarXYZfromColumnRow(det, cluj_x, cluj_y);
	
	float zgg = (fabs(clu_energy - cluj_energy)) / (clu_energy + cluj_energy);
	StThreeVectorD xyzj = mFcsDb->getStarXYZfromColumnRow(det, cluj->x(), cluj->y());
	StLorentzVectorD pj = mFcsDb->getLorentzVector((xyzj), cluj->energy(), 0);
	
	if ((clu_energy + cluj_energy) > bestclu_totalE) {
	  check_fillclu = 1;
	  bestclu_invmass = ((p + pj).m());
	  
	  if (zTPC > -999) 
	    {
	      StThreeVectorD xyzj_Vtpc = xyzj;
	      StLorentzVectorD pj_Vtpc = mFcsDb->getLorentzVector((xyzj_Vtpc), cluj_energy, zTPC);
	      StThreeVectorD xyz_Vtpc = xyz;
	      StLorentzVectorD p_Vtpc = mFcsDb->getLorentzVector((xyz_Vtpc), clu_energy, zTPC);
	      bestclu_invmass_Vtpc = ((p_Vtpc + pj_Vtpc).m());
	      bestclu_invmass_Vz0tpc = bestclu_invmass;
	    }
	  if (zBBC > -999) 
	    {
	      StThreeVectorD xyzj_Vbbc = xyzj;
	      StLorentzVectorD pj_Vbbc = mFcsDb->getLorentzVector((xyzj_Vbbc), cluj_energy, zBBC);
	      StThreeVectorD xyz_Vbbc = xyz;
	      StLorentzVectorD p_Vbbc = mFcsDb->getLorentzVector((xyz_Vbbc), clu_energy, zBBC);
	      bestclu_invmass_Vbbc = ((p_Vbbc + pj_Vbbc).m());
	      if (zTPC > -999) {bestclu_invmass_Vbbctpc = bestclu_invmass_Vbbc;}
	    }
	  if (zVPD > -999) 
	    {
	      StThreeVectorD xyzj_Vvpd = xyzj;
	      StLorentzVectorD pj_Vvpd = mFcsDb->getLorentzVector((xyzj_Vvpd), cluj_energy, zVPD);
	      StThreeVectorD xyz_Vvpd = xyz;
	      StLorentzVectorD p_Vvpd = mFcsDb->getLorentzVector((xyz_Vvpd), clu_energy, zVPD);
	      bestclu_invmass_Vvpd = ((p_Vvpd + pj_Vvpd).m());
	      if (zTPC > -999) {bestclu_invmass_Vvpdtpc = bestclu_invmass_Vvpd;}
	    }
	  
	  bestclu_totalE = (clu_energy + cluj_energy);
	  bestclu_dgg = (sqrt((xyz[0] - xyzj[0]) * (xyz[0] - xyzj[0]) + (xyz[1] - xyzj[1]) * (xyz[1] - xyzj[1]) + (xyz[2] - xyzj[2]) * (xyz[2] - xyzj[2])));
	  bestclu_Zgg = fabs((clu_energy - cluj_energy) / (clu_energy + cluj_energy));
	  bestclu_opening_angle = acos((xyz[0] * xyzj[0] + xyz[1] * xyzj[1] + xyz[2] * xyzj[2]) / (sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]) * sqrt(xyzj[0] * xyzj[0] + xyzj[1] * xyzj[1] + xyzj[2] * xyzj[2])));
	  
	    int pnti_nTowers = clu->nTowers();
	    int pntj_nTowers = cluj->nTowers();
	    TRefArray* clui_hits = clu->hits();
	    TRefArray* cluj_hits = cluj->hits();
	    float max_energy_tower = -1;
	    unsigned short tower_id = 0;
	    for( int k=0; k<pnti_nTowers; k++ ){
	      StMuFcsHit* clushit = (StMuFcsHit*)clui_hits->At(k);
	      if( clushit->energy() > max_energy_tower ){
		max_energy_tower = clushit->energy();
		tower_id = clushit->id();
		best_tower_det_cluster = det;
	      }
	    }
	    best_tower_id1_cluster = tower_id;
	    max_energy_tower = -1;
	    for( int k = 0; k<pntj_nTowers; k++){
	      StMuFcsHit* clushit = (StMuFcsHit*)cluj_hits->At(k);
	      if( clushit->energy() > max_energy_tower) {
		max_energy_tower = clushit->energy();
		tower_id = clushit->id();
	      }
	    }
	    best_tower_id2_cluster = tower_id;
	  }
	}
      }
    TClonesArray* points = mMuFcsColl->getPointArray();
    if( points!=0 ){
      int np = mMuFcsColl->numberOfPoints(det);

      
      for (int i = mMuFcsColl->indexOfFirstPoint(det); i < np; i++) {
	StMuFcsPoint* pnt = (StMuFcsPoint*)points->At(i);
	float pnt_x = pnt->x();
	float pnt_y = pnt->y();
	float pnt_energy = pnt->energy();
	if( pnt_energy > E_min ){ n_EcalPoint_cut++; }
	StThreeVectorD poiPos = mFcsDb->getStarXYZfromColumnRow(det, pnt_x, pnt_y);
	float poiPos_x = poiPos.x();
	float poiPos_y = poiPos.y();
	h2_point_position->Fill(poiPos_x, poiPos_y);
	h1_each_point_energy->Fill(pnt_energy);
	StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(det, pnt_x, pnt_y);
	StLorentzVectorD p = mFcsDb->getLorentzVector((xyz), pnt_energy, 0);
	if (i == np - 1) continue;
	for (int j = i + 1; j < np; j++) {
	  StMuFcsPoint* pntj = (StMuFcsPoint*)points->At(j);
	  float pntj_energy = pntj->energy();

	  h1_two_point_energy_nocut->Fill(pnt_energy + pntj_energy);
	  float zgg = (fabs(pnt_energy - pntj_energy)) / (pnt_energy + pntj_energy);
	  h1_Zgg_nocut_point->Fill(zgg);
	  StThreeVectorD xyzj = mFcsDb->getStarXYZfromColumnRow(det, pntj->x(), pntj->y());
	  StLorentzVectorD pj = mFcsDb->getLorentzVector((xyzj), pntj->energy(), 0);
	  h1_inv_mass_point_nocut->Fill((p + pj).m());

	  if( pntj->energy() < E_min ){ continue; }
	  if( zgg >= 0.7 ){ continue; }
	  if ((pnt_energy + pntj_energy) > bestpnt_totalE) {
	    check_fillpnt = 1;
	    bestpnt_invmass = ((p + pj).m());
	    bestpnt_totalE = (pnt_energy + pntj_energy);
	    bestpnt_dgg = (sqrt((xyz[0] - xyzj[0]) * (xyz[0] - xyzj[0]) + (xyz[1] - xyzj[1]) * (xyz[1] - xyzj[1]) + (xyz[2] - xyzj[2]) * (xyz[2] - xyzj[2])));
	    bestpnt_Zgg = fabs((pnt_energy - pntj_energy) / (pnt_energy + pntj_energy));
	    bestpnt_opening_angle = acos((xyz[0] * xyzj[0] + xyz[1] * xyzj[1] + xyz[2] * xyzj[2]) / (sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]) * sqrt(xyzj[0] * xyzj[0] + xyzj[1] * xyzj[1] + xyzj[2] * xyzj[2])));
	  }
	}
      }
    }

    h1_nCluster->Fill(total_nc);
    h1_nclu_good->Fill(n_EcalClust_cut);
    h1_nPoint->Fill(total_np);
    h1_npoi_good->Fill(n_EcalPoint_cut);
    h2_EcalMult_vs_TofMult->Fill(tofMult, n_Ecal_cut);
    if (n_Ecal_cut > 40) { return kStOK; }
    if (check_fillclu == 1) {
      h1_inv_mass_cluster->Fill(bestclu_invmass);
      h1_two_cluster_energy_allcut->Fill(bestclu_totalE);
      h2_cluster_invmass_vs_dgg->Fill(bestclu_dgg, bestclu_invmass);
      h2_cluster_invmass_vs_Zgg->Fill(bestclu_Zgg, bestclu_invmass);
      h1_Zgg_cluster->Fill(bestclu_Zgg);
      h1_opening_angle_cluster->Fill(bestclu_opening_angle);
      h1_dgg_cluster->Fill(bestclu_dgg);
      h2_cluster_dgg_vs_E1pE2->Fill(bestclu_totalE, bestclu_dgg);
      if (best_tower_det_cluster == 0) {
	h1list_mass_by_Ntower[best_tower_id1_cluster]->Fill(bestclu_invmass);
	h1list_mass_by_Ntower[best_tower_id2_cluster]->Fill(bestclu_invmass);
      }
      if (best_tower_det_cluster == 1) {
	h1list_mass_by_Stower[best_tower_id1_cluster]->Fill(bestclu_invmass);
	h1list_mass_by_Stower[best_tower_id2_cluster]->Fill(bestclu_invmass);
      }
      if (bestclu_invmass_Vbbc > -1) {h1_inv_mass_cluster_Vbbc->Fill(bestclu_invmass_Vbbc);}
      if (bestclu_invmass_Vtpc > -1) {h1_inv_mass_cluster_Vtpc->Fill(bestclu_invmass_Vtpc);}
      if (bestclu_invmass_Vvpd > -1) {h1_inv_mass_cluster_Vvpd->Fill(bestclu_invmass_Vvpd);}
      if (bestclu_invmass_Vbbctpc > -1) {h1_inv_mass_cluster_Vbbctpc->Fill(bestclu_invmass_Vbbctpc);}
      if (bestclu_invmass_Vz0tpc > -1) {h1_inv_mass_cluster_Vz0tpc->Fill(bestclu_invmass_Vz0tpc);}
      if (bestclu_invmass_Vvpdtpc > -1) {h1_inv_mass_cluster_Vvpdtpc->Fill(bestclu_invmass_Vvpdtpc);}
    }

    if (check_fillpnt == 1) {
      h1_inv_mass_point->Fill(bestpnt_invmass);
      h1_two_point_energy_allcut->Fill(bestpnt_totalE);
      h2_point_invmass_vs_dgg->Fill(bestpnt_dgg, bestpnt_invmass);
      h2_point_invmass_vs_Zgg->Fill(bestpnt_Zgg, bestpnt_invmass);
      h1_Zgg_point->Fill(bestpnt_Zgg);
      h1_opening_angle_point->Fill(bestpnt_opening_angle);
      h1_dgg_point->Fill(bestpnt_dgg);
      h2_point_dgg_vs_E1pE2->Fill(bestpnt_totalE, bestpnt_dgg);
    }
  }
*/
  return kStOK;
}
