/*
  AUTHOR
  David Kapukchyan

  SYNOPSIS 
  The purpose of this class is to run QA and find Pi0s in FCS data using STAR's MuDst files.

  LOG
  @[January 12, 2023] > Copied from *StFcsPi0FinderForEcal* written by Xilin Liang and adapted to work with STAR's MuDst files.

*/


#ifndef STMUFCSPI0MAKER_HH
#define STMUFCSPI0MAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

//STAR Headers
#include "StEnumerations.h"
#include "StMaker.h"
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
//#include "StEventTypes.h" (not in header)
//#include "StEvent/StEvent.h" (not in header)
//#include "StEvent/StFmsCollection.h" (not in header)
//#include "StEvent/StFmsPoint.h" (not in header)
//#include "StEvent/StFmsPointPair.h" (not in header)
#include "StEvent/StTriggerData.h"
#include "StEvent/StTriggerId.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "Stypes.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFcsHit.h"
#include "StMuDSTMaker/COMMON/StMuFcsCluster.h"
#include "StMuDSTMaker/COMMON/StMuFcsPoint.h"


class StMuFcsPi0Maker : public StMaker {
public:

  StMuFcsPi0Maker(const Char_t* name = "MuFcsPi0Maker");
  ~StMuFcsPi0Maker();
  virtual Int_t InitRun(int runnumber);
  virtual Int_t Make();
  virtual Int_t Finish();
  void setOutFileName(const char* name) { mFilename = name; }
  //{ filename = "StFcsPi0invariantmass"+std::string(run2)+".root"; }
  
private:
  
  //UInt_t mEvent;  //Keeps track of number of events
  //int mTrig = -1;   //Found trigger index in vector 'mTargetTrig'
  //int mXing = 0;    //Bunch Crossing Id
  StMuDstMaker* mMuDstMkr = 0;
  StMuDst* mMuDst = 0;
  StMuEvent* mMuEvent = 0;
  const StTriggerData* mTrigData = 0;
  StRunInfo* mRunInfo = 0;

  StFcsDb* mFcsDb = 0;
  StMuFcsCollection* mMuFcsColl = 0;
  //StFmsEventInfo* mFmsEventInfo;//Since friend functions are not inherited this pointer insures inherited classes still have access to event info
  //StSpinDbMaker* mSpinDbMkr;

  //const char* mOutputName = "fmsMipQa.root";
  //TFile* mFile=0;//July 28, 2021 use pointer in 'mFmsAna_Data'
  //int mPrint=0;
  
  std::vector<int> mTargetTrig;  //For Target Trigger ID
  bool mIgnoreTrig;   //flag to check if ignoring triggers or not
  bool mReadMuDst;    //flag to check if reading from mudst or not (This can be used to turn off populating event info)
  bool mReadSim;     //flag to check if reading mudst from simulations

  //TTree* mDataTree = 0; //An internal tree that you can use to add any desired branches.  It will get written to the output file if it exists.
  //StFmsAna* mFmsAna_Data;//Holds data structures and file to write to
  
   TH1F* h1_num_entries = 0;                //h1_num_entries:# of entries
   TH1F* h1_inv_mass_cluster = 0;           //h1_inv_mass_cluster:invariant mass
   TH1F* h1_inv_mass_cluster_Vtpc = 0;           //h1_inv_mass_cluster_Vtpc:invariant mass with TPC vertex
   TH1F* h1_inv_mass_cluster_Vbbc = 0;           //h1_inv_mass_cluster_Vbbc:invariant mass with BBC vertex
   TH1F* h1_inv_mass_cluster_Vvpd = 0;           //h1_inv_mass_cluster_Vvpd:invariant mass with VPD vertex
   TH1F* h1_inv_mass_cluster_Vz0tpc = 0;           //h1_inv_mass_cluster_Vz0tpc:invariant mass (vertex z = 0) with TPC vertex
   TH1F* h1_inv_mass_cluster_Vbbctpc = 0;           //h1_inv_mass_cluster_Vbbctpc:invariant mass (BBC vertex) with TPC vertex
   TH1F* h1_inv_mass_cluster_Vvpdtpc = 0;           //h1_inv_mass_cluster_Vvpdtpc:invariant mass (VPD vertex) with TPC vertex
   TH1F* h1_Zgg_cluster = 0;                //h1_Zgg:Zgg
   TH1F* h1_opening_angle_cluster = 0;      //h1_opening_angle:opening angle
   TH1F* h1_each_cluster_energy = 0;        //h1_each_cluster_energy:each cluster energy(no cut)
   TH1F* h1_two_cluster_energy_nocut = 0;   //h1_two_cluster_energy_nocut:2 cluster energy(no cut)
   TH1F* h1_two_cluster_energy_allcut = 0;  //h1_two_cluster_energy_allcut:2 cluster energy(all cut)
   TH1F* h1_dgg_cluster = 0;                //h1_dgg_cluster:2 dgg(all cut) distance bewteen two clusters at the det.
   TH1F* h1_Zgg_nocut_cluster = 0;          //h1_Zgg_nocut_cluster:Zgg without cut
   TH1I* h1_nCluster = 0;                   //h1_nCluster: number of clusters
   TH1F* h1_inv_mass_cluster_nocut = 0;     //h1_inv_mass_cluster:invariant mass no cut
   TH1I* h1_nclu_good = 0;                  //h1_nclu_good: number of good clusters
   TH1I* h1_clu_nTowers = 0;                //h1_clu_nTowers: number of towers in cluster

   TH1F* h1_inv_mass_point = 0;           //h1_inv_mass_point:invariant mass
   TH1F* h1_Zgg_point = 0;                //h1_Zgg:Zgg
   TH1F* h1_opening_angle_point = 0;      //h1_opening_angle:opening angle
   TH1F* h1_each_point_energy = 0;        //h1_each_point_energy:each point energy(no cut)
   TH1F* h1_two_point_energy_nocut = 0;   //h1_two_point_energy_nocut:2 point energy(no cut)
   TH1F* h1_two_point_energy_allcut = 0;  //h1_two_point_energy_allcut:2 point energy(all cut)
   TH1F* h1_dgg_point = 0;                //h1_dgg_point:2 dgg(all cut) distance bewteen two points at the det.
   TH1F* h1_Zgg_nocut_point = 0;          //h1_Zgg_nocut_point:Zgg without cut:
   TH1F* h1list_mass_by_Ntower[748];      //h1list_mass_by_Ntower: invariant mass sorted by highest energy tower[64]
   TH1F* h1list_mass_by_Stower[748];      //h1list_mass_by_Stower: invariant mass sorted by highest energy tower[64]
   TH1F* h1list_NEtower[748];             //h1list_NEtower: energy spectrum for north Ecal tower (no cut)
   TH1F* h1list_SEtower[748];             //h1list_SEtower: energy spectrum for south Ecal tower (no cut)
   TH1F* h1_EcalMult_E1 = 0;              //h1_EcalMult_E1: Ecal Milt (E>1)
   TH1I* h1_nPoint = 0;                   //h1_nPoint: number of point
   TH1I* h1_npoi_good = 0;                //h1_npoi_good: number of good points
   TH1F* h1_inv_mass_point_nocut = 0;     //h1_inv_mass_point:invariant mass no cut

   TH1D* h1_zVtpc=0;			//h1_zVtpc: TPC vertex z
   TH1D* h1_zVbbc=0;			//h1_zVbbc: BBC vertex z
   TH1D* h1_zVvpd=0;			//h1_zVvpd: VPD vertex z

   TH2F* h2_EcalMult_vs_TofMult = 0;     //h2_EcalMult_vs_TofMult
   TH2F* h2_cluster_position = 0;        //h2_cluster_position
   TH2F* h2_cluster_position_cut = 0;    //h2_cluster_position_cut
   TH2F* h2_point_position = 0;          //h2_point_position
   TH2F* h2_cluster_invmass_vs_dgg = 0;  //h2_cluster_invmass_vs_dgg
   TH2F* h2_cluster_invmass_vs_Zgg = 0;  //h2_cluster_invmass_vs_Zgg
   TH2F* h2_cluster_dgg_vs_E1pE2 = 0;    //h2_cluster_dgg_vs_E1+E2
   TH2F* h2_point_invmass_vs_dgg = 0;    //h2_cluster_invmass_vs_dgg
   TH2F* h2_point_invmass_vs_Zgg = 0;    //h2_cluster_invmass_vs_Zgg
   TH2F* h2_point_dgg_vs_E1pE2 = 0;      //h2_cluster_dgg_vs_E1+E2

   int mDebug = 0;
   int mFilter = 0;
   int mNEvents = -1;
   int mNAccepted = 0;
   int mMaxEvents = 30000;
   int bins = 150;
   TString mFilename;
   float m_low = 0;
   float m_up = 0.4;
   float E_up = 10;
   float E_low = 8;
   float E_min = 1;

#ifndef SKIPDefImp
  ClassDef(StMuFcsPi0Maker, 1);
#endif
};

#endif
