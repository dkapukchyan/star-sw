/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To generate a ROOT file of histograms related to the FCS and event data in the produced MuDsts from RHIC Run 22 that can be used to do quality assurance (QA) on the data.

  DESCRIPTION
  This class inherits from StMaker and contains many histograms to be used for Quality Assurance (QA) of MuDst files that contain Forward Calorimeter System (FCS) data. It uses #LoadHists() together with the functions #AddH1F(), #AddH2F(), #AddH1FArr(), #AddH2FArr() to ease histogram creation and management. These functions should be used exclusively in this and inherited classes to automate histogram management. #FillEventInfo() is used (and can be re-implemented) to fill QA histograms related to event information. #FillFcsInfo() is used (and can be re-implemented) to fill QA histograms related to FCS information. The idea is that inherited classes don't need to implement the #Init() function only the #LoadHists() function.

  LOG
  @[May 24, 2024 .. June 4, 2024] > Implemented basic variables needed to read data from MuDst. QA histograms of event information: number of events, spin, vertex, trigger, and bunch crossing. Fcs hit information: ADC vs. TB, Energy, Multiplicity, NPeaks, Peak location, Fit Chi^2/NDF, Total Energy. EPD DEP Adc and time peak location to QT adc and time peak location. Fcs Cluster information: Multiplicity (Tower, Neighbor, Points), Energy, Location, SigmaMax and SigmaMin, theta,Chi^2/NDF for 1 and 2 Photon fits. Cluster Pi0 reconstruction with highest energy clusters: Invariant Mass, Angle, Energy Sum, dgg, zgg, High vs. Low Energy. Point information: Multiplicity, Energy, Location. Point Pi0 reconstructions with highest energy points: Invariant Mass, Angle, Energy Sum, dgg, zgg, High vs. Low Energy. Implemented functions #AddH1F(), #AddH2F(), #AddH1FArr(), #AddH2FArr() to ease histogram creation and management. #LoadHists() can be modified by inherited classes to create separate QA histograms. Internal #mAllHists will own and handle all the created histograms automatically when functions #AddH1F(), #AddH2F(), #AddH1FArr(), #AddH2FArr() are used. #FillEventInfo() is used to fill the event information histograms. #FillFcsInfo() fills the others. Also implemented #mSpinRndm to generate random spin patterns for testing, as long #mSpinDbMkr equals 0 it will generate random spins.

  @[June 7, 2024] > Modified the ranges of some histograms. Change #mH1F_Epd_NHits to store total epd hits, and added #mH1F_Epd_NHitsWest to hold information for just the EPD west hits. Added #SetOwner() which will set #mAllHists as owner of all the histograms also changed how histograms are created so that #TObjArray::SetOwner() is only called once as it should be since one call loops over the entire array. Added #mFileOutput for saving ROOT files and is created in #Init() before the histograms to comply with ROOT framework philosophy. In #Make() get the hit, cluster, and point arrays outside the idet loop since each idet call needs the same hit/cluster/point #TClonesArray. Fixed how to loop over hits, clusters, and points since the number of hits/clusters/points doesn't take into the account the first index so need to add the number of hits/clusters/points to the first index to get the correct maximum index of the loop. Added drawing functions for the histograms.

  @[June 25, 2024] > Added booleans #mFcsAdcTbOn, #mEpdAdcQaOn, #mEpdTacQaOn to control turning on and off the histograms for the Fcs ADC vs. TB histograms, the FCS Dep sums vs. Adc from Epd Qa, and the Fcs peak location vs. EPD TAC values respectively. These histograms are all object arrays and will take huge amounts of space and are not essential usually to assess data quality. Added #mH1F_BbcTimeDiff histogram to look at BBC time difference which is used for vertex. Added #mH1F_VertexZdc for checking z vertex from ZDC. Added #mH2F_Mult_tofVecal to check Fcs multiflicity against TOF multiplicity since the reference multiplicity was missing from data. Modified #DrawEventInfo() to plot the new histograms.

  @[July 10, 2024] > Checking if able to retrieve StEpdHits from Trig data. It seems StMuEvent does have trigger data but StEvent does not. The trigger data from MuDsts can be used to fill the EPD hits, which I discovered after running code similar to that in EPD hit maker. modified StEpdHitMaker to read trigger data from MuDsts. Implmented code that will either grab StMuEpdHitCollection from MuDst or from StEpdHitMaker whichever is available, respective of that order. Changed and added some more plot options.

  @[July 15, 2025] > Checking if MuDst contains a #TClonesArray of EPD hits always returns a non-zero value. This means to check if EPD data is in MuDsts it is better to check the size of the #TClonesArray, which is what the code does now.

  Do DEP calib of EPD chs, bunch xing analysis for spin. Change some plots so they use logz and move/remove the stats box for some of hte 2d histograms when plotting. Show on the fly EPD MIP peak locations and valleys
 */

#ifndef STFCSRUN22QAMAKER_HH
#define STFCSRUN22QAMAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TRandom3.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

//STAR Headers
#include "StEnumerations.h"
#include "StMaker.h"
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
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
#include "StEpdDbMaker/StEpdDbMaker.h"
#include "StEpdHitMaker/StEpdHitMaker.h"

class StFcsRun22QaMaker : public StMaker
{
 public:
  StFcsRun22QaMaker(const char* name = "FcsRun22Qa");
  ~StFcsRun22QaMaker();

  virtual Int_t Init();
  virtual Int_t InitRun(int runnumber);
  virtual Int_t Make();
  virtual Int_t Finish();

  void setOutFileName(const char* name){ mFileName = name; }
  //static void LoadDataFromFile(TFile* file, TObjArray* arr);
  void setRandomSeed(ULong_t seed){ mSpinRndm.SetSeed(seed); }
  UInt_t getRandomSeed(){ return mSpinRndm.GetSeed(); }
  Short_t getRandomSpin(); //!< If <0.5 spin down, if >=0.5 spin up
  const char* getFileName(){ return mFileName.Data(); }

  virtual UInt_t LoadHists(TFile* file);
  void SetOwner(Bool_t enable=kTRUE){ mAllHists->SetOwner(enable); }
  static void spinFrom4BitSpin( int spin4bit, int& bpol, int& ypol ); //@[June 3, 2024] > Taken from https://drupal.star.bnl.gov/STAR/blog/oleg/spin-patterns-and-polarization-direction
  void setFcsAdcTbOn(bool value=true)  { mFcsAdcTbOn = value; }
  void setEpdAdcQaOn(bool value=true)  { mEpdAdcQaOn = value; }
  void setEpdTacQaOn(bool value=true)  { mEpdTacQaOn = value; }

  //virtual void Paint(Option_t opt="");
  void DrawEventInfo(TCanvas* canv, const char* savename);

  void DrawVertex(TCanvas* canv, const char* savename);
  void DrawBxId(TCanvas* canv, const char* savename);
  void DrawFcsHitSingle(TCanvas* canv, unsigned int det, const char* savename);
  void DrawFcsTotalE(TCanvas* canv, const char* savename);
  void DrawFcsClusterSingle(TCanvas* canv, unsigned int det, const char* savename);
  void DrawFcsPointSingle(TCanvas* canv, unsigned int det, const char* savename);
  
  
  void DrawAdcVTb(TCanvas* canv, const char* savename);
  void DrawFcsHitQa(TCanvas* canv, const char* savename);
  void DrawEpdHitQa(TCanvas* canv, const char* savename);
  void DrawFcsClusterQa(TCanvas* canv, const char* savename);
  void DrawFcsClusterPi0(TCanvas* canv, const char* savename);
  void DrawFcsPointQa(TCanvas* canv, const char* savename);
  void DrawFcsPointPi0(TCanvas* canv, const char* savename);
  
protected:
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
  //TClonesArray* mMuEpdHits = 0;
  StSpinDbMaker* mSpinDbMkr = 0;     //!< @[May 27, 2024] > Doesn't have proper spin database simply a placeholder
  //StEpdGeom* mEpdGeo=0;
  StEpdHitMaker* mEpdHitMkr = 0;     //@[July 5, 2024] > Checking for afterburner of EPD hits in trig data
  
  //Data to save
  TString mFileName;
  //TFile* mFile_Output = 0; //!< TFile to save all the data

  virtual Int_t FillEventInfo();
  virtual Int_t FillFcsInfo();

  //Machinery to make managing and creating a large number of histograms easier
  UInt_t AddH1F(TFile* file, TH1*& h, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh); //!< This functions should be used to make 1D histograms so that the internal #AllHists obj array can hold a copy to it which will make it easier to write and delete the histograms. Returns 1 if histogram was created/loaded from file, 0 otherwise
  UInt_t AddH1FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh); //!< This functions should be used to make #nobjs number of the same 1D histogram (names will be incremented from 0 to nobjs) and store them in the TObjArray given by #arr. Those histogram pointers will also be copied into #AllHists which will own the object. Returns number of histograms created/loaded from file
  UInt_t AddH2F(TFile* file, TH1*& h, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh);//!< This function should be used to make 2D histograms so that the internal #AllHists obj array can hold a copy to it which will make it easier to write and delete the histograms. Returns 1 if histogram was created/loaded from file, 0 otherwise
  UInt_t AddH2FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh);//!< This functions should be used to make #nobjs number of the same 2D histogram (names will be incremented from 0 to nobjs) and store them in the TObjArray given by #arr. Those histogram pointers will also be copied into #AllHists which will own the object. Returns number of histograms created/loaded from file
  
  TH1* mH1F_Entries = 0;              //!< Number of events processed no cuts (i.e. "Make" calls)
  TH1* mH1F_Triggers = 0;             //!< Triggers in the events
  TH1* mH1F_VertexPrimZ = 0;          //!< Vertex histograms from Primary Vertex
  TH1* mH1F_VertexVpd = 0;            //!< Vertex histograms from VPD
  TH1* mH1F_VertexBbc = 0;            //!< Vertex histograms from BBC
  TH1* mH1F_BbcTimeDiff = 0;          //!< BBC Time difference used to compute the vertex
  TH1* mH1F_VertexZdc = 0;            //!< Vertex from ZDC
  TH1* mH2F_BxId_7V48 = 0;            //!< Bunch crossing Id 7 bit vs. 48 bit
  TH1* mH2F_Mult_tofVref = 0;         //!< Tof multiplicty vs. Reference multiplicity
  TH1* mH2F_Mult_tofVecal = 0;        //!< Tof multiplicity vs. Fcs Ecal multiplicity
  TH1* mH1F_Spin = 0;                 //!< Spin info distribution
  
  TObjArray* mH2F_Hit_adcVtb[kFcsNDet];  //!< Adc vs. tb for all channels
  TH1* mH2F_Hit_enVid[kFcsNDet];         //!< Energy vs. channel Id
  TH1* mH2F_Hit_fitpeakVid[kFcsNDet];    //!< Timebin of peaks vs. channel id
  TH1* mH2F_Hit_chi2Vid[kFcsNDet];       //!< chi^2 of fitted peaks (npeaks>1 only) vs. channel id
  TH1* mH2F_Hit_npeaksVid[kFcsNDet];     //!< number peaks in fit vs. channel id
  TH1* mH1F_Hit_NHits[kFcsNDet];         //!< Hit multiplicity in FCS
  TH1* mH1F_Hit_ESum[3];                 //!< Total energy sum in ecal[0], hcal[1], pres[2]
  //TH1* mH2F_Hit_colVrow[3];            //!< @[May 28, 2024] > Trying to emulate mHitMap in StFcsQaMaker

  TH1* mH1F_Epd_NHits = 0;              //!< Number of hits from EPD collection
  TH1* mH1F_Epd_NHitsWest = 0;          //!< Number of hits from EPD collection only west side
  TObjArray* mH2F_HitPres_depVqt[2];    //!< Special for checking EPD ADC Qt vs. DEP sum split by North[0], South[1] Fcs designation
  TObjArray* mH2F_HitPres_peakVtac[2];  //!< Special for checking EPD TAC values vs. found peak time from FCS split by North[0], South[1] Fcs designation
  
  TH1* mH1F_NClusters[kFcsNDet];              //!< Cluster multiplicity
  TH1* mH1F_Clu_NTowers[kFcsNDet];            //!< Number towers in a cluster
  TH1* mH1F_Clu_NNei[kFcsNDet];               //!< Number of neighbor clusters
  TH1* mH1F_Clu_NPoints[kFcsNDet];            //!< Number points in a cluster
  TH1* mH1F_Clu_En[kFcsNDet];                 //!< Cluster energy
  TH1* mH2F_Clu_yVx[kFcsNDet];                //!< Cluster reconstruction location in local x,y space (i.e. row,column space)
  TH1* mH2F_Clu_sigmaxVsigmin[kFcsNDet];      //!< Cluster sigma max vs. sigma min
  TH1* mH1F_Clu_theta[kFcsNDet];              //!< Cluster theta (angle in x-y plane that defines direction of least second sigma)
  TH1* mH2F_Clu_Chi2NdfPhoton_2V1[kFcsNDet];  //!< Chi^2/NDF for 2 photon fit vs. 1 photon fit
  TH1* mH2F_CluHigh_angleVesum = 0;           //!< opening angle in highest two clusters vs. energy sum of the two clusters
  TH1* mH2F_CluHighEn_lowVhigh = 0;           //!< Highet two cluster energies with highest energy cluster being x-axis and lower one being y-axis
  TH1* mH2F_CluHigh_dggVesum = 0;             //!< Highest two cluster energies, dgg vs. esum
  TH1* mH2F_CluHigh_invmassVesum = 0;         //!< Highest two cluster energies, invariant mass vs. energy sum of the two clusters
  TH1* mH2F_CluHigh_invmassVdgg = 0;          //!< highest 2 clusters dgg vs. invariant mass 
  TH1* mH2F_CluHigh_invmassVzgg = 0;          //!< highest 2 clusters zgg vs. invariant mass

  TH1* mH1F_NPoints[kFcsNDet];          //!< Point Multiplicity
  TH1* mH1F_Poi_En[kFcsNDet];           //!< Point Energy
  TH1* mH1F_Poi_NCluPhotons[kFcsNDet];  //!< number of photons in parent cluster
  TH1* mH2F_Poi_yVx[kFcsNDet];          //!< Point reconstruction location in local x,y space (i.e. row, column space)
  TH1* mH2F_PoiHigh_angleVesum = 0;     //!< point opening angle of two highest points
  TH1* mH2F_PoiHighEn_lowVhigh = 0;     //!< point energy of 2 highest energy points, higher energy point on x-axis, lower energy point on y-axis
  TH1* mH2F_PoiHigh_dggVesum = 0;       //!< 2 highest energy points, dgg vs. Sum of the energy of the 2 points
  TH1* mH2F_PoiHigh_invmassVesum = 0;   //!< Highest two point energies, invariant mass vs. energy sum of the two points
  TH1* mH2F_PoiHigh_invmassVdgg = 0;    //!< highest 2 points
  TH1* mH2F_PoiHigh_invmassVzgg = 0;    //!< highest 2 points

  //TGraph* EpdNmips;                 //!< At first/last event in a run fill with all nMIP values for west side epd channels
  float mEnCut = 1.0;

  bool mFcsAdcTbOn = true;             //!< For turning on/off Adc V tb histograms for the FCS
  bool mEpdAdcQaOn = true;             //!< For turning on/off Qt V Dep histograms from the EPD data
  bool mEpdTacQaOn = true;             //!< For turning on/off Tac V PeakX histograms from the EPD data

private:
  TFile* mFileOutput = 0;              //!< For saving histograms not loading
  TRandom3 mSpinRndm;

  TObjArray* mAllHists = 0;            //!< Store all histogram pointers in this array to make it easier to write and delete them. It will also own all the histograms to make things easier
#ifndef SKIPDefImp
  ClassDef(StFcsRun22QaMaker,2);
#endif
};

#endif
