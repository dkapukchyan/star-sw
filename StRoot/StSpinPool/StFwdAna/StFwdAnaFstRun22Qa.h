/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To generate a ROOT file of histograms related to the FST in the produced MuDsts from RHIC Run 22 that can be used to do quality assurance (QA) on the data.

  DESCRIPTION
  This class inherits from StFwdAnaVirtual and contains many histograms to be used for Quality Assurance (QA) of MuDst files that contain the Forward Silicon Tracker (FST) data. It uses #LoadHists() to ease histogram creation and management. #DoMake() is used to fill QA histograms related to event information.

  LOG
  @[July 27, 2026] > First instance that mostly copied from #StFstQAMaker
*/

#ifndef STFWDANA_STFWDANAFSTRUN22QA_HH
#define STFWDANA_STFWDANAFSTRUN22QA_HH

//C/C++ Headers
#include <iostream>
#include <sstream>

//STAR Headers
#include "StEvent/StFstConsts.h"

//Custom headers in this folder
#include "StFwdAnaVirtual.h"

class StFwdAnaFstRun22Qa : public StFwdAnaVirtual
{
 public:
  StFwdAnaFstRun22Qa();
  ~StFwdAnaFstRun22Qa();
  
  virtual UInt_t LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata);
  virtual Int_t DoMake(StFwdAnaData* anadata);

protected:
  //position
  TH1* mH2S_RawHitStrip_rVphi[kFstNumSensors];   ///< raw hit phistrip vs. rstrip per sensor
  TH1* mH2S_HitStripMean_rVphi[kFstNumSensors];  ///< hit mean phistrip vs. mean rstrip per sensor
  TH1* mH2S_Hit_rVphi[kFstNumDisk];              ///< hit map in r vs. phi per disk
  TH1* mH2S_Hit_apvVgeoid[kFstNumDisk];          ///< hit map in APV geometry Id vs. module geometry Id per disk
  TH1* mH2S_HitGlobal_yVx[kFstNumDisk];          ///< hit global x vs. y per disk
  TH1* mH2S_HitGlobal_rVphi[kFstNumDisk];        ///< hit global r vs. phi per disk
  //Charge
  TH1* mH2S_RawHit_adcVgeoid[kFstNumTimeBins];   ///< Charge (ADC) vs channel ID over all time bins
  TH1* mH2S_RawHit_adcerrVgeoid = 0;             ///< RMS noise vs channel ID
  TH1* mH2S_RawHit_maxtbVapv = 0;                ///< Raw hit max ADC time bin vs APV electronics ID [48*(ARC-1)+16*ARM+APV]
  TH1* mH2S_Hit_adcVid = 0;                      ///< Charge vs sensorID
  TH1* mH2S_Hit_adcerrVid = 0;                   ///< Charge uncertainty vs sensorID
  TH1* mH2S_Hit_maxtbVid = 0;	                 ///< hit max ADC time bin vs sensorID
  //hit or raw hit number
  TH1* mH2S_nrawhitsVid = 0;                     ///< number of raw hits vs sensor Id
  TH1* mH2S_nhitsVid = 0;	                 ///< number of hits vs sensor Id
  //TProfile* numOfRawHits_EventId[kFstNumSensors];
  //cluster size
  TH1* mH2S_Hit_nrawhitsVid = 0;                 ///< hit cluster size
  TH1* mH2S_Hit_nrawhitsrVid = 0;                ///< hit cluster size in R direction
  TH1* mH2S_Hit_nrawhitsphiVid = 0;              ///< hit cluster size in Phi direction  

  ClassDef(StFwdAnaFstRun22Qa,1);
};

#endif
  
