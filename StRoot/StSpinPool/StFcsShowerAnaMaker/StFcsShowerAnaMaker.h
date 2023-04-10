/*
  The purpose of this StMaker is to understand the shower shape function used by StFcsPointMaker to develop improvements for it
  
  @[April 3, 2023](David Kapukchyan) > First instance
 */

#ifndef StFcsShowerAnaMaker_H
#define StFcsShowerAnaMaker_H

//ROOT headers
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

//STAR headers
#include "StMaker.h"

class StFcsDb;
class StFcsCollection;

class StFcsShowerAnaMaker : public StMaker
{
  
 public:
  StFcsShowerAnaMaker(const char* name="StFcsShowerAnaMaker");
  ~StFcsShowerAnaMaker();
  
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();

  void setFileName(const char* name){ mFileName = name; }    //!< set file name that #Finish() will use to write to a file

  //static void ReadHists(const char* filename); //!< Load from file all the histograms??
  //void WriteHists(const char* filename, const char* mode="RECREATE");      //!< Write all histograms to a file
  
 protected:
  StFcsDb* mFcsDb = 0;               //!< FCS db object
  StFcsCollection* mFcsColl = 0;     //!< FCS collection of hits, clusters, and points

  TString mFileName = "";        //!< name of output file name
  
  TH1F* mH1F_PointXLocal = 0;    //!< Histogram of point local x value but only the decimal value is stored. Local means that deltax=1 is equivalent to 1 cell, e.g. x=3.4 means 3 towers away from beam pipe plus 40% extra into that tower
  TH1F* mH1F_PointYLocal = 0;    //!< Histogram of point local y value but only the decimal value is stored. Local means that deltay=1 is equivalent to 1 cell, e.g. y=2.6 means 2 towers from bottom plus 60% extra into that tower

  TH1F* mH1F_ClusSigMax = 0;  //!< Histogram of cluster sigma max
  TH1F* mH1F_ClusSigMin = 0;  //!< Histogram of cluster sigma min

  TH2F* mH2F_ClusSigMaxEn = 0; //!< Histogram of cluster sigma max vs. cluster energy
  TH2F* mH2F_ClusSigMinEn = 0; //!< Histogram of cluster sigma min vs. cluster energy
  
  void InitHists();     //!< Create histograms from this class. If reading create=false. If writing create=true
  void WriteHists();    //!< Write all histograms to file
  void CleanHists();    //!< Delete all histograms
  
  ClassDef(StFcsShowerAnaMaker,1);
};

#endif
