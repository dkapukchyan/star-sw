/*
  The purpose of this StMaker is to understand the shower shape function used by StFcsPointMaker to develop improvements for it

  @[April 3, 2023](David Kapukchyan) > First instance
  
  @[April 14, 2023](David Kapukchyan) > Added a energy cut to both clusters and points
  
  @[September 18, 2023](David Kapukchyan) > Added some more histograms and started using the #LoadH1() and #LoadH2() functions in *Rtools* to create/load histograms. This way I can just copy the histograms and the #LoadHistograms() function to easily make plots in another macro. Also added the analysis histograms from 'anaShowerAna.cc' to here.
 */

#ifndef StFcsShowerAnaMaker_H
#define StFcsShowerAnaMaker_H

//ROOT headers
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

//STAR headers
#include "StMaker.h"

//Custom Headers
#include "StFcsPico.h"
#include "ClonesArrTree.h"

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
  void setEnCut( double value )  { mEnCut = value; }

  //static void ReadHists(const char* filename); //!< Load from file all the histograms??
  //void WriteHists(const char* filename, const char* mode="RECREATE");      //!< Write all histograms to a file

  bool LoadHistograms(TObjArray* arr, TFile* file=0);
  
 protected:
  StFcsDb* mFcsDb = 0;               //!< FCS db object
  StFcsCollection* mFcsColl = 0;     //!< FCS collection of hits, clusters, and points

  TString mFileName = "";        //!< name of output file name
  
  TH1* mH1F_PointXLocal = 0;    //!< Histogram of point local x value but only the decimal value is stored. Local means that deltax=1 is equivalent to 1 cell, e.g. x=3.4 means 3 towers away from beam pipe plus 40% extra into that tower
  TH1* mH1F_PointYLocal = 0;    //!< Histogram of point local y value but only the decimal value is stored. Local means that deltay=1 is equivalent to 1 cell, e.g. y=2.6 means 2 towers from bottom plus 60% extra into that tower

  TH1* mH2F_PointXProjX = 0;    //! Histogram of reconstructed point x-value vs. track projected x-value
  TH1* mH2F_PointYProjY = 0;   //! Histogram of reconstructed point y-value vs. track projected y-value

  TH1* mH1F_ClusSigMax = 0;  //!< Histogram of cluster sigma max
  TH1* mH1F_ClusSigMin = 0;  //!< Histogram of cluster sigma min

  TH1* mH2F_ClusSigMaxEn = 0; //!< Histogram of cluster sigma max vs. cluster energy
  TH1* mH2F_ClusSigMinEn = 0; //!< Histogram of cluster sigma min vs. cluster energy

  TH1* mH1F_primid = 0;       //!< Histogram of primary track geant id
  TH1* mH1F_parentid = 0;     //!< Histogram of parent track geant id
  TH1* mH2F_npoiVnclus = 0;   //!< Histogram of number of points vs number of clusters
  TH1* mH2F_cluseVlore = 0;   //!< Histogram of cluster energy vs. cluster lorentz energy (This was to check if cluser->energy() and cluster->fourmomentum()->e() give same result which they do)
  TH1* mH2F_trkeVpoie = 0;    //!< Histogram of primary track energy vs. point energy
  TH1* mH1F_invmasspoi = 0;   //!< Histogram of invariant mass with two highest energy points
  TH1* mH1F_invmasstrk = 0;   //!< Histogram of invariant mass of two highest tracks
  TH1* mH1F_dpoitrk = 0;      //!< Histogram of magnitude of difference between reconstructed point and track
  TH1* mH2F_massVdgg = 0;     //!< Histogram of 2 highest energy points invariant vs. d_{gg} the distance between them
  TH1* mH2F_trkmassVdgg = 0;  //!< Histogram of 2 tracks corresponding to 2 highest energy points invariant mass vs. d_{gg} the distance between them
  
  //void InitHists();     //!< Create histograms from this class. If reading create=false. If writing create=true
  void WriteHists();    //!< Write all histograms to file
  void CleanHists();    //!< Delete all histograms

private:
  double mEnCut = -1;   //!< energy to cut on for both clusters and points

  Rtools::ClonesArrTree* mDataTree = 0; //! Tree to hold information
  TObjArray* mHistsArr = 0;             //! Array to hold histograms
  
  ClassDef(StFcsShowerAnaMaker,1);
};

#endif
