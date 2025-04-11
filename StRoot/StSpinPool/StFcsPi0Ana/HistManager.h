/*
  AUTHOR
  David Kapukchyan
  
  PURPOSE
  To hold and manage all created histograms when running a ROOT macro so that reading from a file, writing to a file and deleting histograms at the end of a macro can happen in one line. It is best to treat this collection as the owners of the histogram and leave the creation and destruction of histograms to the collection.

  DESCRIPTION
  This class inherits from #TObjArray hence is a kind of #TObjArray. Histograms should be added to the array with the corresponding *AddH* functions. Also, has an internal file pointer because creating histograms before a #TFile leads to "bad"[1] behavior when reading and writing to files. The idea is that you pass a file pointer and a histogram pointer into the function and it will check if the histogram can be loaded from the file first; if so then set histogram pointer to the histogram from the file. If not, then create a new histogram (with the rest of the arguments) and add it to the array and set histogram pointer to that new historam. This way you can have one function whose job is to load histograms and you can re use that to both write and read the same histograms.

  [1]The "bad" behavior mostly has to do with how ROOT works in assigning histograms to a #TDirectory. If no file exists no directory is set and thus writing and reading don't work as expected.
  
  LOG
  @[August 20, 2024] > First instance. It was copied from the relevant parts of #StMuFcsRun22QaMaker.
  @[September 9, 2024] > Added #getFileName()
  @[December 9, 2024] > Added #Draw() to help in drawing a specific named histogram in the collection.
  @[January 4, 2025] > Changed how histograms are created and deleted to address a bug found when working with multiple #TFiles that would delete the histograms on closing so that the status bits were rendered unusable. The collection will now delete pointers rather than reassigning them. Also, fixed how #AddH1FArr() and #AddH2FArr() work so that they delete the histograms
  @[January 29, 2025] > Added #AddH2F() using variable y bins. Added #AddH3F() for creating 3D histograms with variable bin sizes
  @[February 6, 2025] > Added #AddH1D() for 1D double histograms. Note that when histograms are not found when loading from multiple files, It doesn't properly clean things up for the next file
 */

#ifndef HISTMANAGER_HH
#define HISTMANAGER_HH

#include <iostream>
#include <sstream>

//ROOT Headers
#include "Compression.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


class HistManager : public TObjArray
{
public:
  HistManager();
  virtual ~HistManager();
  const char* getFileName(){ if( mOutputFile!=0 ){ return mOutputFile->GetName(); }else{return 0;} }
  
  TFile* InitFile(const char* fname = "", Option_t* option = "", const char* ftitle="", Int_t compress=101 ); //!< initialize the file to store all the histograms in
  //Using 101 since *ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault* is not found when including Compression.h and that is what the value is equal to according to ROOT documentation
  
  //Machinery to make managing and creating a large number of histograms easier
  UInt_t AddH1F(TFile* file, TH1*& h, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh); //!< This functions should be used to make 1D histograms so that the internal #AllHists obj array can hold a copy to it which will make it easier to write and delete the histograms. Returns 1 if histogram was created/loaded from file, 0 otherwise
  UInt_t AddH1D(TFile* file, TH1*& h, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh); //!< This functions should be used to make 1D histograms so that the internal #AllHists obj array can hold a copy to it which will make it easier to write and delete the histograms. Returns 1 if histogram was created/loaded from file, 0 otherwise
  UInt_t AddH1FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh); //!< This functions should be used to make #nobjs number of the same 1D histogram (names will be incremented from 0 to nobjs) and store them in the TObjArray given by #arr. Those histogram pointers will also be copied into #AllHists which will own the object. Returns number of histograms created/loaded from file
  UInt_t AddH2F(TFile* file, TH1*& h, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh);//!< This function should be used to make 2D histograms so that the internal #AllHists obj array can hold a copy to it which will make it easier to write and delete the histograms. Returns 1 if histogram was created/loaded from file, 0 otherwise
  UInt_t AddH2F(TFile* file, TH1*& h, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t* ybins);//!< This function should be used to make 2D histograms with variable bins so that the internal #AllHists obj array can hold a copy to it which will make it easier to write and delete the histograms. Returns 1 if histogram was created/loaded from file, 0 otherwise
  UInt_t AddH2FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh);//!< This functions should be used to make #nobjs number of the same 2D histogram (names will be incremented from 0 to nobjs) and store them in the TObjArray given by #arr. Those histogram pointers will also be copied into #AllHists which will own the object. Returns number of histograms created/loaded from file
  UInt_t AddH3F(TFile* file,  TH1*& h, const char* name, const char* title,  Int_t nbinsx, Double_t* xbins, Int_t nbinsy, Double_t* ybins, Int_t nbinsz, Double_t *zbins);

  virtual void Draw(Option_t* option=""); ///< Plot a single histogram by its name. Add drawing option using semicolon

private:
  TFile* mOutputFile = 0;  //!< For using a common file to store all the histgrams

  ClassDef(HistManager,1)
};

#endif //HISTMANAGER_HH

