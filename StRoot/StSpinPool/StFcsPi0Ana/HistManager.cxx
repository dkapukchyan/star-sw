#include "HistManager.h"

ClassImp(HistManager)

HistManager::HistManager():TObjArray()
{
}

HistManager::~HistManager()
{
  delete mOutputFile;
}

TFile* HistManager::InitFile(const char* fname, Option_t* option, const char* ftitle, Int_t compress )
{
  TString name(fname);
  if( name.Length()!=0 && mOutputFile==0 ){ std::cout << "NEWFILE" << std::endl; mOutputFile = new TFile(fname,option,ftitle,compress); }
  return mOutputFile;
}

UInt_t HistManager::AddH1F(TFile* file, TH1*& h1, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh)
{
  UInt_t status = 0;
  if( h1!=0 ){
    //Here I am using bit 22 since that is unused by TH1
    if( h1->TestBit(22) ){ h1=0; } //It is true that the file was loaded so safe to change pointer without delete
    else{ Add(h1); return status; } //Object was new so stop and return 0 and re-add to array
  }
  if( file!=0 ){ h1 = (TH1*)file->Get(name); }
  if( h1==0 ){
    //h1 = (TH1*)FindObject(name);
    //if( h1==0 ){
      h1 = new TH1F(name,title,nbins,xlow,xhigh);
      h1->Sumw2();
      //}
      //else{ return 1; }
  }
  else{
    h1->SetBit(22);
    ++status;
  }
  h1->SetTitle(title);
  Add(h1);
  return status;//1 if histogram loaded or exists, 0 otherwise
}

UInt_t HistManager::AddH1FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbins, Double_t xlow, Double_t xhigh)
{
  UInt_t status = 0;
  for( UInt_t iobj = 0; iobj<nobjs; ++iobj ){
    std::stringstream ss_name;
    ss_name << name << "_" << iobj;
    TH1F* h1 = 0;
    if( file!=0 ){ h1 = (TH1F*)file->Get(ss_name.str().c_str()); }
    if( h1==0 ){
      h1 = new TH1F(ss_name.str().c_str(),title,nbins,xlow,xhigh);
      h1->Sumw2();
    }
    else{
      h1->SetBit(22);
      ++status;
    }
    h1->SetTitle(title);
    arr->Add(h1);
    Add(h1);
  }
  return status;
}

UInt_t HistManager::AddH2F(TFile* file, TH1*& h2, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh)
{
  UInt_t status = 0;
  if( h2!=0 ){
    //Here I am using bit 22 since that is unused by TH1
    if( h2->TestBit(22) ){ h2=0; } //It is true that the file was loaded so safe to change pointer without delete
    else{ Add(h2); return status; } //Object was new so stop and return 0 and re-add to array
  }
  if( file!=0 ){ h2 = (TH1*)file->Get(name); }
  if( h2==0 ){
    //h2 = (TH2*)FindObject(name);
    //if( h2==0 ){

    h2 = new TH2F(name,title, nbinsx,xlow,xhigh, nbinsy,ylow,yhigh);
    h2->Sumw2();
    //}
    //else{ return 1; }
  }
  else{
    h2->SetBit(22);
    ++status;
  }
  h2->SetTitle(title);
  Add(h2);
  return status;//1 if histogram loaded, 0 if new
}

UInt_t HistManager::AddH2FArr(TFile* file, TObjArray*& arr, UInt_t nobjs, const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xhigh, Int_t nbinsy, Double_t ylow, Double_t yhigh)
{
  UInt_t status = 0;
  for( UInt_t iobj = 0; iobj<nobjs; ++iobj ){
    std::stringstream ss_name;
    ss_name << name << "_" << iobj;
    TH2* h2 = 0;
    if( file!=0 ){ h2 = (TH2F*)file->Get(ss_name.str().c_str()); }
    if( h2==0 ){
	h2 = new TH2F(ss_name.str().c_str(),title,nbinsx,xlow,xhigh, nbinsy,ylow,yhigh);
	h2->Sumw2();
    }
    else{
      h2->SetBit(22);
      ++status;
    }
    h2->SetTitle(title);
    arr->Add(h2);
    Add(h2);
  }
  return status; 
}
