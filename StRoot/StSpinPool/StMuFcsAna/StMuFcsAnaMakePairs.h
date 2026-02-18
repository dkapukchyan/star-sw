/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpoe of this class is to make pairs of #FcsPhotonCandidate in #StMuFcsAnaData::mPhArr and store them into #StMuFcsAnaData::mPairArr. These pairs serve as the basline for analysis that requires reconstructing a particle by looking at its two decay particles

  DESCRIPTION
  The analysis module loops over all the #FcsPhotonCandidate in #StMuFcsAnaData::mPhArr and creates #FcsPi0Candidate and stores them into #StMuFcsAnaData::mPairArr. It will also check how different EPD nMIP cuts effect the invariant mass of the pairs 

  LOG
  @[January 14, 2026] > First instance where relevant functionality was copied from #StMuFcsTreeMaker

*/


#ifndef STMUFCSANAMAKEPAIRS_HH
#define STMUFCSANAMAKEPAIRS_HH

#include "StMuFcsVirtualAna.h"

class StMuFcsAnaMakePairs : public StMuFcsVirtualAna
{
public:
  StMuFcsAnaMakePairs();
  ~StMuFcsAnaMakePairs();

  virtual UInt_t LoadHists(TFile* file, HistManager* histman, StMuFcsAnaData* anadata);
  virtual Int_t DoMake(StMuFcsAnaData* mufcsdata);

  void PaintEnergy(TCanvas* canv, const char* savename) const;
  void PaintEpdNmipCuts(TCanvas* canv, const char* savename) const;
  
protected:
  TH1* mH2F_Energy_ph1Vph2 = 0;         ///< Histogram of two photons energy used in reconstruction
  TH1* mH1F_NBadEpdProj = 0;            ///< Number of points that did not have a valid projection to an EPD tile in a given event
  TH1* mH1F_NBadEpdProjVcut = 0;        ///< Number of points that did not have a valid projection to an EPD tile in a given event with cut |vertex|<150

  static const short NEPDCUTS = 8;
  TObjArray* mH1F_InvMassEpdCuts[2];    ///< Invariant Mass using different epd nmip cuts and all triggers or only EM triggers
  
  ClassDef(StMuFcsAnaMakePairs,1)
};

#endif

