/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpoe of this class is to match the photon candidates to the EPD hits and properly set the variables of the matched photons. Also contains static functions to generate EPD adjacency tile maps

  DESCRIPTION
  Loops over all the #FcsPhotonCandidate in #StMuFcsAnaData::mPhArr and for each #FcsPhotonCandidate projects the #FcsPhotonCandidate onto the EPD plane; for each projected position it loops over all EPD tiles and checks if this position lies inside any EPD tile; if yes set nMIP to zero and store the tiles usual EPD key (100*pp+tt) into the #FcsPhotonCandidate's data. It then loops over all the EPD hits and sets the nMIP value for the #FcsPhotonCandidate to the matching intersecting EPD tile. The finding of the matched tiles is done in the static function #CheckInsideEpdTile(). It is static so that it can be re-used in other #StMuFcsAna analysis modules.

  LOG
  @[January 22, 2026] > Copied static methods of finding adjacenct EPD tiles from #StMuEpdRun22QaMaker
  @[January 14, 2026] > First instance where relevant functionality was copied from #StMuFcsTreeMaker

*/


#ifndef STMUFCSANAEPDMATCH_HH
#define STMUFCSANAEPDMATCH_HH

#include <map>

#include "StMuFcsVirtualAna.h"

class StMuFcsAnaEpdMatch : public StMuFcsVirtualAna
{
public:
  StMuFcsAnaEpdMatch();
  ~StMuFcsAnaEpdMatch();

  virtual UInt_t LoadHists(TFile* file, HistManager* histman, StMuFcsAnaData* data);
  virtual Int_t DoMake(StMuFcsAnaData* mufcsdata);

  static void CheckInsideEpdTile(StEpdGeom* epdgeo, FcsPhotonCandidate* photon, Double_t projx, Double_t projy);
  
  static std::vector<Int_t> GetAdjacentEpdIds(Int_t pp,Int_t tt);
  static void GetEpdPPandTTFromId(Int_t id, Int_t& pp, Int_t& tt);

  /* Only checked for West EPD
     Chosen so that #TPolyLine can be made to encompass a larger area
     Adjacency functions as a sqaure (CW=ClockWise,CCW=CounterClockWise)
     |----------------------------|
     | outerCCW | outer | outerCW |
     |----------|-------|---------|
     |   CCW    |  main |   CW    |
     |----------|-------|---------|
     | innerCCW | inner | innerCW |
     |----------------------------|
     
     
     X (Beam line)
  */
  static void GetEpdTileOuter(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);       ///! Get the tile going away from beam line from the given tile
  static void GetEpdTileOuterCCW(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);    ///! Get the tile going away and counterclockwise (increasing pp)
  static void GetEpdTileCCW(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);         ///! Get the tile going counterclockwise (increasing pp)
  static void GetEpdTileInnerCCW(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);    ///! Get the tile going in and counterclockwise (increasing pp)
  static void GetEpdTileInner(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);       ///! Get the tile going towards the beam line (in) from the given tile
  static void GetEpdTileInnerCW(Int_t pp, Int_t tt, Int_t &newpp, Int_t& newtt);     ///! Get the tile in and clockwise (decreasing pp)
  static void GetEpdTileCW(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);          ///! Get the tile going clockwise (decreasing pp)
  static void GetEpdTileOuterCW(Int_t pp, Int_t tt, Int_t& newpp, Int_t& newtt);     ///! Get the tile going away and clockwise (decreasing pp)
  
  //void GetAdjacentEpdTile(int pp, int tt, int& pp_adj, int& tt_adj) const;
  static TPolyLine* EpdTilePoly(StEpdGeom* epdgeo, short pp, short tt);       ///< Returns a new polyline using corners from StEpdGeom
  static TPolyLine* EpdCCWOuterCorner(StEpdGeom* epdgeo, short pp, short tt);  ///< Return a new polyline  using adjacency of outer CCW
  static TPolyLine* EpdCWOuterCorner(StEpdGeom* epdgeo, short pp, short tt);
  static TPolyLine* EpdCWInnerCorner(StEpdGeom* epdgeo, short pp, short tt);
  static TPolyLine* EpdCCWInnerCorner(StEpdGeom* epdgeo, short pp, short tt);

  void PaintEpdProjections(TCanvas* canv, const char* savename)   const;
  
protected:
  TH1* mH2F_EpdProjHitMap = 0;          ///< Distribution of x,y projections of photon candidates onto STAR EPD plane in x,y space
  TH1* mH2F_EpdProjHitMap_Vcut = 0;     ///< Distribution of x,y projections of photon candidates onto STAR EPD plane in x,y space with cut |vertex|<150cm
  TH1* mH2F_EpdNmip = 0;                ///< Nmip distributions for EPD matched projected clusters (x-axis bin 1) and points (y-axis bin 2)
  
  static std::map<Int_t,TPolyLine*> mEpdCcwLines;  ///< Map of EPD key to TPolyline to keep track of created polyline objects
  
  ClassDef(StMuFcsAnaEpdMatch,1)
};

#endif

