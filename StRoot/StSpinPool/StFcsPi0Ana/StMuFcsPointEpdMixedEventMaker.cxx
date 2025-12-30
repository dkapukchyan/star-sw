#include "StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StEventTypes.h"
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

#include "StMuFcsPointEpdMixedEventMaker.h"

ClassImp(StMuFcsPointEpdMixedEventMaker)

StMuFcsPointEpdMixedEventMaker::StMuFcsPointEpdMixedEventMaker(const Char_t* name) : StMuFcsPi0TreeMaker(name)
{
}

StMuFcsPointEpdMixedEventMaker::~StMuFcsPointEpdMixedEventMaker()
{
  delete mMixedPhArr;
}

UInt_t StMuFcsPointEpdMixedEventMaker::LoadHists( TFile* file )
{
  UInt_t loaded = StMuFcsPi0TreeMaker::LoadHists(file);

  //std::cout << "StMuFcsPointEpdMixedEventMaker::LoadHists()" << std::endl;

  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldx,"H2F_PointProj_nmipValldx","nMIP vs. FCS projected point to EPD x minus EPD x of all hits;dX (cm);nmip", 200,-100,100, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldy,"H2F_PointProj_nmipValldy","nMIP vs. FCS projected point to EPD y minus EPD y of all hits;dY (cm);nmip", 200,-100,100, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldr,"H2F_PointProj_nmipValldr","nMIP vs. FCS projected point to EPD, polar distance to all other hits;dR (cm);nmip", 200,-100,100, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipValldphi,"H2F_PointProj_nmipValldphi","nMIP vs. FCS projected point to EPD, angle difference to all other hits;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipValldr,"H2F_MixedPointProj_nmipValldr","Mixed event nMIP vs. FCS projected point to EPD, polar distance to all other hits;dR (cm);nmip", 200,-100,100, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipValldphi,"H2F_MixedPointProj_nmipValldphi","Mixed event nMIP vs. FCS projected point to EPD, angle difference to all other hits;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);

  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledx,"H2F_PointProj_nmipVtiledx","nMIP vs. FCS projected point to EPD x minus EPD x of proj tile;dX (cm);nmip", 80,-20,20, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledy,"H2F_PointProj_nmipVtiledy","nMIP vs. FCS projected point to EPD y minus EPD y of proj tile;dY (cm);nmip", 80,-20,20, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledr,"H2F_PointProj_nmipVtiledr","nMIP vs. FCS projected point to EPD, polar distance to projected hit;dR (cm);nmip", 80,-20,20, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_nmipVtiledphi,"H2F_PointProj_nmipVtiledphi","nMIP vs. FCS projected point to EPD, angle difference to proj tile;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipVtiledr,"H2F_MixedPointProj_nmipVtiledr","Mixed event nMIP vs. FCS projected point to EPD, polar distance to projected hit;dR (cm);nmip", 80,-20,20, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_MixedPointProj_nmipVtiledphi,"H2F_MixedPointProj_nmipVtiledphi","Mixed event nMIP vs. FCS projected point to EPD, angle difference to proj tile;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);

  loaded += mHists->AddH2F(file,mH2F_PointProj_LowMult_nmipValldr,"H2F_PointProj_LowMult_nmipValldr","Low Mult nMIP vs. FCS projected point to EPD, polar distance to all other hits;dR (cm);nmip", 200,-100,100, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_PointProj_LowMult_nmipValldphi,"H2F_PointProj_LowMult_nmipValldphi","Low Mult nMIP vs. FCS projected point to EPD, angle difference to all other hits;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_MixedPointProj_LowMult_nmipValldr,"H2F_MixedPointProj_LowMult_nmipValldr","Low Mult Mixed event nMIP vs. FCS projected point to EPD, polar distance to projected hit;dR (cm);nmip", 200,-100,100, 70,0,7);
  loaded += mHists->AddH2F(file,mH2F_MixedPointProj_LowMult_nmipValldphi,"H2F_MixedPointProj_LowMult_nmipValldphi","Low Mult Mixed event nMIP vs. FCS projected point to EPD, angle difference to proj all;d#phi;nmip", 100,-TMath::Pi(),TMath::Pi(), 70,0,7);

  return loaded;
}

Int_t StMuFcsPointEpdMixedEventMaker::Init()
{
  //This is still needed in Make so even if you are not writing the tree still need these objects
  mMixedPhArr  = new TClonesArray("FcsPhotonCandidate");

  return StMuFcsPi0TreeMaker::Init();
}

Int_t StMuFcsPointEpdMixedEventMaker::InitRun(int runnumber)
{
  return StMuFcsPi0TreeMaker::InitRun(runnumber);
}

//-----------------------
Int_t StMuFcsPointEpdMixedEventMaker::Finish()
{
  return StMuFcsPi0TreeMaker::Finish();
}

//----------------------
Int_t StMuFcsPointEpdMixedEventMaker::Make()
{
  Int_t result = Make_LoadEvent();
  if( result==kStErr ){ return kStErr; }
  result = Make_Polarization();
  if( result == kStErr ){ return kStErr; }
  result = Make_CheckAndSetTrig();
  if( !mValidTrigFound ){ mNTrig=0; return kStSkip; } //Reset trigger array size before going to next event
  if( result!=kStOk ){ return result; }
  result = Make_SpinInfo();
  if( result!=kStOk ){ return result; }

  result = this->Make_GetEpdColl();
  if( result==kStWarn ){ return kStWarn; }

  result = this->Make_VertexInfo();
  if( result!=kStOk ){ return kStErr; }

  result = this->Make_FillFcsClusPoint();
  result = this->Make_CheckAndSetEpdHit();
  //result = this->Make_PointPairs(mPi0Arr);
  //result = this->Make_TssaAna(mPi0Arr);

  if( mPi0Tree!=0 ){ mPi0Tree->Fill(); }
  //Do this here since clear gets called before Make so that MixedPhArr doesn't get filled with junk before first pass
  mNOldPoints = mPhArr->GetEntriesFast();
  mOldVertex = mUseVertex;  //Set old vertex for mixed event analysis on next event
  mMixedPhArr->Clear("C");
  //mMixedPhArr->AbsorbObjects(mPhArr); //@[September 4, 2025] > This may giving a memory leak adding my own for loop to copy
  for( int i=0; i<mPhArr->GetEntriesFast(); ++i ){
    FcsPhotonCandidate* mixedph = (FcsPhotonCandidate*)mMixedPhArr->ConstructedAt(i);
    FcsPhotonCandidate* ph = (FcsPhotonCandidate*)mPhArr->ConstructedAt(i);
    ph->Copy(*mixedph); //Copy function copies ph into mixedph
    //ph->Print();
    //mixedph->Print();
  }

  return kStOk;
}

Int_t StMuFcsPointEpdMixedEventMaker::Make_CheckAndSetEpdHit()
{
  //Check photon candidates if they have any hits in the EPD. Use a separate loop so that this information could be used in the pi0 checking loop if needed. In future may also want to check against FCS preshower (EPD) hits
  Int_t noldhits = mMixedPhArr->GetEntriesFast();
  Int_t nnewhits = mPhArr->GetEntriesFast();
  Int_t npoints = nnewhits - mEvtInfo->mClusterSize;
  Int_t ntotal = noldhits+nnewhits;
  unsigned int nepdhits = 0;
  StSPtrVecEpdHit* epdhits = 0;
  if( mMuEpdHits!=0 ){ nepdhits = mMuEpdHits->GetEntriesFast(); }
  else if( mEpdColl!=0 ){
    epdhits = &(mEpdColl->epdHits());
    nepdhits = epdhits->size();
  }
  else{ LOG_ERROR << "StMuEpdRun22QaMaker::FillEpdinfo() - If you see this error then there is a bug that is setting EPD hits improperly" << endm; return kStErr; }
  Int_t nepdwesthits = 0;
  for( Int_t iph = 0; iph<ntotal; ++iph ){
    //std::cout << "|iph:"<<iph << "|iphnew:"<<iph-noldhits << std::endl;
    FcsPhotonCandidate* ph = 0;
    if( iph>=noldhits ){
      ph = (FcsPhotonCandidate*) mPhArr->UncheckedAt(iph-noldhits);
    }
    else{
      ph = (FcsPhotonCandidate*) mMixedPhArr->UncheckedAt(iph);
    }
    if( ph==0 ){ std::cout << "==========I=CANNOT=BE=ZERO==========" << std::endl; return kStErr; }
    std::vector<Double_t> epdproj;
    if( iph>=noldhits ){
      epdproj = ProjectToEpd(ph->mX,ph->mY,ph->mZ,mUseVertex);
    }
    else{
      //iph<noldhits
      epdproj = ProjectToEpd(ph->mX,ph->mY,ph->mZ,mOldVertex);
    }
    if( iph>=noldhits ){
      //Only fill for the points from current event
      mH2F_EpdProjHitMap->Fill( epdproj.at(0),epdproj.at(1) );
      if( mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh ){ mH2F_EpdProjHitMap_Vcut->Fill(epdproj.at(0),epdproj.at(1)); }
      //std::cout << " + |phx:"<<ph->mX << "|phy:"<<ph->mY << "|phz:"<<ph->mZ << "|v:"<<mUseVertex << std::endl;
      //std::cout << " + |epdx:"<<epdproj.at(0) << "|epdy:"<<epdproj.at(1) << "|epdz:"<<epdproj.at(2) << std::endl;
      //std::cout << " ** |iph:"<< iph-noldhits << std::endl;
      CheckInsideEpdTile(ph,epdproj.at(0),epdproj.at(1));  //Function that will check which EPD tiles photon candidate overlaps with and sets the appropriate variables for it
      /*
      if( ph->mEpdMatch[0]==0 ){ //If no intersection found it would be -1 so now check all the CCW adjacencies
	std::cout << "     ** |iph:"<<iph-noldhits <<"|projx:"<<epdproj.at(0) << "|projy:"<<epdproj.at(1) << "|nmip:"<< ph->mEpdHitNmip[0] << "|epdkey:"<<ph->mEpdMatch[0] << std::endl;
	}*/
    }
    //loop over all hits and if an nmip value exists set for the point
    StMuEpdHit* muepdhit = 0;
    StEpdHit* epdhit = 0;
    for(unsigned int i=0; i<nepdhits; ++i ){
      if( mMuEpdHits!=0 ){ muepdhit = (StMuEpdHit*)mMuEpdHits->UncheckedAt(i); } //To match similar in StMuDstMaker->epdHit(int i)
      else if( epdhits!=0 ){ epdhit = (StEpdHit*)((*epdhits)[i]); }
      else{ LOG_ERROR << "IF YOU SEE THIS ERROR THEN THERE IS A VERY SERIOUS BUG IN THE CODE" << endm; return kStErr; } 
      //std::cout << "|i:"<<i << "|muepdhit:"<<muepdhit << "|epdhit:"<<epdhit << std::endl;
      int ew    = muepdhit!=0 ? muepdhit->side()    : epdhit->side();      //east=-1, west=1
      if( ew==-1 ){ continue; }
      if( iph==noldhits){ ++nepdwesthits; }
      int epdpp = muepdhit!=0 ? muepdhit->position(): epdhit->position();  //Supersector runs [1,12]
      int epdtt = muepdhit!=0 ? muepdhit->tile()    : epdhit->tile();      //Tile number [1,31]
      float nmip = muepdhit!=0 ? muepdhit->nMIP(): epdhit->nMIP();         //The ADC value of the hit divided by the MIP peak position; e.g. if nmip==1 then adc value sits at the MIP peak
      TVector3 epdhitxyz = mEpdGeo->TileCenter(epdpp,epdtt,ew);
      Double_t dx = epdproj.at(0)-epdhitxyz[0];
      Double_t dy = epdproj.at(1)-epdhitxyz[1];
      double rpoint = sqrt(epdproj.at(0)*epdproj.at(0) + epdproj.at(1)*epdproj.at(1));
      double rhit = sqrt(epdhitxyz[0]*epdhitxyz[0] + epdhitxyz[1]*epdhitxyz[1]);
      Double_t phipoint = TMath::ATan2(epdproj.at(1),epdproj.at(0));
      Double_t phihit = TMath::ATan2(epdhitxyz[1],epdhitxyz[0]);
      Double_t diffphi = phipoint-phihit;
      if( diffphi>TMath::Pi() ){ diffphi = diffphi - TMath::Pi(); }
      if( diffphi<(-1.0*TMath::Pi()) ){ diffphi = diffphi + TMath::Pi(); }
      //std::cout << "|epdz:"<<epdhitxyz[2] << std::endl;
      if( ! ph->mFromCluster ){
	//if( mTrigEm2==3 && mTrigEm0<0 && mTrigEm1<0 ){
	  if( iph>=noldhits ){
	    if( mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh ){
	      mH2F_PointProj_nmipValldx->Fill(dx,nmip);
	      mH2F_PointProj_nmipValldy->Fill(dy,nmip);
	      mH2F_PointProj_nmipValldr->Fill(rpoint-rhit,nmip);
	      mH2F_PointProj_nmipValldphi->Fill(diffphi,nmip);
	      if( npoints<=5 ){
		mH2F_PointProj_LowMult_nmipValldr->Fill(rpoint-rhit,nmip);
		mH2F_PointProj_LowMult_nmipValldphi->Fill(diffphi,nmip);
	      }
	    }
	  }
	  else{
	    //iph<noldhits
	    if( mVertexCutLow<=mOldVertex && mOldVertex<=mVertexCutHigh ){
	      mH2F_MixedPointProj_nmipValldr->Fill(rpoint-rhit,nmip);
	      mH2F_MixedPointProj_nmipValldphi->Fill(diffphi,nmip);
	      /*Int_t overflowbin = mH2F_MixedPointProj_nmipValldr->GetBin(201,71);
		if( mH2F_MixedPointProj_nmipValldr->GetBinContent(overflowbin)!=0 ){
		std::cout << "Filled overflow bin:(" << overflowbin << "," << mH2F_MixedPointProj_nmipValldr->GetBinContent(overflowbin)<<")";
		std::cout << "|rpoint:"<<rpoint << "|rhit:"<<rhit << "|nmip:"<<nmip;
		std::cout << "|rpoint-rhit:"<<rpoint-rhit;
		std::cout << "|vert:"<<mOldVertex<< "|("<<epdproj.at(0) << ","<<epdproj.at(1) << ","<<epdproj.at(2)<<")";
		std::cout << std::endl;
		exit(0);
		}*/
	      if( mNOldPoints<=5 ){
		mH2F_MixedPointProj_LowMult_nmipValldr->Fill(rpoint-rhit,nmip);
		mH2F_MixedPointProj_LowMult_nmipValldphi->Fill(diffphi,nmip);
	      }
	    }
	  }
	  //}
      }
      if( ph->mEpdMatch[0] == (100*epdpp+epdtt)  ){  //Check match above
	if( iph>=noldhits ){
	  if( mVertexCutLow<=mUseVertex && mUseVertex<=mVertexCutHigh ){
	    if( ! ph->mFromCluster ){
	      mH2F_PointProj_nmipVtiledx->Fill(dx,nmip);
	      mH2F_PointProj_nmipVtiledy->Fill(dy,nmip);
	      mH2F_PointProj_nmipVtiledr->Fill(rpoint-rhit,nmip);
	      mH2F_PointProj_nmipVtiledphi->Fill(diffphi,nmip);	    
	    }
	  }
	  //int adc = muepdhit!=0 ? muepdhit->adc() : epdhit->adc();
	  ph->mEpdHitNmip[0] = nmip;
	  //std::cout << "|epdpp:"<<epdpp <<"|epdtt:"<<epdtt <<"|nmip:"<<nmip << std::endl;
	}
	else{
	  //iph<noldhits
	  if( mVertexCutLow<=mOldVertex && mOldVertex<=mVertexCutHigh ){
	    if( ! ph->mFromCluster ){
	      mH2F_MixedPointProj_nmipVtiledr->Fill(rpoint-rhit,nmip);
	      mH2F_MixedPointProj_nmipVtiledphi->Fill(diffphi,nmip);
	    }
	  }
	}
      }
    }
    if( iph>=noldhits ){ mH2F_EpdNmip->Fill(ph->mFromCluster,ph->mEpdHitNmip[0]); }
  }
  //std::cout << "|nold:"<<noldhits << "|nnew:"<<nnewhits << "|ntotal:"<<ntotal << "|oldvert:"<<mOldVertex << "|newvert:"<<mUseVertex<< "|nepdhits:"<<nepdwesthits <<"|noldpoints:"<<mNOldPoints << "|npoints:"<<npoints << std::endl;

  //std::cout << "|clustersize:"<<clustersize << "|ncandidates:"<<ncandidates << "|npoints:"<<npoints << std::endl;
  return kStOk;
}

void StMuFcsPointEpdMixedEventMaker::Clear(Option_t* option)
{
  StMuFcsPi0TreeMaker::Clear();
}



void StMuFcsPointEpdMixedEventMaker::PaintEpdAllDistQa(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1)->SetLogz();
  mH2F_PointProj_nmipValldx->Draw("colz");
  mH2F_PointProj_nmipValldx->GetYaxis()->SetRangeUser(0,3);
  canv->cd(4)->SetLogz();
  mH2F_PointProj_nmipValldy->Draw("colz");
  mH2F_PointProj_nmipValldy->GetYaxis()->SetRangeUser(0,3);
  canv->cd(2)->SetLogz();
  mH2F_PointProj_nmipValldr->Draw("colz");
  mH2F_PointProj_nmipValldr->GetYaxis()->SetRangeUser(0,3);
  canv->cd(3)->SetLogz();
  mH2F_PointProj_nmipValldphi->Draw("colz");
  mH2F_PointProj_nmipValldphi->GetYaxis()->SetRangeUser(0,3);

  canv->cd(5)->SetLogz();
  mH2F_MixedPointProj_nmipValldr->Draw("colz");
  mH2F_MixedPointProj_nmipValldr->GetYaxis()->SetRangeUser(0,3);
  canv->cd(6)->SetLogz();
  mH2F_MixedPointProj_nmipValldphi->Draw("colz");
  mH2F_MixedPointProj_nmipValldphi->GetYaxis()->SetRangeUser(0,3);

  canv->Print(savename);
}

void StMuFcsPointEpdMixedEventMaker::PaintEpdDistAnaQa(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  canv->Divide(3,2);

  TH1* pointprojdr = ((TH2*)mH2F_PointProj_nmipValldr)->ProjectionX("pointprojdr",1,mH2F_PointProj_nmipValldr->GetNbinsY(),"e");
  TH1* pointprojdphi = ((TH2*)mH2F_PointProj_nmipValldphi)->ProjectionX("pointprojdphi",1,mH2F_PointProj_nmipValldphi->GetNbinsY(),"e");
  TH1* mixedpointprojdr = ((TH2*)mH2F_MixedPointProj_nmipValldr)->ProjectionX("mixedpointprojdr",1,mH2F_MixedPointProj_nmipValldr->GetNbinsY(),"e");
  TH1* mixedpointprojdphi = ((TH2*)mH2F_MixedPointProj_nmipValldphi)->ProjectionX("mixedpointprojdphi",1,mH2F_MixedPointProj_nmipValldphi->GetNbinsY(),"e");
  canv->cd(1);
  pointprojdr->DrawCopy("hist e");
  canv->cd(4);
  mixedpointprojdr->Draw("hist e");
  canv->cd(2);
  pointprojdphi->DrawCopy("hist e");
  canv->cd(5);
  mixedpointprojdphi->Draw("hist e");
  
  pointprojdr->Divide(mixedpointprojdr);
  pointprojdphi->Divide(mixedpointprojdphi);
  
  canv->cd(3);
  pointprojdr->Draw("hist e");
  pointprojdr->GetYaxis()->SetRangeUser(0.95,1.15);
  canv->cd(6);
  pointprojdphi->Draw("hist e");

  canv->Print(savename);
  
}

void StMuFcsPointEpdMixedEventMaker::PaintEpdAllDistQaLowMult(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  canv->Divide(3,2);

  canv->cd(1)->SetLogz();
  mH2F_PointProj_LowMult_nmipValldr->Draw("colz");
  mH2F_PointProj_LowMult_nmipValldr->GetYaxis()->SetRangeUser(0,3);
  canv->cd(2)->SetLogz();
  mH2F_PointProj_LowMult_nmipValldphi->Draw("colz");
  mH2F_PointProj_LowMult_nmipValldphi->GetYaxis()->SetRangeUser(0,3);
  canv->cd(4)->SetLogz();
  mH2F_MixedPointProj_LowMult_nmipValldr->Draw("colz");
  mH2F_MixedPointProj_LowMult_nmipValldr->GetYaxis()->SetRangeUser(0,3);
  canv->cd(5)->SetLogz();
  mH2F_MixedPointProj_LowMult_nmipValldphi->Draw("colz");
  mH2F_MixedPointProj_LowMult_nmipValldphi->GetYaxis()->SetRangeUser(0,3);

  TH1* pointprojlowmultdr = ((TH2*)mH2F_PointProj_LowMult_nmipValldr)->ProjectionX("pointprojlowmultdr",1,mH2F_PointProj_LowMult_nmipValldr->GetNbinsY(),"e");
  TH1* pointprojlowmultdphi = ((TH2*)mH2F_PointProj_LowMult_nmipValldphi)->ProjectionX("pointprojlowmultdphi",1,mH2F_PointProj_LowMult_nmipValldphi->GetNbinsY(),"e");
  TH1* mixedpointprojlowmultdr = ((TH2*)mH2F_MixedPointProj_LowMult_nmipValldr)->ProjectionX("mixedpointprojlowmultdr",1,mH2F_MixedPointProj_LowMult_nmipValldr->GetNbinsY(),"e");
  TH1* mixedpointprojlowmultdphi = ((TH2*)mH2F_MixedPointProj_LowMult_nmipValldphi)->ProjectionX("mixedpointprojlowmultdphi",1,mH2F_MixedPointProj_LowMult_nmipValldphi->GetNbinsY(),"e");
  /*
  canv->cd(5);
  pointprojdr->DrawCopy("hist e");
  canv->cd(6);
  mixedpointprojdr->Draw("hist e");
  canv->cd(7);
  pointprojdphi->DrawCopy("hist e");
  canv->cd(8);
  mixedpointprojdphi->Draw("hist e");
  */

  pointprojlowmultdr->Divide(mixedpointprojlowmultdr);
  pointprojlowmultdphi->Divide(mixedpointprojlowmultdphi);
  
  canv->cd(3);
  pointprojlowmultdr->Draw("hist e");
  pointprojlowmultdr->GetYaxis()->SetRangeUser(0.85,0.95);
  canv->cd(6);
  pointprojlowmultdphi->Draw("hist e");


  canv->Print(savename);

}

void StMuFcsPointEpdMixedEventMaker::PaintEpdTileDistQa(TCanvas* canv, const char* savename) const
{
  canv->Clear();
  canv->Divide(3,2);
  canv->cd(1)->SetLogz();
  mH2F_PointProj_nmipVtiledx->Draw("colz");
  mH2F_PointProj_nmipVtiledx->GetXaxis()->SetRangeUser(-20,20);
  mH2F_PointProj_nmipVtiledx->GetYaxis()->SetRangeUser(0,3);
  canv->cd(4)->SetLogz();
  mH2F_PointProj_nmipVtiledy->Draw("colz");
  mH2F_PointProj_nmipVtiledy->GetXaxis()->SetRangeUser(-20,20);
  mH2F_PointProj_nmipVtiledy->GetYaxis()->SetRangeUser(0,3);
  canv->cd(2)->SetLogz();
  mH2F_PointProj_nmipVtiledr->Draw("colz");
  mH2F_PointProj_nmipVtiledr->GetXaxis()->SetRangeUser(-20,20);
  mH2F_PointProj_nmipVtiledr->GetYaxis()->SetRangeUser(0,3);
  canv->cd(3)->SetLogz();
  mH2F_PointProj_nmipVtiledphi->Draw("colz");
  mH2F_PointProj_nmipVtiledphi->GetXaxis()->SetRangeUser(-0.3,0.3);
  mH2F_PointProj_nmipVtiledphi->GetYaxis()->SetRangeUser(0,3);
  canv->cd(5)->SetLogz();
  mH2F_MixedPointProj_nmipVtiledr->Draw("colz");
  mH2F_MixedPointProj_nmipVtiledr->GetXaxis()->SetRangeUser(-20,20);
  mH2F_MixedPointProj_nmipVtiledr->GetYaxis()->SetRangeUser(0,3);
  canv->cd(6)->SetLogz();
  mH2F_MixedPointProj_nmipVtiledphi->Draw("colz");
  mH2F_MixedPointProj_nmipVtiledphi->GetXaxis()->SetRangeUser(-0.3,0.3);
  mH2F_MixedPointProj_nmipVtiledphi->GetYaxis()->SetRangeUser(0,3);

  canv->Print(savename);
}


