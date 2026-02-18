#include "StMuFcsAnaMakePairs.h"

ClassImp(StMuFcsAnaMakePairs)

StMuFcsAnaMakePairs::StMuFcsAnaMakePairs()
{
  memset(mH1F_InvMassEpdCuts,0,sizeof(mH1F_InvMassEpdCuts));
}

StMuFcsAnaMakePairs::~StMuFcsAnaMakePairs()
{
  delete mH1F_InvMassEpdCuts[0];
  delete mH1F_InvMassEpdCuts[1];  
}

UInt_t StMuFcsAnaMakePairs::LoadHists(TFile* file, HistManager* histman, StMuFcsAnaData* anadata)
{
  UInt_t loaded = 0;
  if( histman==0 ){ return loaded; }

  loaded += histman->AddH2F(file,mH2F_Energy_ph1Vph2,"H2F_Energy_ph1Vph2","Energy of photon 1 vs. photon 2 no EPD cuts;Energy Highest Photon (GeV);Energy Next Highest (GeV)", 1000,0,200, 1000,0,200);
  loaded += histman->AddH1F(file,mH1F_NBadEpdProj,"H1F_NBadEpdProj","Number of points in an event that did not project back to a valid EPD tile;;",15,0,15);
  loaded += histman->AddH1F(file,mH1F_NBadEpdProjVcut,"H1F_NBadEpdProjVcut","Number of points in an event that did not project back to a valid EPD tile with |vertex|<150;;",15,0,15);

  if( mH1F_InvMassEpdCuts[0]==0 ){ mH1F_InvMassEpdCuts[0] = new TObjArray(); }
  loaded += histman->AddH1FArr(file,mH1F_InvMassEpdCuts[0],NEPDCUTS,"H1F_InvMassEpdCuts_AllTrig","Different EPD NMIP cuts all triggers",500,0,1);
  if( mH1F_InvMassEpdCuts[1]==0 ){ mH1F_InvMassEpdCuts[1] = new TObjArray(); }
  loaded += histman->AddH1FArr(file,mH1F_InvMassEpdCuts[1],NEPDCUTS,"H1F_InvMassEpdCuts_EmTrig","Different EPD NMIP cuts EM triggers",500,0,1);
  
  return loaded;
}

Int_t StMuFcsAnaMakePairs::DoMake(StMuFcsAnaData* anadata)
{
  Int_t clustersize = anadata->getEvtInfo()->mClusterSize;
  Int_t npi0candidate = 0;
  TClonesArray* mpharr = anadata->getPhArr();
  TClonesArray* pointpairs = anadata->getPhPairArr();
  Double_t usevertex = anadata->mUseVertex;
  Double_t vertexcutlow = anadata->mVertexCutLow;
  Double_t vertexcuthigh = anadata->mVertexCutHigh;
  //Filling cluster pi0s and cluster photon/elecron epd nmip cut. For clusters only store best pair to speed up code
  for( Int_t ic = 0; ic<clustersize; ++ic ){
    FcsPhotonCandidate* iclus = (FcsPhotonCandidate*) mpharr->UncheckedAt(ic);
    if( !(iclus->mFromCluster) ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
    //std::cout << "|ic:"<<ic << std::endl;
    //std::cout << "  + ";
    //iclus->Print();
    if( ic==(clustersize-1) ){ continue; }

    ///if( iclus->mEpdHitNmip>-0.1){ //Only include candidates who have their nmip value set
    //if( iclus->mEpdHitNmip<mEpdNmipCut ){ goodclusphotonsidx.emplace_back(ic); }
    //else{ goodcluselectronsidx.emplace_back(ic); }
    //}
    for( Int_t jc=ic+1; jc<clustersize; jc++ ){
      FcsPhotonCandidate* jclus = (FcsPhotonCandidate*) mpharr->UncheckedAt(jc);
      if( !(jclus->mFromCluster) ){ std::cout << "MAJOR ERROR - cluster size of array found a point crashing" << std::endl; exit(0); }
      TLorentzVector pi0Vert_LV = iclus->lvVert() + jclus->lvVert();
      if( ic==0 && jc==ic+1 ){ //Since we have a sorted photon array highest two energies are the first two entries
	FcsPi0Candidate* pi0c = (FcsPi0Candidate*) pointpairs->ConstructedAt(npi0candidate++);
	pi0c->mFromCluster = true;
	pi0c->mFromPh = 0;
	pi0c->mPhoton1Idx = ic;
	pi0c->mPhoton2Idx = jc;
	
	pi0c->mPx = pi0Vert_LV.Px();
	pi0c->mPy = pi0Vert_LV.Py();
	pi0c->mPz = pi0Vert_LV.Pz();
	pi0c->mEn = pi0Vert_LV.E();
	
	pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
	pi0c->mDgg     = FcsPi0Candidate::dgg(*iclus,*jclus);
	pi0c->mZgg     = FcsPi0Candidate::zgg(*iclus,*jclus);
	pi0c->mAlpha   = FcsPi0Candidate::alpha(*iclus,*jclus);
	pi0c->mInvMass = pi0Vert_LV.Mag();
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|clustermass:"<<pi0c->mInvMass <<  std::endl;
	break;
      }
      /*
	else{
	//std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
	mH1F_AllPointPairMass->Fill(pi0Vert_LV.Mag());
	}*/
    }
  }

  //Filling point pi0s and cluster photon/elecron epd nmip cut
  //std::cout << "===== EventId:"<< mEvtInfo->mEvent <<" =====" << std::endl;
  FcsPhotonCandidate* firstphotoncut[NEPDCUTS] = {0};
  FcsPhotonCandidate* secondphotoncut[NEPDCUTS] = {0};
  Int_t n_noepdproj = 0;
  Int_t n_noepdproj_vcut = 0;
  for( Int_t ip = clustersize; ip<mpharr->GetEntriesFast(); ++ip ){
    FcsPhotonCandidate* ipoi = (FcsPhotonCandidate*) mpharr->UncheckedAt(ip);
    if( ipoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
    //std::cout << "|ip:"<<ip << std::endl;
    //std::cout << "  + ";
    //ipoi->Print();

    if( ipoi->mEpdHitNmip[0]>-0.1){ //Only include candidates who have their nmip value set
      for( short i=0; i<NEPDCUTS; ++i ){
	if( ipoi->mEpdHitNmip[0] < 0.2+0.1*static_cast<double>(i) ){
	  if( firstphotoncut[i]==0 ){ firstphotoncut[i]=ipoi; }
	  else{ if( secondphotoncut[i]==0 ){ secondphotoncut[i]=ipoi; } }
	}
      }
    }
    else{
      ++n_noepdproj;
      if( vertexcutlow<=usevertex && usevertex<=vertexcuthigh ){ ++n_noepdproj_vcut; }
    }
    
    if( ip==(mpharr->GetEntriesFast()-1) ){ continue; }
    for( Int_t jp=ip+1; jp<mpharr->GetEntriesFast(); ++jp ){
      FcsPhotonCandidate* jpoi = (FcsPhotonCandidate*) mpharr->UncheckedAt(jp);
      if( jpoi->mFromCluster ){ std::cout << "MAJOR ERROR - point size of array found a cluster crashing" << std::endl; exit(0); }
      TLorentzVector pi0Vert_LV = ipoi->lvVert() + jpoi->lvVert();
      //if( ip==clustersize && jp==ip+1 ){ //Since we have a sorted photon array highest two energies are the first two entries
      FcsPi0Candidate* pi0c = (FcsPi0Candidate*) pointpairs->ConstructedAt(npi0candidate++);
      pi0c->mFromCluster = false;
      pi0c->mFromPh = 0;
      pi0c->mPhoton1Idx = ip;
      pi0c->mPhoton2Idx = jp;
      
      pi0c->mPx = pi0Vert_LV.Px();
      pi0c->mPy = pi0Vert_LV.Py();
      pi0c->mPz = pi0Vert_LV.Pz();
      pi0c->mEn = pi0Vert_LV.E();
      
      pi0c->mEta     = pi0Vert_LV.PseudoRapidity();
      pi0c->mDgg     = FcsPi0Candidate::dgg(*ipoi,*jpoi);
      pi0c->mZgg     = FcsPi0Candidate::zgg(*ipoi,*jpoi);
      pi0c->mAlpha   = FcsPi0Candidate::alpha(*ipoi,*jpoi);
      pi0c->mInvMass = pi0Vert_LV.Mag();
      //std::cout << "|idx1:"<<pi0c->mPhoton1Idx << "|idx2:"<<pi0c->mPhoton2Idx << "|pointmass:"<<pi0c->mInvMass <<  std::endl;
      mH2F_Energy_ph1Vph2->Fill(ipoi->mEn,jpoi->mEn);
    }
  }
  mH1F_NBadEpdProj->Fill(n_noepdproj);
  mH1F_NBadEpdProjVcut->Fill(n_noepdproj_vcut);

  //Handle the epd cuts
  for( short i=0; i<NEPDCUTS; ++i ){
    if( firstphotoncut[i]!=0 && secondphotoncut[i]!=0 ){
      TLorentzVector pi0Vert_LV = firstphotoncut[i]->lvVert() + secondphotoncut[i]->lvVert();
      if( anadata->mEmTrigFound ){ ((TH1*)mH1F_InvMassEpdCuts[1]->UncheckedAt(i))->Fill(pi0Vert_LV.Mag()); }
      ((TH1*)mH1F_InvMassEpdCuts[0]->UncheckedAt(i))->Fill(pi0Vert_LV.Mag());
    }
  }
  return kStOk;
}

void StMuFcsAnaMakePairs::PaintEnergy(TCanvas* canv, const char* savename) const
{
  canv->Clear();

  canv->Divide(2,1);

  canv->cd(1)->SetLogz();
  mH2F_Energy_ph1Vph2->Draw("colz");
  canv->cd(2)->SetLogy();
  mH1F_NBadEpdProj->Draw("hist e");
  mH1F_NBadEpdProjVcut->Draw("hist e same");

  canv->Print(savename);
}

void StMuFcsAnaMakePairs::PaintEpdNmipCuts(TCanvas* canv, const char* savename ) const
{
  canv->Clear();
  canv->Divide(2,1);
  for( UInt_t i=0; i<2; ++i ){
    canv->cd(i+1);
    for( Int_t icut=0, ipad=1; icut<NEPDCUTS; ++icut,++ipad ){
      ((TH1*)mH1F_InvMassEpdCuts[i]->UncheckedAt(icut))->SetLineColor(icut+1);//Hack to get rainbow colors since I know I only have 8 histograms
      if( icut==0 ){ ((TH1*)mH1F_InvMassEpdCuts[i]->UncheckedAt(icut))->Draw("hist e"); }
      else{  ((TH1*)mH1F_InvMassEpdCuts[i]->UncheckedAt(icut))->Draw("hist e same"); }
    }
  }
  canv->Print(savename);
}

			  
