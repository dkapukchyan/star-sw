#include "StFwdAnaFstRun22Qa.h"

StFwdAnaFstRun22Qa::StFwdAnaFstRun22Qa()
{
  memset(mH2S_RawHitStrip_rVphi,  0, sizeof(mH2S_RawHitStrip_rVphi));
  memset(mH2S_HitStripMean_rVphi, 0, sizeof(mH2S_HitStripMean_rVphi));
  memset(mH2S_Hit_rVphi,          0, sizeof(mH2S_Hit_rVphi));
  memset(mH2S_Hit_apvVgeoid,      0, sizeof(mH2S_Hit_apvVgeoid));
  memset(mH2S_HitGlobal_yVx,      0, sizeof(mH2S_HitGlobal_yVx));
  memset(mH2S_HitGlobal_rVphi,    0, sizeof(mH2S_HitGlobal_rVphi));

  memset(mH2S_RawHit_adcVgeoid,   0, sizeof(mH2S_RawHit_adcVgeoid));
}

StFwdAnaFstRun22Qa::~StFwdAnaFstRun22Qa()
{}

UInt_t StFwdAnaFstRun22Qa::LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata)
{
  UInt_t nloaded = 0;
  std::stringstream ss_histname;
  std::stringstream ss_histtitle;
  for( int iwedge=0; iwedge<kFstNumWedges; ++iwedge ){
    for( int isensor=0; isensor<kFstNumSensorsPerWedge; ++isensor ){
      ss_histname.str("");
      ss_histtitle.str("");
      
      ss_histname << "H2S_RawHitStrip_rVphi_"<< iwedge*3+isensor+1;
      ss_histtitle << "Raw Hit rstrip vs. #phi strip: Wedge "<<iwedge+1<< " Sensor " <<  isensor << ";#phi Strip;r Strip";
      nloaded += histman->AddH2S(file,mH2S_RawHitStrip_rVphi[iwedge*3+isensor],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumPhiSegPerWedge,0,kFstNumPhiSegPerWedge, kFstNumRStripsPerWedge,0,kFstNumRStripsPerWedge);

      //ss_histname.str("");
      //ss_histtitle.str("");
      //sprintf(buffer,"numOfRawHitsVsEventId_Sensor%d", iwedge*3+isensor+1);
      //sprintf(histname, "Number of raw hits vs. EventID: Wedge %d Sensor %d", iwedge+1, isensor);
      //numOfRawHits_EventId[iwedge*3+isensor] = new TProfile(buffer, histname, 10000, 0, 10000);
      //numOfRawHits_EventId[iwedge*3+isensor]->GetXaxis()->SetTitle("EventID (Time)");
      //numOfRawHits_EventId[iwedge*3+isensor]->GetYaxis()->SetTitle("<rawhits>");

      ss_histname.str("");
      ss_histtitle.str("");
      ss_histname << "H2S_HitStripMean_rVphi_"<<iwedge*3+isensor+1;
      ss_histtitle << "Hit mean rstrip vs. #phi strip: Wedge "<<iwedge+1<< " Sensor " << isensor <<";Mean #phi strip;Mean r strip";
      nloaded += histman->AddH2S(file,mH2S_HitStripMean_rVphi[iwedge*3+isensor],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumPhiSegPerWedge,0,kFstNumPhiSegPerWedge, kFstNumRStripsPerWedge,0,kFstNumRStripsPerWedge);
    }		
  }
  
  for( int idisk = 0; idisk<kFstNumDisk; ++idisk ){
    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_Hit_rVphi_Disk"<<idisk+1;
    ss_histtitle << "FST Hit Map for Disk "<<idisk+1 << " r vs #phi;#phi strip;r strip";
    nloaded += histman->AddH2S(file,mH2S_Hit_rVphi[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumPhiSegPerWedge*kFstNumWedgePerDisk,0,kFstNumPhiSegPerWedge*kFstNumWedgePerDisk, kFstNumRStripsPerWedge,0,kFstNumRStripsPerWedge);

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_Hit_apvVgeoid_Disk"<<idisk;
    ss_histtitle << "FST hit map in APV vs Wedge Idx for Disk "<<idisk <<";WedgeIdx [wedge-(disk-1)*WedgePerDisk+1];APV Geometry ID";
    nloaded += histman->AddH2S(file,mH2S_Hit_apvVgeoid[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(), kFstNumWedgePerDisk,1,kFstNumWedgePerDisk+1, kFstApvsPerWedge,0,kFstApvsPerWedge);

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_HitGlobal_yVx_Disk"<<idisk;
    ss_histtitle << "FST Hit Map in STAR coordinates y vs. x for Disk "<<idisk+1 <<";STAR X (cm);STAR Y(cm)";
    nloaded += histman->AddH2S(file,mH2S_HitGlobal_yVx[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(), 140,-35,35, 140,-35,35);

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_HitGlobal_rVphi_Disk"<<idisk;
    ss_histtitle << "FST Hit Map in STAR coordinates r vs. #phi for Disk "<<idisk+1 <<";STAR #phi (rad.);STAR r (cm)";
    nloaded += histman->AddH2S(file,mH2S_HitGlobal_rVphi[idisk],ss_histname.str().c_str(),ss_histtitle.str().c_str(),kFstNumPhiSegPerWedge*kFstNumWedgePerDisk/8, -TMath::Pi(), TMath::Pi(), kFstNumRStripsPerWedge, 5, 28);
  }

  for( int itb=0; itb<kFstNumTimeBins; ++itb ){
    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_RawHit_adcVgeoid_tb"<<itb;
    ss_histtitle << "ADC of raw hits at timebin "<<itb << "vs. channel geometry ID;Channel Geometry ID;ADC of Raw Hits";
    nloaded += histman->AddH2S(file,mH2S_RawHit_adcVgeoid[itb],ss_histname.str().c_str(),ss_histtitle.str().c_str(),288, 0, 36864, 512, 0, kFstMaxAdc);    
  }

    nloaded += histman->AddH2S(file,mH2S_RawHit_adcerrVgeoid,"H2S_RawHit_adcerrVgeoid","RMS noise of raw hits vs. channel geometry ID;Channel Geometry ID;RMS noise of Raw Hits",288, 0, 36864, 128, 0, 64);
    nloaded += histman->AddH2S(file,mH2S_RawHit_maxtbVapv,"H2S_Hit_maxtbVapv","Max time bin of hit vs. APV ID;APV ID [48*(RDO-1)+16*ARM+APV];Max Time Bin Index",kFstNumApvs, 0, kFstNumApvs, kFstNumTimeBins, 0, kFstNumTimeBins);
    nloaded += histman->AddH2S(file,mH2S_Hit_adcVid,"H2S_Hit_adcVid","ADC of hits vs. Sensor ID;Sensor ID;ADC of Hits",kFstNumSensors,0,kFstNumSensors, 512,0,kFstMaxAdc);
    nloaded += histman->AddH2S(file,mH2S_Hit_adcerrVid,"H2S_Hit_adcerrVid","RMS noise of hits vs Sensor ID;Sensor ID;RMS noise of Hits", kFstNumSensors,0,kFstNumSensors, 128,0,64);
    nloaded += histman->AddH2S(file,mH2S_Hit_maxtbVid,"H2S_Hit_maxtbVid","Max time bin of hits vs. Sensor ID;Sensor ID;Max Time Bin Index",kFstNumSensors,0,kFstNumSensors, kFstNumTimeBins,0,kFstNumTimeBins);
  
  nloaded += histman->AddH2S(file,mH2S_nrawhitsVid,"H2S_nrawhitsVid","Number of Raw Hits vs. Sensor ID;Sensor ID;Number of Raw Hits", kFstNumSensors,0,kFstNumSensors, 128,0,128);
  nloaded += histman->AddH2S(file,mH2S_nhitsVid,"H2S_nhitsVid","The number of hits vs. Sensor ID;Sensor ID,Number of Hits", kFstNumSensors,0,kFstNumSensors, 128, 0, 128);
  
  nloaded += histman->AddH2S(file,mH2S_Hit_nrawhitsVid,"H2S_Hit_nrawhitsVid","Number of Raw Hits in a Cluster(Hit) vs. Sensor ID;Sensor ID;Number of Raw Hits in Cluster",kFstNumSensors,0,kFstNumSensors, 20,0,20);
  nloaded += histman->AddH2S(file,mH2S_Hit_nrawhitsrVid,"H2S_Hit_nrawhitsrVid","Hit size in R of raw hits vs Sensor ID;Sensor ID;Hit Size in R of Raw Hits",kFstNumSensors,0,kFstNumSensors, 20,0,20);
  nloaded += histman->AddH2S(file,mH2S_Hit_nrawhitsphiVid,"H2S_Hit_nrawhitsphiVid","Hit size in #phi of raw hits vs. Sensor ID;Sensor ID;Hit Size in #phi of Raw Hits", kFstNumSensors,0,kFstNumSensors, 20,0,20);

  return nloaded;
}


Int_t StFwdAnaFstRun22Qa::DoMake(StFwdAnaData* anadata)
{
  //StMuFstCollection* mufstcoll = anadata->muFstColl();
  TClonesArray* tc_mufstraw = anadata->muFstRawHitColl();
  TClonesArray* tc_mufsthit = anadata->muFstHitColl();
  if( tc_mufstraw==0 ){
    LOG_WARN << "StFwdAnaFstRun22Qa::No FST Raw Hit Collection" << endm;
    return kStWarn;
  }
  if( tc_mufsthit==0 ){
    LOG_WARN << "StFwdAnaFstRun22Qa::No FST Hit Collection" << endm;
    return kStWarn;
  }

  //loop over raw hits
  unsigned int rawcounter[kFstNumSensors]; //raw hit multiplicity per sensor per event
  memset(rawcounter,0,sizeof(rawcounter));
  for( int iraw = 0; iraw<tc_mufstraw->GetEntriesFast(); ++iraw ){
    StMuFstRawHit* rawhit = (StMuFstRawHit*)tc_mufstraw->At(iraw);

    //Cast to unsigned int to avoid bad char to int conversions
    unsigned int wedge = rawhit->getWedge();
    unsigned int sensor = rawhit->getSensor();
    unsigned int maxTimeBin = rawhit->getMaxTimeBin();
    unsigned int sensorId = (wedge-1)*static_cast<unsigned int>(kFstNumSensorsPerWedge) + sensor;
    ++rawcounter[sensorId];

    for( unsigned char timeBin = 0; timeBin < kFstNumTimeBins; ++timeBin ) {
      mH2S_RawHit_adcVgeoid[timeBin]->Fill(rawhit->getGeoId(), (int)rawhit->getCharge( timeBin ));
    }
    mH2S_RawHit_adcerrVgeoid->Fill(rawhit->getGeoId(), (int)(rawhit->getChargeErr( maxTimeBin )+0.5));
    mH2S_RawHitStrip_rVphi[sensorId]->Fill((int)rawhit->getPhiStrip(), (int)rawhit->getRStrip());
    mH2S_RawHit_maxtbVapv->Fill(((int)rawhit->getRdo()-1)*48+(int)rawhit->getArm()*16+(int)rawhit->getApv(), (int)maxTimeBin);

  }

  unsigned int hitcounter[kFstNumSensors];  //hit multiplicity per sensor
  memset(hitcounter,0,sizeof(hitcounter));
  for( int ihit = 0; ihit<tc_mufsthit->GetEntriesFast(); ++ihit ){
    StMuFstHit* hit = (StMuFstHit*)tc_mufsthit->At(ihit);
    if( hit!=0 ){
      const TVector3& hitpos = hit->xyz();
      int sensorIdxTemp = ((int)hit->getWedge()-1)*kFstNumSensorsPerWedge + (int)hit->getSensor(); // 0-107
      ++hitcounter[sensorIdxTemp];
      int diskIdxTemp = ((int)hit->getWedge()-1)/kFstNumWedgePerDisk + 1; // 1-3
      int wedgeIdxTemp = (int)hit->getWedge() - (diskIdxTemp-1)*kFstNumWedgePerDisk; // 1-12
      int phiIdxTemp = (wedgeIdxTemp-1)*kFstNumPhiSegPerWedge+(int)(hit->getMeanPhiStrip()+0.5);
      int rIdxTemp = (int)(hit->getMeanRStrip()+0.5);

      mH2S_HitStripMean_rVphi[sensorIdxTemp]->Fill((int)(hit->getMeanPhiStrip()+0.5), (int)(hit->getMeanRStrip()+0.5));
      mH2S_Hit_rVphi[diskIdxTemp-1]->Fill(phiIdxTemp, rIdxTemp);
      mH2S_Hit_apvVgeoid[diskIdxTemp-1]->Fill(wedgeIdxTemp, (int)hit->getApv()+(wedgeIdxTemp%2-1)*kFstApvsPerWedge);
      mH2S_HitGlobal_yVx[diskIdxTemp-1]->Fill((float)hitpos.X(), (float)hitpos.Y());
      mH2S_HitGlobal_rVphi[diskIdxTemp-1]->Fill((float)hitpos.Phi(), (float)hitpos.Perp());

      mH2S_Hit_adcVid->Fill(sensorIdxTemp, (int)hit->getCharge());
      mH2S_Hit_adcerrVid->Fill(sensorIdxTemp, (int)(hit->getChargeErr()+0.5));
      mH2S_Hit_maxtbVid->Fill(sensorIdxTemp, (int)hit->getMaxTimeBin());
      mH2S_Hit_nrawhitsVid->Fill(sensorIdxTemp, (int)hit->getNRawHits());
      mH2S_Hit_nrawhitsrVid->Fill(sensorIdxTemp, (int)hit->getNRawHitsR());
      mH2S_Hit_nrawhitsphiVid->Fill(sensorIdxTemp, (int)hit->getNRawHitsPhi());
    }
  }

  for( int iS=0; iS<kFstNumSensors; ++iS ){
    mH2S_nrawhitsVid->Fill(iS, rawcounter[iS]);
    mH2S_nhitsVid->Fill(iS, hitcounter[iS]);
    //numOfRawHits_EventId[iS]->Fill((int)eventPtr->id()/100+1, counter[iS]);
  }
  
  return kStOk;
}


