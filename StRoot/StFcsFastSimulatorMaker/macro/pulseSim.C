/*
  Macro that uses "FcsPulseSim" to create a pulse in all the modes available.
  
  @[March 29, 2023](David Kapukchyan) > First instance.
*/
  
int pulseSim(const char* savename="testPulse.png", int seed=1)
{
  gROOT->Macro("Load.C");
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gSystem->Load("StEventMaker");
  gSystem->Load("StFcsDbMaker");
  gSystem->Load("StFcsFastSimulatorMaker");
  
  
  StFcsDbPulse* dbpulse = new StFcsDbPulse();
  dbpulse->setTail(2);
  StFcsPulseSim* pulsesim = new StFcsPulseSim();
  pulsesim->setDbPulse(dbpulse); //Needs a StFcsDbPulse object to paramaterize the pulse shape
  pulsesim->setSeed(seed);
  
  TCanvas* c1=new TCanvas("PulseSim","Pulse Sim",1400,800);
  int npadx=3, npady=3;
  c1->Divide(npadx,npady);

  for(int imode=0; imode<StFcsPulseSim::MaxMode; imode++){
    TGraph* G_Temp = pulsesim->pulseSim(imode,100);
    c1->cd(imode+1);
    G_Temp->Draw("APL");
    G_Temp->SetBit(kCanDelete);
  }
  
  c1->SaveAs(savename);

  delete c1;
  delete pulsesim;
  delete dbpulse;

  return 0;
}
