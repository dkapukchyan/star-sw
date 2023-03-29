#include "StFcsPulseSim.h"

ClassImp(StFcsPulseSim);

StFcsPulseSim::StFcsPulseSim(std::string name)
{
  mName = name;
  mRND = 0;
  mDbPulse = 0;
  //InitVars();
}
/*
StFcsPulseSim::StFcsPulseSim(std::string name, ULong_t seed)
{
  mName = name;
  //InitVars(seed);
}
*/
StFcsPulseSim::~StFcsPulseSim()
{
  delete mRND;
}

void StFcsPulseSim::InitVars(ULong_t seed)
{
  mDEBUG=false;
  memset(ddata,0,sizeof(ddata));
  memset(idata,0,sizeof(idata));
  memset(timebin,0,sizeof(timebin));

  mRND = new TRandom2(seed);

  //mDbPulse = static_cast<StFcsDbPulse*>(StMaker::GetDataSet("fcsPulse"));//Needs to inherit from StMaker to work??
  if( mDbPulse==0 ){std::cout << "WARNING:No FcsDbPulse data found" << std::endl;}
  else{mDbPulse->setTail(2);}//Default for 2020
  
  //Input control variables 
  mGSigmaSigma = 0;  // distibution of signal sigma
  mPed            = 0;   // pedestal
  mPedSig         = 0.8; // pedestal sigma [adc counts]
  //mBeamLengthSig  = 10/mDbPulse->nsecPerTB();   // [timebin]
  //TBMax          = 50;
  mTBTrg          = 25;
  mPulseMin       = 100;
  mPulseMax       = 4000;
  //MaxMode        = 9;
  ZS             = int(mPedSig*3);

  mMinTb = 0;

}

/*
//this one is just shower shape function
double StFcsPulseSim::pulseShape(double* x, double* p)
{
  double ret =  p[0] + p[1]*exp(-0.5*pow((x[0]-p[2])/p[3],2));
  double x1 = x[0] - p[2] - Xoff1;
  if(x1>0){
    double a0 = p[1] * p[3];
    ret += a0*A1/Tau1/Tau1*pow(x1,P1)*exp(-x1/Tau1);
    double x2 = x[0] - p[2] - Xoff2;
    if(x2>0){
      ret += a0*A2/Tau2/Tau2*pow(x2,P2)*exp(-x2/Tau2);
    }
  }
  return ret;
}
*/
/*
TF1* StFcsPulseSim::pulseFunc(Double_t xmin, Double_t xmax)
{
  std::string FuncName = "F1_"+mName+"_PulseFunc";
  TF1* f1 = new TF1(FuncName.c_str(),StFcsPulseSim::pulseShape,xmin,xmax,4);
  return f1;
}
*/
/*
//computing actual integral of the gausian part
double StFcsPulseSim::pulseShapeIntegral(double* x, double* p)
{
  double a0 = p[1] * p[3];
  double ret =  p[0];
  ret += sqrtpi*a0/nsecPerTB/sqrt(2)*
    (TMath::Erf((x[2]-p[2])/p[3]/sqrt(2))
     -TMath::Erf((x[1]-p[2])/p[3]/sqrt(2)));
  //    return ret;
  double x1 = x[0] - p[2] - Xoff1;
  if(x1>0){
    ret += a0*A1/Tau1/Tau1*pow(x1,P1)*exp(-x1/Tau1);
    double x2 = x[0] - p[2] - Xoff2;
    if(x2>0){
      ret += a0*A2/Tau2/Tau2*pow(x2,P2)*exp(-x2/Tau2);
    }
  }
  return ret;
}
*/

void StFcsPulseSim::addPulse(double pulseTime, double pulseHeight)
{
  if( mDEBUG ){printf("AddStart|pulseTime:%4f|pulseHeight:%4f \n",pulseTime,pulseHeight);}
  double t[2];
  double p[4];
  p[3] = mPed;//April 27, 2021 pulseShape in StFcsDbPulse has no pedestal so this access will never happen
  if(pulseHeight==0){
    p[0] = mRND->Uniform(mPulseMin,mPulseMax);
  }else{
    p[0] = pulseHeight;
  }
  p[1] = pulseTime*mDbPulse->nsecPerTB();//Don't need to multiply by nsecPerTB() if pulseTime is timebin value not nanosec. This may be problematic in future for now correct this in macro not here??
  p[2] = mDbPulse->GSigma() + mRND->Uniform(0,mGSigmaSigma);
  if( mDEBUG ){printf("Peak=%4f Time=%4f Sigma=%4.2f Ped=%2.2f\n",p[0],p[1],p[2],p[3]);}
  for(int tb=0; tb<TBMax; tb++){
    t[0] = (tb+0.5)*mDbPulse->nsecPerTB();//See comment above
    double pulse = mDbPulse->pulseShape(t,p)+mPed;
    //printf("%6.2f|%6.1f\n",t[0],pulse);
    //if(tb%20==19) printf("\n");
    ddata[tb] += pulse;
  }
}

/*
void StFcsPulseSim::addPulse(double pulseTime, double pulseHeight, int method)
{
  double t[4];
  double p[4];
  double s1=0,s2=0;
  p[0] = mPed;
  if(pulseHeight==0){
    p[1] = mRND->Uniform(mPulseMin,mPulseMax);
  }else{
    p[1] = pulseHeight;
  }
  p[2] = pulseTime*nsecPerTB;
  p[3] = GSigma + mRND->Uniform(0,mGSigma);
  //res[0]=p[1]*p[3]/nsecPerTB*sqrt(2)*sqrtpi;
  //res[1]=p[2];
  //res[2]=p[3];
  printf("Integral=%10.2f Peak=%8.2f Time=%5.1f Sigma=%4.2f\n",res[0],p[1],res[1],res[2]);
  for(int tb=0; tb<TBMax; tb++){
    t[0] = (tb+0.5)*nsecPerTB;
    t[1] = (tb+0.0)*nsecPerTB;
    t[2] = (tb+1.0)*nsecPerTB;
    double pulse =  pulseShape(t,p);
    double pulse2 = pulseShapeIntegral(t,p);
    s1 += pulse;
    s2 += pulse2;
    //printf("%10.4f %10.4f\n",pulse,pulse2);
    //printf("%6.1f ",pulse); if(tb%20==19) printf("\n");
    if(method==0){
      ddata[tb] += pulse;
    }else{
      ddata[tb] += pulse2;
    }
  }
  //printf("sum(middle point) = %10.1f\n",s1);
  //printf("sum(integral)     = %10.1f\n",s2);
}
*/

void StFcsPulseSim::addNoiseDigitize()
{
  for(int tb=0; tb<TBMax; tb++){
    timebin[tb]=mMinTb+tb;
    ddata[tb] += mRND->Gaus(0,mPedSig);
    int adc = int(ddata[tb]);
    if(adc<=ZS) adc=0;
    if(adc>4000) adc=4000;
    idata[tb]=adc;
    if(mDEBUG){
      printf("%6d ",adc);
      if(tb%20==19) printf("\n");
    }
  }
  if(mDEBUG) printf("\n");	
}


TGraphAsymmErrors* StFcsPulseSim::pulseSim(int mode, double pulseHeight )
{
  double pulseTime = 0;
  memset(ddata,0,sizeof(ddata));
  memset(idata,0,sizeof(idata));
  switch(mode){
  case 0: //single pulse
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    std::cout << "|TBTrg:"<<mTBTrg<<"|pulseTime:" << pulseTime << "|nsecTB:"<<mDbPulse->nsecPerTB() << std::endl;
    addPulse(pulseTime,pulseHeight);
    break;
  case 1: //pulse in pre=1 and triggered xing
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);
    pulseTime = mRND->Gaus(mTBTrg-mDbPulse->TBPerRC(),mDbPulse->BeamLengthSig());
    addPulse(pulseTime);
    break;
  case 2: //pulse in post=1 and triggered xing
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);	    
    pulseTime = mRND->Gaus(mTBTrg+mDbPulse->TBPerRC(),mDbPulse->BeamLengthSig());
    addPulse(pulseTime);	    
    break;
  case 3: //no pulse in triggered xing, one in pre=1
    pulseTime = mRND->Gaus(mTBTrg-mDbPulse->TBPerRC(),mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);	    
    break;
  case 4: //no pulse in triggered xing, one in ppst=1
    pulseTime = mRND->Gaus(mTBTrg+mDbPulse->TBPerRC(),mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);	    
    break;
  case 5: //pulses in triggered, pre1 and post1
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    addPulse(pulseTime);
    pulseTime = mRND->Gaus(mTBTrg-mDbPulse->TBPerRC(),mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);	    
    pulseTime = mRND->Gaus(mTBTrg+mDbPulse->TBPerRC(),mDbPulse->BeamLengthSig());
    addPulse(pulseTime);	    
    break;
  case 6: //single pulse that too big and saturate
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,mRND->Uniform(4500,8000));
    break;
  case 7: //two pulses in triggered xing
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);	    
    pulseTime = mRND->Gaus(mTBTrg,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,pulseHeight);	    
    break;
  case 8: //two pulse in triggered xing, separated by time
    pulseTime = mRND->Gaus(mTBTrg-2,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,mRND->Uniform(mPulseMin,mPulseMax/2.0));
    pulseTime = mRND->Gaus(mTBTrg+2,mDbPulse->BeamLengthSig());
    addPulse(pulseTime,mRND->Uniform(mPulseMin,mPulseMax/2.0));
    break;
  default: 
    std::cout << "Specify mode"<<std::endl;
    return (TGraphAsymmErrors*)0;
  }//switch(mode)

  addNoiseDigitize();

  TGraphAsymmErrors* G_Pulse =  new TGraphAsymmErrors(TBMax,timebin,idata);
  for( UInt_t i=0;i<TBMax;++i)
    {
      StFcsDbPulse::setTGraphAsymmErrors(G_Pulse,i,idata[i],mPedSig,1000.0);
    }
  return G_Pulse;
}

double StFcsPulseSim::Sum(UInt_t StartTb, UInt_t EndTb)
{
  double Total=0;
  if( StartTb==0 && EndTb==0 )
    {
      StartTb = mTBTrg-mDbPulse->TBPerRC();
      EndTb   = mTBTrg+mDbPulse->TBPerRC();
    }
  if( StartTb<0 || EndTb<0 || StartTb>EndTb || StartTb>TBMax || EndTb>TBMax )
    {
      std::cout << "ERROR:Invalid Start timebin("<<StartTb <<") or End timebin("<<EndTb <<") MaxTb("<<TBMax<<"); Return 0" << std::endl;
      return Total;
    }
  for( UInt_t itb=StartTb; itb<=EndTb; ++itb )
    {
      Total += ddata[itb];
    }
  return Total;
}

