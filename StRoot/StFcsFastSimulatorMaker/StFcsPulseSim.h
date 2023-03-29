//Author:David Kapukchyan
//Date: Aug. 6, 2020
//Class implementation of FCS pulse simulator in 'pulseSim.cpp'

#ifndef STROOT_STFCSFASTSIMULATOR_STFCSPULSESIM_H_
#define STROOT_STFCSFASTSIMULATOR_STFCSPULSESIM_H_

//C/C++ Headers
#include <iostream>
#include <string>

//ROOT Headers
#include "TObject.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TGraphAsymmErrors.h"

//STAR Headers
//#include "StMaker.h"
#include "StFcsDbMaker/StFcsDbPulse.h"

class StFcsPulseSim : public TObject
{
 public:
  //StFcsPulseSim();
  StFcsPulseSim(std::string name="fcspeaksim");
  ~StFcsPulseSim();

  virtual const char* GetName()const{return mName.c_str();}
  std::string GetStrName()const{return mName;}

  void setDbPulse(StFcsDbPulse* p){mDbPulse = p;  mDbPulse->Print(); }
  void setDebug(bool val=true){mDEBUG=val;}//Sets debug option on for printing
  void setSeed(ULong_t seed){mRND->SetSeed(seed);}   //Setting seed
  //Basic Constants
  void setPulseMin(double val){mPulseMin = val;}
  void setPulseMax(double val){mPulseMax = val;}
  void setPed(double val){mPed=val;}
  void setPedSig(double val){mPedSig=val; ZS = int(mPedSig*3);}
  void setMinTb(UShort_t mintb){mMinTb=mintb;}
  void setTBTrg(double val){mTBTrg=val;}
  virtual void InitVars(ULong_t seed);

  //Customize pulse shape
  void setTail(int tail){if( mDbPulse ){mDbPulse->setTail(tail);}}
  
  //TGraph* pulseSim(int mode=0 );//Nine total modes
  TGraphAsymmErrors* pulseSim(int mode=0,double pulseHeight=0);//Nine total modes (pulseHeight only for triggered crossing in certain modes)??
  //static TF1* pulseFunc();//{return TF1* f1 = new TF1(pulseShape());};Maybe make static function so we can reuse for fits??
  static const UShort_t TBMax        = 50;        //Sim Pulse Max timebin
  static const int    MaxMode        = 9;         //Sim Pulse modes

  double TBTrg(){return mTBTrg;}                  //Sim Pulse Peak
  double Ped(){return mPed;}                      //Sim Pulse Pedestal
  double PedSig(){return mPedSig;}                //Sim Pulse Pedestal sigma [adc counts]
  //double BeamLengthSig(){return mBeamLengthSig;}  //[timebin]
  double PulseMin(){return mPulseMin;}            //Sim Pulse min
  double PulseMax(){return mPulseMax;}            //Sim Pulse max

  void setGSigmaSigma(double v){mGSigmaSigma = v;}
  double GSigmaSigma(){return mGSigmaSigma;}

  double Sum(UInt_t StartTb=0, UInt_t EndTb=0);

 protected:
  StFcsDbPulse* mDbPulse = 0;

  //double pulseShapeIntegral(double* x, double* p);
  virtual void addPulse(double pulseTime, double pulseHeight=0);//Add pulse to data arrays with predefined values
  //virtual void addPulse(double pulseTime, double pulseHeight, int method=1)//Newer/Alternate way to add pulse
  virtual void addNoiseDigitize();//Adds noise and converts double data to integers (digital numbers)
  //void Reset();??
  bool mDEBUG;
  TRandom2* mRND;  // Random generator
  double mTBTrg;
  UShort_t mMinTb;
  //UShort_t mMaxTb;

  //generated data
  double ddata[TBMax];
  double idata[TBMax];
  double timebin[TBMax];
  
 private:
  std::string mName;
  
  //Input control variables 
  double mPed;             // pedestal
  double mPedSig;          // pedestal sigma [adc counts]
  //double mBeamLengthSig;   // [timebin]
  //static const UInt_t TBMax          = 50;
  double mPulseMin;
  double mPulseMax;
  //static const int    MaxMode        = 9;
  int    ZS;
  
  double mGSigmaSigma;  // distibution of signal sigma
  
  ClassDef( StFcsPulseSim, 1);

};

#endif

