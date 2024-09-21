#include "StMuFcsPi0Data.h"

ClassImp(FcsEventInfo);

FcsEventInfo::FcsEventInfo()
{}

FcsEventInfo::~FcsEventInfo()
{}

/*Not sure if needed below code is not right
Double_t FcsEventInfo::EpdTacDiffEarly()const
{
  if( mEpdEarlyW>50 && mEpdEarlyE>0 ){ return mEpdEarlyW-mEpdEarlyE;}
  else{ return -999; }
}
Double_t FcsEventInfo::EpdTacDiffAvg()const{ return mEpdAvgW-mEpdAvgE; }  
*/

ClassImp(FcsPhotonCandidate)

FcsPhotonCandidate::FcsPhotonCandidate()
{}

FcsPhotonCandidate::~FcsPhotonCandidate()
{}

TLorentzVector FcsPhotonCandidate::lvRaw()
{
  TLorentzVector v;
  v.SetPxPyPzE(mPxRaw,mPyRaw,mPzRaw,mEn);
  return v;
}

TLorentzVector FcsPhotonCandidate::lvVert()
{
  TLorentzVector v;
  v.SetPxPyPzE(mPxVert,mPyVert,mPzVert,mEn);
  return v;
}

Double_t FcsPhotonCandidate::magPosition()
{ return sqrt( mX*mX + mY*mY + mZ*mZ ); }

ClassImp(FcsPi0Candidate)

FcsPi0Candidate::FcsPi0Candidate()
{}

FcsPi0Candidate::~FcsPi0Candidate()
{}

Double_t FcsPi0Candidate::eta()
{ if( mEta<0 ){ return asinh(mPz/pt()); }else{ return mEta; } }

Double_t FcsPi0Candidate::phi()
{ return atan2(mPy,mPx); }

Double_t FcsPi0Candidate::pt()
{ return sqrt( mPx*mPx + mPy*mPy ); }

Double_t FcsPi0Candidate::ptot()
{ return sqrt( mPx*mPx + mPy*mPy + mPz*mPz ); }

Double_t FcsPi0Candidate::theta()
{ return 2.0*atan(exp(-1.0*eta())); }

Double_t FcsPi0Candidate::mass()
{ return sqrt(mEn*mEn - ptot()*ptot()); }

Double_t FcsPi0Candidate::zgg(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2)
{return fabs(ph1.mEn-ph2.mEn)/(ph1.mEn+ph2.mEn);}

Double_t FcsPi0Candidate::dgg(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2)
{ return sqrt( (ph1.mX-ph2.mX)*(ph1.mX-ph2.mX) + (ph1.mY-ph2.mY)*(ph1.mY-ph2.mY) + (ph1.mZ-ph2.mZ)*(ph1.mZ-ph2.mZ) ); }

Double_t FcsPi0Candidate::alpha(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2)
{
  Double_t ph1dotph2 = ph1.mX*ph2.mX + ph1.mY*ph2.mY + ph1.mZ*ph2.mZ; //dot product of vectors for the current cluster position and cluster j position
  Double_t ph1mag = ph1.magPosition(); //magnitude of position vector for candidate 1
  Double_t ph2mag = ph2.magPosition(); //magnitude of position vector for candidate 2
  return acos( ph1dotph2 / (ph1mag*ph2mag) );
}


