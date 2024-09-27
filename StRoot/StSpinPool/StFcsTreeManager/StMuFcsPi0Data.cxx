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

Short_t FcsEventInfo::BlueSpin()
{
  if( mSpin==10 || mSpin==9 ){ return 1; }
  else if( mSpin==6 || mSpin==5 ){ return -1; }
  else{ return 0; }
}

Short_t FcsEventInfo::YellowSpin()
{
  if( mSpin==10 || mSpin==6 ){ return 1; }
  else if( mSpin==9 || mSpin==5 ){ return -1; }
  else{ return 0; }
}

void FcsEventInfo::Clear(Option_t* option)
{
  mRunTime = -1;
  mRunNum = -1;
  mFill = 0;
  mEvent = -1;
  mBx48Id = -1;
  mBx7Id = -1;
  mSpin = 0;

  mTofMultiplicity = -1;
  
  mVpdVz = -999;
  mBbcVz = -999;
  mBbcTacDiff = 0;
  mEpdTacEarlyW = 0;
  mEpdTacEarlyE = 0;
  mEpdAvgW = 0;
  mEpdAvgE = 0;
  mEpdVz = -999;
  mZdcVz = -999;
  mFoundVertex = 0;

  mClusterSize = 0;
}

void FcsEventInfo::Print(Option_t* option) const
{
  std::cout <<"## FcsEventInfo|ClusterSize:"<<mClusterSize << std::endl;
  std::cout << " + |RunTime:"<<mRunTime << "|RunNum:"<<mRunNum << "|mFill:"<<mFill << "|mEvent:"<<mEvent << "|Bx48Id:"<<mBx48Id << "|Bx7Id:"<<mBx7Id << "|mSpin:"<<mSpin << "|mTofMult:"<<mTofMultiplicity << std::endl;

  std::cout << " + |VertexBit:0x0"<<std::hex<<mFoundVertex<<std::dec << "|VpdVz:"<<mVpdVz << "|ZdcVz:"<<mZdcVz << "|BbcTDiff:"<<mBbcTacDiff << "|BbcVz:"<<mBbcVz << std::endl;
  std::cout << " + |EpdTacEarlyW:"<<mEpdTacEarlyW << "|EpdTacEarlyE:"<<mEpdTacEarlyE << "|EpdTacAvgW:"<<mEpdAvgW << "|EpdTacAvgE:"<<mEpdAvgE << "|EpdVz:"<<mEpdVz << std::endl;
}

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

void FcsPhotonCandidate::Clear(Option_t* opt)
{
  mFromCluster = false;
  mDetId = -1;
  mX = 0;
  mY = 0;
  mZ = 0;
  
  mEn = 0;
  mPxRaw = 0;
  mPyRaw = 0;
  mPzRaw = 0;
  
  mPxVert = 0;
  mPyVert = 0;
  mPzVert = 0;
  
  mEpdHitNmip = 0;
}

void FcsPhotonCandidate::Print(Option_t* opt) const
{
  std::cout << "|Clus:"<<mFromCluster << "|Pos:("<<mX<<","<<mY<<","<<mZ<<")|En:"<<mEn << "|PRaw:"<<mPxRaw<<","<<mPyRaw<<","<<mPzRaw<<","<<"|PVert:"<<mPxVert<<","<<mPyVert<<","<<mPzVert << std::endl;
}

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

Bool_t FcsPi0Candidate::IsEqual(const TObject* obj) const
{
  FcsPi0Candidate* other = (FcsPi0Candidate*) obj;
  if( obj==0 ){ return kFALSE; }
  //else{ return (this->mInvMass==other->mInvMass); }
  else{ return ( fabs(this->mInvMass-Pi0Mass())==fabs(other->mInvMass-Pi0Mass()) ); }
}

Int_t FcsPi0Candidate::Compare(const TObject* obj) const
{
  FcsPi0Candidate* other = 0;
  other = (FcsPi0Candidate*)obj;
  if( other==0 ){ std::cout << "ERROR - Not an FcsPi0Candidate" << std::endl; return 0; }
  if(      fabs(this->mInvMass-Pi0Mass()) < fabs(other->mInvMass-Pi0Mass()) ){ return -1; }
  else if( fabs(this->mInvMass-Pi0Mass()) > fabs(other->mInvMass-Pi0Mass()) ){ return  1; }
  //if(      (this->mEn) < (other->mEn) ){ return -1; }
  //else if( (this->mEn) > (other->mEn) ){ return  1; }
  else { return 0; }
}

void FcsPi0Candidate::Clear(Option_t* opt)
{
  mFromCluster = false;
  mPhoton1Idx = -1;
  mPhoton2Idx = -1;
  mPx = 0;
  mPy = 0;
  mPz = 0;
  mEn = 0;  
  mEta = -1;
  mDgg = 0;
  mZgg = 0;
  mAlpha = 0;
  mInvMass = -1;  
}

void FcsPi0Candidate::Print(Option_t* opt) const
{
  std::cout << "|Clus:"<<mFromCluster << "|Ph1Idx:"<<mPhoton1Idx << "|Ph2Idx:"<<mPhoton2Idx << "|En:"<<mEn << "|P:("<<mPx<<","<<mPy<<","<<mPz<<")|Eta:"<<mEta << "|Dgg:"<<mDgg << "|Zgg:"<<mZgg << "|Alpha:"<<mAlpha << "|InvMass:"<<mInvMass << std::endl;
}

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


