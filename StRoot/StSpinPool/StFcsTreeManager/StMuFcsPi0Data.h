/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The Classes related to holding the data in the Pi0TTree in #StMuFcsPi0TreeMaker

  DESCRIPTION
  Contains the classes #FcsEventInfo, #FcsPhotonCandidate, and #FcsPi0Candidate that is used to store information for the Pi0 Transverse Single Spin Asymmetry Analysis (TSSA). The #TTree is created and used in #StMuFcsPi0TreeMaker. Putting these classes in a separate folder was done so that the library can be loaded without loading other STAR libraries. This greatly simplifies the rootmap file needed as well. 

  LOG
  @[September 20, 2024] > Copied from *StMuFcsPi0TreeMaker*

  @[September 21, 2024] > Changed the comment style from Doxygen friendly to ROOT friendly so that the variables show up in the dictionary. This means changing __//!<__  to __///<__ according to [here](https://root.cern/for_developers/doxygen/)

  @[September 26, 2024] > Added #BlueSpin() and #YellowSpin() to #FcsEventInfo to correctly get the blue and yellow beam polarization from the #mSpin. Implemented #FcsEventInfo::Clear() and #FcsEventInfo::Print() as well as #FcsPhotonCandidate::Clear(), #FcsPhotonCandidate::Print(), #FcsPi0Candidate::Clear(), and #FcsPi0Candidate::Print(). Also added comparison functions to #FcsPi0Candidate so that they can be sorted in the #TClonesArray in #StMuFcsPi0TreeMaker. They will be sorted by their distance to the known pi0 mass; e.g. a pi0 candidate with mass 0.14 is less than a candidate with mass 0.1 because 0.14 is closer to the knwon pi0 mass. Added static const #FcsPi0Candidate::Pi0Mass() to just return the mass of the pi0 particle.

*/


#ifndef STMUFCSPI0DATA_HH
#define STMUFCSPI0DATA_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TLorentzVector.h"

class FcsEventInfo : public TObject
{
public:
  FcsEventInfo();
  ~FcsEventInfo();

  Int_t mRunTime = -1;       ///< Time of the run
  Int_t mRunNum = -1;        ///< Run number for event
  UInt_t mFill = 0;          ///< Fill number for event
  UInt_t mEvent = -1;        ///< STAR Event Id
  Int_t mBx48Id = -1;        ///< 48 bit bunch Id for event
  Int_t mBx7Id = -1;         ///< 7 bit bunch Id for event
  UShort_t mSpin = 0;        ///< [Spin bit, this is source polarization](https://drupal.star.bnl.gov/STAR/blog/oleg/spin-patterns-and-polarization-direction)
  //Short_t spinFrom4BitSpin(); ///< Correctly accounts for the spin flip when working with STAR data
  Short_t BlueSpin();        ///< Blue beam polarization at STAR +1 for B+ and -1 for B-
  Short_t YellowSpin();      ///< Yelllow beam polarization at STAR +1 for Y+ and -1 for Y-

  Int_t mTofMultiplicity = -1; ///< TOF Multiplicity
  
  Double_t mVpdVz = -999;      ///< VPD z Vertex
  Double_t mBbcVz = -999;      ///< BBC z Vertex
  Double_t mBbcTacDiff = 0; ///< BBC TAC difference
  Double_t mEpdTacEarlyW = 0;  ///< Earliest EPD TAC for West with cuts 1<adcnmip<15 && TAC>50
  Double_t mEpdTacEarlyE = 0;  ///< Earliest EPD TAC for East with cuts 1<adcnmip<15 && TAC>50
  Double_t mEpdAvgW = 0;    ///< Average EPD TAC for West with cuts 1<adcnmip<15 && TAC>50
  Double_t mEpdAvgE = 0;    ///< Average EPD TAC for East with cuts 1<adcnmip<15 && TAC>50
  //Double_t EpdTacDiffEarly();
  //Double_t EpdTacDiffAvg();
  Double_t mEpdVz = -999;      ///< EPD z Vertex
  Double_t mZdcVz = -999;      ///< ZDC z Vertex
  Short_t mFoundVertex = 0;    ///< Bit vector encoding for which vertex was best; 0 means no vertex, 1=Vpd,2=Epd,4=Bbc

  //This will be used to indicate how many clusters are in the #TClonesArray of #FcsPhotonCandidate. Everything from this number to the size of the array will be points for a given detector Id. I did it this way so I don't have to create a separate branch holding these two numbers and there should only be one #FcsEventInfo object. Also, didn't want a seperate class for clusters and points since they will store the same information. It is kind of a hack since I know that I am only looping up to detector id 2.
  Int_t mClusterSize = 0;       ///< Size of clusters in #mPhArr in #StMuFcsPi0TreeMaker. This means 0 to <#mClusterSize is cluster photon candidates

  virtual void Clear(Option_t* opt="");          ///< Resets all variables to default
  virtual void Print(Option_t* opt="") const;    ///< Prints all values no options

  ClassDef( FcsEventInfo, 1 )
};

//Class to hold basic particle info from which pi0s can be reconstructed
class FcsPhotonCandidate : public TObject
{
public:
  FcsPhotonCandidate();
  ~FcsPhotonCandidate();

  bool mFromCluster = false;  ///< True if from an FCS cluster
  Short_t mDetId = -1;        ///< Detector Id where candidate was found

  Double_t mX = 0;           ///< STAR global x coordinate
  Double_t mY = 0;           ///< STAR global y coordinate
  Double_t mZ = 0;           ///< STAR global z coordinate

  Double_t mEn = 0;          ///< Energy
  Double_t mPxRaw = 0;       ///< X momentum assuming 0,0,0 vertex
  Double_t mPyRaw = 0;       ///< Y momentum assuming 0,0,0 vertex
  Double_t mPzRaw = 0;       ///< Z momentum assuming 0,0,0 vertex

  Double_t mPxVert = 0;       ///< X momentum using best found vertex
  Double_t mPyVert = 0;       ///< Y momentum using best found vertex
  Double_t mPzVert = 0;       ///< Z momentum using best found vertex

  Double_t mEpdHitNmip = 0;   ///< NMIP value from EPD hit

  TLorentzVector lvRaw();        ///< TLorentz vector for this condidate with 0,0,0 vertex momentum
  TLorentzVector lvVert();       ///< TLorentz vector for this candidate with vertex momentum
  Double_t magPosition();        ///< Magnitude of postiion vector i.e. sqrt(#mX^2+#mY^2+#mZ^2)

  virtual void Clear(Option_t* opt="");          ///< Resets all variables to defaults
  virtual void Print(Option_t* opt="") const;    ///< Print all variables no options

  ClassDef( FcsPhotonCandidate, 1 )
};

//Class to hold basic info for reconstructed pi0 candidates
class FcsPi0Candidate : public TObject
{
public:
  FcsPi0Candidate();
  ~FcsPi0Candidate();

  bool mFromCluster = false;  ///< Pi0 reconstructed from clusters or points
  Int_t mPhoton1Idx = -1;     ///< Index in #TClonesArray of #FcsPhotonCandidate 1 that was used to reconstruct this Pi0
  Int_t mPhoton2Idx = -1;     ///< Index in #TClonesArray of #FcsPhotonCandidate 2 that was used to reconstruct this Pi0

  //#StMuFcsPi0TreeMaker will only store information for the Lorentz vector, and other data from the best vertex. To switch to 0,0,0 vertex; use the photon index
  Double_t mPx = 0;            ///< X-Momentum from Lorentz vector of two reconstructed candidates
  Double_t mPy = 0;            ///< Y-Momentum from Lorentz vector of two reconstructed candidates
  Double_t mPz = 0;            ///< Z-Momentum from Lorentz vector of two reconstructed candidates
  Double_t mEn = 0;            ///< Energy from Lorentz vector of two reconstructed candidates
  
  Double_t mEta = -1;          ///< Pseudorapidity from the lorentz vector
  Double_t eta();              ///< Needed since in simulations the stored eta was not there when reconstructing from simulated photons
  Double_t phi();              ///< Angle in x,y plane
  Double_t pt();               ///< Transverse momentum
  Double_t ptot();             ///< Total momentum
  Double_t theta();            ///< azimuthal angle (angle from z-axis to x-y plane)
  Double_t mass();             ///< Invariant mass of the Pi0
  //Need to project using momentum
  //Double_t mStarX = 0;     ///< Global STAR x postion from best vertex
  //Double_t mStarY = 0;     ///< Global STAR y postion from best vertex
  //Double_t mStarZ = 0;     ///< Global STAR z postion from best vertex
  Double_t mDgg = 0;        ///< Euclidean Distance between the two particles (cm) 
  Double_t mZgg = 0;       ///< Energy Asymmetry |E1-E2|/(E1+E2) of pi0
  Double_t mAlpha = 0;     ///< Opening angle of pi0
  Double_t mInvMass = -1;  ///< invariant mass using best vertex as a variable to make it easier for analysis
  
  static Double_t zgg(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2);        ///< Energy asymmetry of pi0
  static Double_t dgg(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2);        ///< Distance between the two particles (cm)
  static Double_t alpha(FcsPhotonCandidate& ph1, FcsPhotonCandidate& ph2);      ///< Opening angle of the two photons

  static const Double_t Pi0Mass(){ return 0.13498; }

  //@[September 26, 2024] > [Need these three functions so that TClonesArray can sort your object](https://root-forum.cern.ch/t/sort-a-tclonesarray/38056)
  virtual Bool_t IsSortable() const {return kTRUE; }  ///< I guess this is flag to indicate to ROOT that object is sortable
  virtual Bool_t IsEqual(const TObject* obj) const;   ///< if both are equal to pi0 mass then return true otherwise false 
  virtual Int_t Compare(const TObject* obj) const;    ///< -1 if distance to pi0 mass of this object is less than the other's distance, 1 if it is greater than, 0 otherwise

  virtual void Clear(Option_t* option="");            ///< Resets all variables to defaults
  virtual void Print(Option_t* option="") const;      ///< Print all variables no options

  ClassDef( FcsPi0Candidate, 1 )
};

#endif
