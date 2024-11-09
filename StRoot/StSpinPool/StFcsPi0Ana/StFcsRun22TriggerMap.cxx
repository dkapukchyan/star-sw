#include "StFcsRun22TriggerMap.h"

ClassImp(StFcsRun22TriggerMap)

StFcsRun22TriggerMap::TrigRange::TrigRange()
{
  mName = "";
  mOfflineTrigId = 0;
  mStartRun = 0;
  mEndRun = 0;
  mPtThreshold = 0;
  mPtThresholdAsym = -1;
}

StFcsRun22TriggerMap::TrigRange::TrigRange(const char* name, Int_t trigid, Int_t startrun, Int_t endrun)
{
  mName = name;
  mOfflineTrigId = trigid;
  mStartRun = startrun;
  mEndRun = endrun;
  mPtThreshold = 0;
  mPtThresholdAsym = -1;
}

bool StFcsRun22TriggerMap::TrigRange::ValidRun(Int_t runnum) const
{
  if( mStartRun<=runnum && runnum<=mEndRun ){ return true; }
  else{ return false; }
}


StFcsRun22TriggerMap::StFcsRun22TriggerMap(const char* name):StMaker(name)
{}

StFcsRun22TriggerMap::~StFcsRun22TriggerMap()
{}

void StFcsRun22TriggerMap::readTxtFile(const char* filename)
{
  std::fstream infile(filename);
  if( !infile.is_open() ){ LOG_ERROR << "StFcsRun22TriggerMap::readTxtFile - Unable to open file" << endm; }
  while( !infile.eof() ){
    std::string name;
    infile >> name;
    if( name.size()==0 ){ continue; }
    bool namefound = false;
    for( unsigned int i=0; i<mListOfTriggers.size(); ++i ){
      if( mListOfTriggers.at(i) == name ){ namefound = true; }
    }
    if( !namefound ){ mListOfTriggers.push_back(name); }
    UInt_t trigid;
    Int_t startrun;
    Int_t endrun;
    infile >> trigid >> startrun >> endrun;
    //@[September 16, 2024] > To use emplace just create the object on the fly as hinted [here](https://cplusplus.com/forum/beginner/283109/)
    //                         + A better solution without neccessary object creation [here](https://www.machinet.net/tutorial-eng/how-to-use-emplace-in-cpp-maps)
    //                         + Same method emphasized [here](https://stackoverflow.com/questions/68645539/how-to-use-emplace-in-map-for-custom-class)
    //                         + and [here](https://stackoverflow.com/questions/6162201/c11-use-case-for-piecewise-construct-of-pair-and-tuple)
    auto inserted = mAllTrigRanges.emplace(std::piecewise_construct,
					   std::forward_as_tuple(trigid),
					   std::forward_as_tuple(name.c_str(),trigid,startrun,endrun)
					   );
    if( !inserted.second ){ LOG_WARN << "StFcsRun22TriggerMap::readTxtFile - Trigger Ids not unique found at least two with ID:"<<trigid << endm; }
    else{ SetTriggerPtThresholds( (*(inserted.first)).second ); }
  }
}

const char* StFcsRun22TriggerMap::nameFromId(Int_t trigidtomatch, Int_t runnumber) const
{
  auto founditr = mAllTrigRanges.find(trigidtomatch);
  if( founditr==mAllTrigRanges.end() ){ return "NF"; }
  else{
    if( founditr->second.ValidRun(runnumber) ){ return founditr->second.mName.c_str(); }
    else{ return "NF"; }
  }
}

void StFcsRun22TriggerMap::SetTriggerPtThresholds(TrigRange& input)
{
  //List order matches the table on the website https://www.star.bnl.gov/protected/spin/akio/fcs/run22trg.html
  if( input.mName=="fcsJP2" )    { input.mPtThreshold = 8; return; }
  if( input.mName=="fcsJPA1" )   { input.mPtThreshold = 6; return; }
  if( input.mName=="fcsJPA0" )   { input.mPtThreshold = 4; return; }
  if( input.mName=="fcsJPBC1" )  { input.mPtThreshold = 6; return; }
  if( input.mName=="fcsJPBC0" )  { input.mPtThreshold = 4; return; }
  if( input.mName=="fcsJPDE1" )  { input.mPtThreshold = 6; return; }
  if( input.mName=="fcsJPDE0" )  { input.mPtThreshold = 4; return; }
  if( input.mName=="fcsDiJP" )   { input.mPtThreshold = 5; input.mPtThresholdAsym = 5; return; }
  if( input.mName=="fcsDiJPAsy" ){ input.mPtThreshold = 5; input.mPtThresholdAsym = 4; return; }
  if( input.mName=="fcsDY" )     { input.mPtThreshold = 1; return; }
  if( input.mName=="fcsJPsi" )   { input.mPtThreshold = 0.7; return; }
  //if( input.mName=="fcsDYNoEpd" ){ input.mPtThreshold = 1; return; } //Not in trigger map
  if( input.mName=="fcsDYAsy" )  { input.mPtThreshold = 1; input.mPtThresholdAsym=0.7; return; }
  if( input.mName=="fcsHad2" )   { input.mPtThreshold = 6; return; }
  if( input.mName=="fcsHad1" )   { input.mPtThreshold = 4; return; }
  if( input.mName=="fcsHad0" )   { input.mPtThreshold = 2; return; }
  if( input.mName=="fcsEM2" )    { input.mPtThreshold = 5; return; }
  if( input.mName=="fcsEM1" )    { input.mPtThreshold = 3.5; return; }
  if( input.mName=="fcsEM0" )    { input.mPtThreshold = 2; return; }
  if( input.mName=="fcsELE2" )   { input.mPtThreshold = 1; return; }
  if( input.mName=="fcsEM3" )    { input.mPtThreshold = 1; return; }
  //TPC triggers match the non TPC counterpart
  if( input.mName=="fcsEM2_tpc" )    { input.mPtThreshold = 5; return; }
  if( input.mName=="fcsEM1_tpc" )    { input.mPtThreshold = 3.5; return; }
  if( input.mName=="fcsEM0_tpc" )    { input.mPtThreshold = 2; return; }

  //Debug triggers or triggers unable to find on website
  if( input.mName=="fcs_led" )       { input.mPtThreshold = 0; return; }
  if( input.mName=="fcs_HcalTot" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="fcs_EcalTot" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS7-HTOT-S" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS6-ETOT-S" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS5-HHT-S" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS4-EHT-S" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS3-EM3-S" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS2-ELE3-S" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS1-ELE1-S" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSS0-ELE0-S" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN7-HTOT-N" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN6-ETOT-N" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN5-HHT-N" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN4-EHT-N" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN3-EM3-N" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN2-ELE2-N" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN1-ELE1-N" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCSN0-ELE0-N" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsHTOT-S" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsHTOT-N" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsHHT-S" )      { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsHHT-N" )      { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsETOT-S" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsETOT-N" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsEHT-S" )      { input.mPtThreshold = 0; return; }
  if( input.mName=="fcsEHT-N" )      { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS15-DiELEA" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS15_DiELEA" )  { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS14-DiJPAsy" ) { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS13-DiJP" )    { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS12-JPDE0" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS11-JPBC0" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS10-JPA0" )    { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS09-JPDE1" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS08-JPBC1" )   { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS07-JPA1" )    { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS06-JP2" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS05-EM2" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS04-EM1" )     { input.mPtThreshold = 0; return; }
  if( input.mName=="FCS03-EM0" )     { input.mPtThreshold = 0; return; }

  LOG_WARN << "Unable to find trigger P_T==E_T threshold for '" << input.mName << "'" << endm;
  return;
}

int StFcsRun22TriggerMap::sizeOfTriggers() const
{
  return mListOfTriggers.size();
}

const char* StFcsRun22TriggerMap::triggerName(int idx) const
{
  return mListOfTriggers.at(idx).c_str();
}

/*
const TrigRange& StFcsRun22TriggerMap::GetTrig(Int_t trigid) const
{
  auto itr=mAllTrigRanges.find(trigid);
  if( itr!=mAllTrigRanges.end() ){ return itr->second; }
  else{ return 0; }
}

TrigRange& StFcsRun22TriggerMap::GetTrig(Int_t trigid) const
{
  auto itr=mAllTrigRanges.find(trigid);
  if( itr!=mAllTrigRanges.end() ){ return itr->second; }
  else{ return 0; }
}
*/

Float_t StFcsRun22TriggerMap::GetPtThr(Int_t trigid) const
{
  auto itr=mAllTrigRanges.find(trigid);
  if( itr!=mAllTrigRanges.end() ){ return itr->second.mPtThreshold; }
  else{ return 999; } //if no trigger id exists then return a large value to avoid including bad triggers
}

Float_t StFcsRun22TriggerMap::GetPtThrAsym(Int_t trigid) const
{
  auto itr=mAllTrigRanges.find(trigid);
  if( itr!=mAllTrigRanges.end() ){ return itr->second.mPtThresholdAsym; }
  else{ return 999; } //if no trigger id exists then return a large value to avoid including bad triggers
}


