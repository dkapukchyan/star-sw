#include "StFcsRun22TriggerMap.h"

ClassImp(StFcsRun22TriggerMap)

StFcsRun22TriggerMap::TrigRange::TrigRange()
{
  mName = "";
  mOfflineTrigId = 0;
  mStartRun = 0;
  mEndRun = 0;
}

StFcsRun22TriggerMap::TrigRange::TrigRange(const char* name, UInt_t trigid, Int_t startrun, Int_t endrun)
{
  mName = name;
  mOfflineTrigId = trigid;
  mStartRun = startrun;
  mEndRun = endrun;
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
  }
}

const char* StFcsRun22TriggerMap::nameFromId(UInt_t trigidtomatch, Int_t runnumber) const
{
  auto founditr = mAllTrigRanges.find(trigidtomatch);
  if( founditr==mAllTrigRanges.end() ){ return "NF"; }
  else{
    if( founditr->second.ValidRun(runnumber) ){ return founditr->second.mName.c_str(); }
    else{ return "NF"; }
  }
}

int StFcsRun22TriggerMap::sizeOfTriggers() const
{
  return mListOfTriggers.size();
}

const char* StFcsRun22TriggerMap::triggerName(int idx) const
{
  return mListOfTriggers.at(idx).c_str();
}


