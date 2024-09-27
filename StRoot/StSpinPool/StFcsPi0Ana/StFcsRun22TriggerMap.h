/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To read a text file that has trigger names and associated runs and offline trigger Ids for FCS. It will use this map so that you can give it a runnumber and an offline trigger id and it will give you the correct trigger name.

  DESCRIPTION
  This class inherits from StMaker to utilize StChain so that it can be found and used in other Makers. It doesn't use the StMaker virtual methods as the function calls since no standard file exists that can be read in and has to be set by the user (Such a file can be found in another repository). It will read a text file with space separated values where the columns are: name, offline trigger id, starting run number, and the ending run number for which that trigger id is valid. From reading the file it generates a vector of unique names. It also genrates a map where the offline trigger id is a key (which should be unique) and the value is a data structure called #TrigRange which holds the four values it read in from the file. With this information utility functions are created to use this information

  LOG
  @[September 13, 2024] > First instance

  @[September 24, 2024] > Changed offline trigger Id in #TrigRange from unsigned int to int and related to changes of all corresponding data to int.
 */

#ifndef STFCSRUN22TRIGGERMAP_HH
#define STFCSRUN22TRIGGERMAP_HH

#include <fstream>

#include "StMessMgr.h"
#include "StMaker.h"

class StFcsRun22TriggerMap : public StMaker
{
 public:
  StFcsRun22TriggerMap(const char* name="fcsRun22TrigMap");
  virtual ~StFcsRun22TriggerMap();

  virtual void readTxtFile(const char* filename);                      //!< Process the text file and fill #mListOfTriggers and #mAllTrigRanges
  int sizeOfTriggers() const;                                          //!< Return size of #mListOfTriggers
  const char* triggerName(int idx) const;                              //!< Return the trigger name in #mListOfTriggers for a given index (#idx)
  const char* nameFromId(Int_t trigidtomatch, Int_t runnumber) const; //!< Check #mAllTrigRanges for a key #trigidtomatch to find any #TrigRange matches. If not found return "NF" (not found) and if found then check if #runnumber is valid, if so return the name of the trigger; otherwise return "NF".

protected:
  //Simple struct to hold the tigger ranges when reading the file
  struct TrigRange{
    TrigRange();
    TrigRange(const char* name, Int_t trigid, Int_t startrun, Int_t endrun);
    std::string mName;                 //!< Name of the trigger
    Int_t mOfflineTrigId;             //!< Offline trigger id
    Int_t mStartRun;                   //!< Starting run number (inclusive)
    Int_t mEndRun;                     //!< Ending run number (inclusive)
    bool ValidRun(Int_t runnum) const; //!< check if a given runnum falls into the range
  };

  std::vector<std::string> mListOfTriggers;       //!< All found names of the triggers. This just makes it easier to recall unique trigger namesz
  std::map<Int_t,TrigRange> mAllTrigRanges;      //!< Map of the unique offline trigger Ids to all the found ranges for those triggers
  
  ClassDef(StFcsRun22TriggerMap,1)
};

#endif

