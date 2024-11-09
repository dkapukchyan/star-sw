/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To read a text file that has trigger names and associated runs and offline trigger Ids for FCS. It will use this map so that you can give it a runnumber and an offline trigger id and it will give you the correct trigger name. It also contains a hard coded map of the trigger thresholds.

  DESCRIPTION
  This class inherits from StMaker to utilize StChain so that it can be found and used in other Makers. It doesn't use the StMaker virtual methods as the function calls since no standard file exists that can be read in and has to be set by the user (Such a file can be found in another repository). It will read a text file with space separated values where the columns are: name, offline trigger id, starting run number, and the ending run number for which that trigger id is valid. From reading the file it generates a vector of unique names. It also genrates a map where the offline trigger id is a key (which should be unique) and the value is a data structure called #TrigRange which holds the four values it read in from the file. With this information utility functions are created to use this information. Also, contains a hard coded map of the transverse energy (Et) (which is equal to transverse momentum (Pt)) threshold for all triggers by name. The map comes from https://www.star.bnl.gov/protected/spin/akio/fcs/run22trg.html.

  LOG
  @[September 13, 2024] > First instance

  @[September 24, 2024] > Changed offline trigger Id in #TrigRange from unsigned int to int and related to changes of all corresponding data to int.

  @[November 25, 2024] > Added Pt==Et trigger thresholds from https://www.star.bnl.gov/protected/spin/akio/fcs/run22trg.html as hard coded values since trigger thresholds did not change during Run 22. This was done by modifying #TrigRange to add hold the Pt trigger threshold value and added get methods for that Pt threshold value. All Pt theresholds are hard coded by name in #SetTriggerPtThresholds() in GeV. Also, made inline comments more ROOT friendly
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

  virtual void readTxtFile(const char* filename);                      ///< Process the text file and fill #mListOfTriggers and #mAllTrigRanges
  int sizeOfTriggers() const;                                          ///< Return size of #mListOfTriggers
  const char* triggerName(int idx) const;                              ///< Return the trigger name in #mListOfTriggers for a given index (#idx)
  const char* nameFromId(Int_t trigidtomatch, Int_t runnumber) const;  ///< Check #mAllTrigRanges for a key #trigidtomatch to find any #TrigRange matches. If not found return "NF" (not found) and if found then check if #runnumber is valid, if so return the name of the trigger; otherwise return "NF".
  Float_t GetPtThr(Int_t trigid) const;                               ///< Get Pt threshold for a given offline trigger id
  Float_t GetPtThrAsym(Int_t trigid) const;                           ///< Get Pt threshold asymetric for a given offline trigger id

protected:
  //Simple struct to hold the tigger ranges when reading the file
  struct TrigRange{
    TrigRange();
    TrigRange(const char* name, Int_t trigid, Int_t startrun, Int_t endrun);
    std::string mName;                 ///< Name of the trigger
    Int_t mOfflineTrigId;              ///< Offline trigger id
    Int_t mStartRun;                   ///< Starting run number (inclusive)
    Int_t mEndRun;                     ///< Ending run number (inclusive)
    Float_t mPtThreshold;              ///< P_T==E_T threshold for trigger
    Float_t mPtThresholdAsym;          ///< Secondary trigger threshold for asymmetric P_T triggers, -1 if it doesn't exist
    bool ValidRun(Int_t runnum) const; ///< check if a given runnum falls into the range
  };

  std::vector<std::string> mListOfTriggers;       ///< All found names of the triggers. This just makes it easier to recall unique trigger namesz
  std::map<Int_t,TrigRange> mAllTrigRanges;       ///< Map of the unique offline trigger Ids to all the found ranges for those triggers

  //const TrigRange& GetTrig(Int_t trigid) const;
  //TrigRange& GetTrig(Int_t trigid) const;

  virtual void SetTriggerPtThresholds(TrigRange& input);  ///< Copied a map of trigger name to its corresponding P_T==E_T thresholds from https://www.star.bnl.gov/protected/spin/akio/fcs/run22trg.html on November 25, 2024. All values will be in GeV. Called when filling #mAllTrigRanges
  
  ClassDef(StFcsRun22TriggerMap,1)
};

#endif

