#ifndef HcalSimAlgos_HcalTimingBit_h
#define HcalSimAlgos_HcalTimingBit_h

#include <bitset>
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalUpgradeTriggerPrimitiveDigi.h"

class HcalTimingBit {
 public:
  
 HcalTimingBit(int version_LLPflag) : version_LLPflag_(version_LLPflag) {};

  std::bitset<6> compute(int ibin, HcalUpgradeTriggerPrimitiveDigi& TPDigi) const;
  
 private:
  int version_LLPflag_;
  // define prompt-delayed TDC range. Note this is offset from depth and ieta by 1
  const int tdc_HE[29][7] = {{8, 14, 15, 17, 0, 0, 0},
                             {8, 14, 15, 17, 0, 0, 0},
                             {8, 14, 14, 17, 0, 0, 0},
                             {8, 14, 14, 17, 0, 0, 0},
                             {8, 13, 14, 16, 0, 0, 0},
                             {8, 13, 14, 16, 0, 0, 0},
                             {8, 12, 14, 15, 0, 0, 0},
                             {8, 12, 14, 15, 0, 0, 0},
                             {7, 12, 13, 15, 0, 0, 0},
                             {7, 12, 13, 15, 0, 0, 0},
                             {7, 12, 13, 15, 0, 0, 0},
                             {7, 12, 13, 15, 0, 0, 0},
                             {7, 11, 12, 14, 0, 0, 0},
                             {7, 11, 12, 14, 0, 0, 0},
                             {7, 11, 12, 14, 0, 0, 0},
                             {7, 11, 12, 7, 0, 0, 0},
                             {0, 12, 10, 0, 0, 0, 0},
                             {0, 9, 10, 9, 10, 0, 0},
                             {16, 9, 9, 9, 11, 10, 0},
                             {17, 8, 9, 8, 9, 10, 0},
                             {9, 7, 7, 7, 9, 6, 0},
                             {8, 7, 7, 6, 6, 6, 0},
                             {8, 6, 6, 6, 7, 7, 0},
                             {7, 6, 6, 6, 7, 6, 0},
                             {7, 6, 6, 6, 6, 6, 0},
                             {6, 6, 6, 6, 6, 6, 0},
                             {6, 5, 6, 6, 6, 7, 10},
                             {9, 9, 9, 5, 5, 6, 6},
                             {0, 0, 0, 0, 0, 0, 0}};
};

#endif
