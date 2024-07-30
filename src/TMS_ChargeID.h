#ifndef _TMS_CHARGEID_H_SEEN_
#define _TMS_CHARGEID_H_SEEN_

#include "TMS_Hit.h"
#include "TMS_Constants.h"

//Charge identification class
class TMS_ChargeID {
  public:
    TMS_ChargeID();
    // Functions to determine in which magnetic field orientation the track is
    bool region1(const TMS_Hit &Hit);
    bool region2(const TMS_Hit &Hit);
    bool region3(const TMS_Hit &Hit);

    // Function for the actual charge identification
    int ID_Track_Charge(const std::vector<TMS_Hit> &Track);

  private:
    TMS_Hit Hit;
    std::vector<TMS_Hit> Track;
};

#endif
