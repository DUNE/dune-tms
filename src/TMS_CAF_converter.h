#ifndef __TMS_CAF_CONVERTER__
#define __TMS_CAF_CONVERTER__

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED

#include <iostream>

// Include the general TMS_Event class
#include "TMS_Event.h"
// The reco class
#include "TMS_Reco.h"

// Include standard record format for TMS
#include "duneanaobj/StandardRecord/SRTMS.h"

// Converter for TMS format to CAF format
class TMS_CAF_converter {
  public:
    TMS_CAF_converter() {};
    ~TMS_CAF_converter() {};

    // The converter
    //caf::SRTMS ConvertEvent(TMS_Event &event);
    caf::SRTMS ConvertEvent();
  private:

};

#endif // End CAF linkage check
#endif // End header guard
