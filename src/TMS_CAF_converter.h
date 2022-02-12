#ifndef __TMS_CAF_CONVERTER__
#define __TMS_CAF_CONVERTER__

// Can only run this if we have linkage to CAF format
#ifdef __DUNEANAOBJ_ENABLED__

#include <iostream>

// Include the general TMS_Event class
#include "TMS_Event.h"

// Include standard record format for TMS
#include "SRTMS.h"

// Converter for TMS format to CAF format
class TMS_CAF_converter {
  public:
    TMS_CAF_converter();
    ~TMS_CAF_converter() {};
    caf::SRTMS ConvertEvent(TMS_event &event);
  private:
}

#endif // End CAF linkage check
#endif // End header guard
