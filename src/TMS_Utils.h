#ifndef __TMS_UTILS_H__
#define __TMS_UTILS_H__

#include <iostream>

// Include the general TMS_Event class
#include "TMS_Event.h"
// The reco class
#include "TMS_Reco.h"

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED
// Include standard record format for TMS
#include "duneanaobj/StandardRecord/SRTMS.h"
#endif

namespace TMS_Utils {

  // Get the sign of a number
  template <typename T> inline int sgn(T value) {
    return (T(0)<value - (value<T(0)));
  }

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED
  caf::SRTMS ConvertEvent();
#endif

}

#endif
