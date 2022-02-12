#ifndef __TMS_CAF_CONVERTER__
#define __TMS_CAF_CONVERTER__

// Can only run this if we have linkage to CAF format
#ifdef __CAF_LINKAGE__

#include <iostream>

// Include CAF
#include

// Converter for TMS format to CAF format
class TMS_CAF_converter {
  public:
    TMS_CAF_converter();
    ConvertEvent();
  private:
}

#endif // End CAF linkage check
#endif // End header guard
