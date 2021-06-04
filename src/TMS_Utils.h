#ifndef __TMS_UTILS_H__
#define __TMS_UTILS_H__

namespace TMS_Utils {
  // Get the sign of a number
  template <typename T> inline int sgn(T value) {
    return (T(0)<value - (value<T(0)));
  }
}

#endif
