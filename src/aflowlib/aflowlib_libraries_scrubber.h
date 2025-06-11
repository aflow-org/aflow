
#ifndef AFLOWLIB_LIBRARIES_SCRUBBER_H
#define AFLOWLIB_LIBRARIES_SCRUBBER_H

#include <string>

#ifndef uint
typedef unsigned uint;
#endif

namespace aflowlib {
  uint LIB2SCRUB(std::string library, bool VERBOSE);
  bool LIB2AUID(std::string entry, bool TEST, bool _VERBOSE);
  uint MOSFET(int mode, bool VERBOSE);
  uint MAIL2SCAN(std::string library, bool VERBOSE);
} // namespace aflowlib

#endif // AFLOWLIB_LIBRARIES_SCRUBBER_H
