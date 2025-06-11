
#ifndef AFLOW_WYCKOFF_H
#define AFLOW_WYCKOFF_H

#include <ostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xvector.h"

#include "flow/aflow_xclasses.h"
#include "structure/aflow_xatom.h"

aurostd::xvector<double> wv(const double& x, const double& y, const double& z);
void wa(_atom& a, xstructure& str);
xstructure WyckoffPOSITIONS(uint spacegroup, xstructure strin);
xstructure WyckoffPOSITIONS(uint spacegroup, uint option, xstructure strin);

// DX20181010 - grouped Wyckoff class - START
//  --------------------------------------------------------------------------
//  ===== GroupedWyckoffPosition Class ===== //
class GroupedWyckoffPosition {
public:
  GroupedWyckoffPosition();
  ~GroupedWyckoffPosition();
  friend std::ostream& operator<<(std::ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition);
  const GroupedWyckoffPosition& operator=(const GroupedWyckoffPosition& b);
  bool operator<(const GroupedWyckoffPosition& b) const;
  GroupedWyckoffPosition(const GroupedWyckoffPosition& b);
  uint type;
  std::string element;
  std::vector<std::string> site_symmetries;
  std::vector<uint> multiplicities;
  std::vector<std::string> letters;

private:
  void free();
  void copy(const GroupedWyckoffPosition& b);
};
// --------------------------------------------------------------------------
// DX20181010 - grouped Wyckoff class - END

#endif // AFLOW_WYCKOFF_H
