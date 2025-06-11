
#ifndef AFLOW_SPACEGROUP_H
#define AFLOW_SPACEGROUP_H

#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"

#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

extern std::string LibrarySPACEGROUP;
namespace spacegroup {
  using std::string;
  class _spacegroup {
  public:
        // constructors/destructors
    _spacegroup();    // do nothing
    _spacegroup(const _spacegroup& b);    // do nothing
    ~_spacegroup();        // do nothing
        // OPERATORS                                                  // --------------------------------------
    const _spacegroup& operator=(const _spacegroup& b);             // some operators
        // CONTENT
    uint number;
    uint option;
    string stroption;
    string name;
    string sginfo;
    std::vector<_sym_op> fgroup;                                  // rotations/inversions + incell_translations operations
    std::vector<_sym_op> pgroup;                                  // rotations/inversions
  private:                                                       // ---------------------------------------
    void free();                                                  // to free everything
    void copy(const _spacegroup& b);                               // the flag is necessary because sometimes you need to allocate the space.
  };

  extern std::vector<_spacegroup> vspacegroups;
  uint SpaceGroupInitialize();
  bool SpaceGroupNumberStructure(xstructure& str);
  bool SpaceGroupOptionRequired(uint spacegroup);
} // namespace spacegroup

#endif // AFLOW_SPACEGROUP_H
