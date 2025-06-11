// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// Corey Oses
// Frisco Rose (AFLUX)
// Marco Esters (AFLOW DB)
// Hagen Eckert (EntryLoader 2022)

#ifndef _AFLOWLIB_H_
#define _AFLOWLIB_H_

#include <string>

#include "aflowlib/aflowlib_database.h"
#include "aflowlib/aflowlib_entry_loader.h"
#include "aflowlib/aflowlib_libraries.h"
#include "aflowlib/aflowlib_libraries_scrubber.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "aflowlib/aflowlib_web_outreach.h"

#ifndef NOSG
#define NOSG std::string("NNN #0")
#endif
#define CHMODWEB false
#define INF 1E9
#define ENERGY_ATOM_ERROR_meV 50
#define PRESSURE_ZERO_ENTHALPY_ENERGY _FLOAT_TOL_ // 1e-6

// ***************************************************************************
#define HTRESOURCE_MODE_NONE 0
#define HTRESOURCE_MODE_PHP_AUTHOR 4
#define HTRESOURCE_MODE_PHP_THRUST 5
#define HTRESOURCE_MODE_PHP_ALLOY 6

#endif //  _AFLOWLIB_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
