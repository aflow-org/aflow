

#ifndef AUROSTD_DEFS_H
#define AUROSTD_DEFS_H

#include <cstdio>
#include <limits>
#include <string>

#ifndef SWAP
#define SWAP(a, b) \
  {                \
    temp = (a);    \
    (a) = (b);     \
    (b) = temp;    \
  }
#endif

#ifndef NNN
#define NNN -123456
#endif

#ifndef AUROSTD_NAN
#define AUROSTD_NAN 1E9
#endif

#ifndef AUROSTD_DEFAULT_PRECISION
#define AUROSTD_DEFAULT_PRECISION 20
#endif

// CO20171215 - more global control
// this is for PRINTING, not for algorithmic zeroing
#ifndef AUROSTD_ROUNDOFF_TOL
#define AUROSTD_ROUNDOFF_TOL 1e-6
#endif

// CO20171215 - more global control
// this is SMALLER than _ZERO_TOL_ in aflow.h
//(aurostd default is less conservative to match read/write roundoff error)
#ifndef AUROSTD_IDENTITY_TOL
#define AUROSTD_IDENTITY_TOL 1e-6
#endif

// CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_INT
#define AUROSTD_MAX_INT std::numeric_limits<int>::max()
#endif

// CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_UINT
#define AUROSTD_MAX_UINT std::numeric_limits<uint>::max()
#endif

// CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_ULLINT
#define AUROSTD_MAX_ULLINT std::numeric_limits<unsigned long long int>::max()
#endif

// CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_DOUBLE
#define AUROSTD_MAX_DOUBLE std::numeric_limits<double>::max()
#endif

// SD20220914
#ifndef AUROSTD_MAX_SIZET
#define AUROSTD_MAX_SIZET std::numeric_limits<size_t>::max()
#endif

// CO20180101 - stream2stream modes
#ifndef DEFAULT_STREAM
#define DEFAULT_STREAM 'D'
#endif
#ifndef FIXED_STREAM
#define FIXED_STREAM 'F'
#endif
#ifndef SCIENTIFIC_STREAM
#define SCIENTIFIC_STREAM 'S'
#endif

#ifndef __STRICT_ANSI__
#define __STRICT_ANSI__
#endif

#define _AUROSTD_XLIBS_ERROR_ std::string("ERROR - AUROSTD_XLIBS++ [***]:")
#define _AUROSTD_XLIBS_WARNING_ std::string("WARNING - AUROSTD_XLIBS++ [***]:")

#define EMPTY_WORDING std::string("blablabla")

// CO20200624 START - adding from Jahnatek
// http://ascii-table.com/ansi-escape-sequences.php
// http://ascii-table.com/ansi-escape-sequences-vt-100.php
#define cursor_moveyx(y, x, fstr) std::fprintf(fstr, "\033[%d;%dH", y, x) // Move cursor to position y,x (rows, columns) with (1,1) as origin
#define cursor_moveup(y, fstr) std::fprintf(fstr, "\033[%dA", y)          // Move cursor up y
#define cursor_movedown(y, fstr) std::fprintf(fstr, "\033[%dB", y)        // Move cursor down y
#define cursor_moveright(x, fstr) std::fprintf(fstr, "\033[%dC", x)       // Move cursor right x
#define cursor_moveleft(x, fstr) std::fprintf(fstr, "\033[%dD", x)        // Move cursor left x
#define cursor_store(fstr) std::fprintf(fstr, "\033[s")                 // Store current cursor position and color
#define cursor_restore(fstr) std::fprintf(fstr, "\033[u")               // Restore cursor position and color from cursor_store()
#define cursor_clear(fstr) std::fprintf(fstr, "\033[2J")                // Clear screen and leave cursor where is
#define cursor_clearline(fstr) std::fprintf(fstr, "\033[K")             // Clear to end of line and leave cursor where is
#define cursor_fore_black(fstr) std::fprintf(fstr, "\033[30m")          // Change foreground color to black
#define cursor_fore_red(fstr) std::fprintf(fstr, "\033[31m")            // Change foreground color to red
#define cursor_fore_green(fstr) std::fprintf(fstr, "\033[32m")          // Change foreground color to green
#define cursor_fore_orange(fstr) std::fprintf(fstr, "\033[33m")         // Change foreground color to orange
#define cursor_fore_blue(fstr) std::fprintf(fstr, "\033[34m")           // Change foreground color to blue
#define cursor_fore_magenta(fstr) std::fprintf(fstr, "\033[35m")        // Change foreground color to magenta
#define cursor_fore_cyan(fstr) std::fprintf(fstr, "\033[36m")           // Change foreground color to cyan
#define cursor_fore_yellow(fstr) std::fprintf(fstr, "\033[33m\033[1m")  // Change foreground color to yellow (add bold to help visibility)
#define cursor_fore_white(fstr) std::fprintf(fstr, "\033[37m")          // Change foreground color to white
#define cursor_back_black(fstr) std::fprintf(fstr, "\033[40m")          // Change background color to black
#define cursor_back_red(fstr) std::fprintf(fstr, "\033[41m")            // Change background color to red
#define cursor_back_green(fstr) std::fprintf(fstr, "\033[42m")          // Change background color to green
#define cursor_back_orange(fstr) std::fprintf(fstr, "\033[43m")         // Change background color to orange
#define cursor_back_blue(fstr) std::fprintf(fstr, "\033[44m")           // Change background color to blue
#define cursor_back_magenta(fstr) std::fprintf(fstr, "\033[45m")        // Change background color to magenta
#define cursor_back_cyan(fstr) std::fprintf(fstr, "\033[46m")           // Change background color to cyan
#define cursor_back_white(fstr) std::fprintf(fstr, "\033[47m")          // Change background color to white
#define cursor_attr_none(fstr) std::fprintf(fstr, "\033[0m")            // Turn off all cursor attributes
#define cursor_attr_bold(fstr) std::fprintf(fstr, "\033[1m")            // Make test bold
#define cursor_attr_underline(fstr) std::fprintf(fstr, "\033[4m")       // Underline text
#define cursor_attr_blink(fstr) std::fprintf(fstr, "\033[5m")           // Supposed to make text blink, usually bolds it instead
#define cursor_attr_reverse(fstr) std::fprintf(fstr, "\033[7m")         // Swap background and foreground colors
// CO20200624 END - adding from Jahnatek

#define AUROSTD_ZIP_BIN "xz"
#define AUROSTD_ZIP_EXT ".xz"

#define __GNU_CPP
#ifndef __xprototype
#ifdef GNU
#define __xprototype __attribute__((const))
#else
#define __xprototype
#endif
#endif

#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)  // CO20200502 - moved from aflow.h

// CO20200502 START - including gettid()
#ifdef __GLIBC__
#define GLIBC_VERSION (__GLIBC__ * 100 + __GLIBC_MINOR__)
#if (GLIBC_VERSION < 230) // CO20200502 - apparently they patched at 230 - https://stackoverflow.com/questions/30680550/c-gettid-was-not-declared-in-this-scope
//[CO20200502 - too many warnings]#warning "defining getid() with syscall(SYS_gettid)"
#include <sys/syscall.h>  //CO20200502 - need for gettid()
#define gettid() syscall(SYS_gettid)
#endif  // GLIBC_VERSION
#endif  //__GLIBC__
// CO20200502 END - including gettid()

#ifndef __xprototype
#define __xprototype __attribute__((const))
#endif

// ----------------------------------------------------------------------------
// ---------------------------------------------------------- physics constants

#define rad2deg (180.0 / 3.14159265358979323846)
#define deg2rad (3.14159265358979323846 / 180.0)
#define angstrom2bohr (1 / 0.529177249)
#define bohr2angstrom (0.529177249)

// ME20181020 - check out https://physics.nist.gov/cuu/Constants/index.html
// #define PI                    3.14159265359
#define PI 3.14159265358979323846
#define pi PI
#define TWOPI 6.28318530717958647692
#define EULERSNUMBER 2.71828182845904523536
#define RTPI (sqrt(PI))
#define C_VACUUM 2.99792458E+8                   // m/s
#define EPS_VACUUM 8.854187817E-12                 // C/(N-m2)
#define MU_VACUUM (4.0 * PI * 1.0E-7)                 // N/A^2 (T^2m^3/J)
#define AMU2KILOGRAM 1.66054E-27
#define KILOGRAM2AMU 6.022137E+26
#define E_ELECTRON 1.60217662E-19                  // C
#define eV2J E_ELECTRON                      // 1eV=E_ELECTRON J //CO20201111
#define J2eV (1.0 / E_ELECTRON)                // CO20201111
#define PLANCKSCONSTANT_h 6.62607E-34                     // J*s
#define PLANCKSCONSTANT_hbar 1.0545718E-34                   // J*s
#define PLANCKSCONSTANTEV_h (PLANCKSCONSTANT_h / E_ELECTRON)  // eV*s
#define PLANCKSCONSTANTEV_hbar (PLANCKSCONSTANTEV_h / TWOPI)     // eV*s
#define KBOLTZ 1.3806504E-23                   // J/K
#define KBOLTZEV (KBOLTZ / E_ELECTRON)             // eV/K //8.617343E-5
#define eV2K (11604.505)                     // 1eV=11604.505 K
#define meV2K (11.604505)                     // 1meV=11.604505 K
#define mol2atom 6.0221408E23                    // 1mol=6.022e23 atoms    //CO20180329
#define atom2mol (1.0 / 6.0221408E23)              // CO20201111
#define eVatom2kJmol (E_ELECTRON * mol2atom / 1.0e3)     // 1eV/atom=96.5kJ/mol    //CO20180329
#define meVatom2kJmol (eVatom2kJmol / 1.0e3)            // 1meV/atom=0.0965kJ/mol //CO20180329
#define hartree2eV 27.2113862459                   // 1hartree=27.211eV      //ME20200206
#define eV2hartree (1.0 / hartree2eV)                // 1eV=0.0367hartree      //SD20230919
#define kcal2eV 4.336443203200000E-002          // 1(kcal/mol) = 4.33644E-2 eV
#define eV2kcal (1.0 / kcal2eV)

// ME20200107 - (A)APL conversion factors
#define THz2Hz 1E12
#define Hz2THz (1.0 / THz2Hz)
#define au2Hz aurostd::sqrt(E_ELECTRON * 1E20 / AMU2KILOGRAM)  // eV/(A^2 amu) -> Hz
#define au2rcm (au2Hz / (100 * C_VACUUM))                         // eV/(A^2 amu) -> cm^-1
#define au2eV (au2Hz * PLANCKSCONSTANTEV_h)                    // eV/(A^2 amu) -> eV
#define eV2Hz (1.0 / PLANCKSCONSTANTEV_h)
#define eV2rcm (1.0 / (PLANCKSCONSTANTEV_h * 100 * C_VACUUM))
#define au2nmTHz ((E_ELECTRON * Hz2THz * Hz2THz * 1E18) / (0.1 * AMU2KILOGRAM))  // eV/(A amu) -> nm * THz^2
#define PLANCKSCONSTANT_h_THz (PLANCKSCONSTANT_h * THz2Hz) // J/THz
#define PLANCKSCONSTANT_hbar_THz (PLANCKSCONSTANT_hbar * THz2Hz) // J/THz
#define PLANCKSCONSTANTAMU_hbar_THz (PLANCKSCONSTANTEV_hbar * THz2Hz * (10 * au2nmTHz))  // amu A^2 THz
#define BEfactor_hbar_THz (PLANCKSCONSTANTEV_hbar / (KBOLTZEV * Hz2THz))  // hbar/kB in K/THz
#define BEfactor_h_THz (PLANCKSCONSTANTEV_h / (KBOLTZEV * Hz2THz))  // h/kB in K/THz

// AS20200427 - QHA-related conversion factors
#define eV2GPa (E_ELECTRON * 1e21)    // [eV/A^3] --> [GPa]
#define GPa2eV (1.0 / eV2GPa)         // [GPa] --> [eV/A^3]
#define eV2kBar (eV2GPa * 10)         // [eV/A^3] --> [kBar]
#define kBar2eV (1.0 / eV2kBar)       // [kBar] --> [eV/A^3]
#define kBar2GPa 0.1                 // [kBar] --> [GPa]
#define atm2Pa 101325

// DX20210111 - GFA factors
#define TEMPERATURE_ROOM 300.0               // K
#define kBT_ROOM (KBOLTZEV * TEMPERATURE_ROOM) // 0.025

#define LIBRARY_NOTHING 256

#define LIBRARY_ALL 100

#define _LOCK_LINK_SUFFIX_ std::string(".init")

// ZERO PRECISION DEFINITIONS - TIGHT (DEFAULT) AND LOOSE
#define _AUROSTD_ZERO_PRECISION_ 10
#define _AUROSTD_ZERO_TOL_ std::pow(10, -_AUROSTD_ZERO_PRECISION_) // DX
#define _AUROSTD_ZERO_PRECISION_LOOSE_ 3
#define _AUROSTD_ZERO_TOL_LOOSE_ std::pow(10, -_AUROSTD_ZERO_PRECISION_LOOSE_) // DX
#define _AUROSTD_DOUBLE_PRECISION_ 8
#define _AUROSTD_DOUBLE_TOL_ std::pow(10, -_AUROSTD_DOUBLE_PRECISION_)
#define _AUROSTD_FLOAT_PRECISION_ 6
#define _AUROSTD_FLOAT_TOL_ std::pow(10, -_AUROSTD_FLOAT_PRECISION_)  // ME20200519 - tolerance for float precision

#endif // AUROSTD_DEFS_H
