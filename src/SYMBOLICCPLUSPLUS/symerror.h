//   SymbolicC++ : An object oriented computer algebra system written in C++
//
//   Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb
//
//   This library is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This library is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


// symerror.h

#ifndef SYMBOLIC_CPLUSPLUS_ERRORS

#include <string>

//DX20200831 [OBSOLETE] #ifdef  SYMBOLIC_FORWARD
//DX20200831 [OBSOLETE] #ifndef SYMBOLIC_CPLUSPLUS_ERRORS_FORWARD
//DX20200831 [OBSOLETE] #define SYMBOLIC_CPLUSPLUS_ERRORS_FORWARD

namespace symbolic{ //DX20200625
  class SymbolicError;
} //namespace symbolic //DX20200625

//DX20200831 [OBSOLETE] #endif
//DX20200831 [OBSOLETE] #endif

#ifdef  SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_ERRORS_DECLARE
#define SYMBOLIC_CPLUSPLUS_ERRORS_DECLARE

namespace symbolic{ //DX20200625
  class SymbolicError
  {
    public: typedef enum {
              IncompatibleNumeric,
                IncompatibleVector,
                NoMatch,
                NotDouble,
                NotInt,
                NotMatrix,
                NotNumeric,
                NotVector,
                UnsupportedNumeric
            } error;

            error errornumber;
            SymbolicError(const error &e);
            string message() const;
  };
} //namespace symbolic //DX20200625

#endif
#endif

//DX20200824 - moved function definitions to symerror.cpp

#endif
