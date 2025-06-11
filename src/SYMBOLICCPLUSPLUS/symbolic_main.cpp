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

#ifndef _SYMBOLIC_MAIN_CPP_
#define _SYMBOLIC_MAIN_CPP_

// This file was created by David Hicks (DX).
// This file includes all of the CPP files (i.e., function definitions) to 
// build the symbolic functionality and include it into a single object file
// Please preserve the order of including the CPP files 
// (the order mirrors the header/forward declaration dependencies)

#include "symbolic_main.h"

#ifndef _SYMBOLIC_CLONING_CPP_
#include "cloning.cpp"
#endif

#ifndef _SYMBOLIC_IDENTITY_CPP_
#include "identity.cpp"
#endif

#ifndef _SYMBOLIC_SYMBOLIC_CPP_
#include "symbolic.cpp"  // SymbolicInterface, Symbolic ...
#endif

#ifndef _SYMBOLIC_EQUATION_CPP_
#include "equation.cpp"  //   Equation : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_NUMBER_CPP_
#include "number.cpp"    //   Numeric  : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_MATRIX_CPP_ //DX20200825 - needs to be before product.cpp
#include "matrix.cpp"
#endif

#ifndef _SYMBOLIC_PRODUCT_CPP_
#include "product.cpp"   //   Product  : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_SUM_CPP_
#include "sum.cpp"       //   Sum      : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_SYMBOL_CPP_
#include "symbol.cpp"    //   Symbol   : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_FUNCTIONS_CPP_
#include "functions.cpp" //     Sin    : Symbol ...
#endif

#ifndef _SYMBOLIC_SYMMATRIX_CPP_
#include "symmatrix.cpp" //   SymbolicMatrix : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_SYMERROR_CPP_
#include "symerror.cpp"  //   SymbolicError  : CloningSymbolicInterface
#endif

#ifndef _SYMBOLIC_CONSTANTS_CPP_
#include "constants.h" 
#endif

#ifndef _SYMBOLIC_INTEGRATE_CPP_
#include "integrate.cpp"
#endif

#ifndef _SYMBOLIC_SOLVE_CPP_
#include "solve.cpp"
#endif

#ifndef _SYMBOLIC_ARRAY_CPP_
#include "array.cpp"
#endif
#ifndef _SYMBOLIC_DERIVE_CPP_
#include "derive.cpp"
#endif
#ifndef _SYMBOLIC_MATNORM_CPP_
#include "matnorm.cpp"
#endif
#ifndef _SYMBOLIC_MULTINOMIAL_CPP_
#include "multinomial.cpp"
#endif
#ifndef _SYMBOLIC_POLYNOMIAL_CPP_
#include "polynomial.cpp"
#endif
#ifndef _SYMBOLIC_QUATERN_CPP_
#include "quatern.cpp"
#endif
#ifndef _SYMBOLIC_RATIONAL_CPP_
#include "rational.cpp"
#endif
#ifndef _SYMBOLIC_VECNORM_CPP_
#include "vecnorm.cpp"
#endif
#ifndef _SYMBOLIC_VECTOR_CPP_
#include "vector.cpp"
#endif
#ifndef _SYMBOLIC_VERYLONG_CPP_
#include "verylong.cpp"
#endif
#ifndef _SYMBOLIC_CPLUSPLUS_CPP_
#include "symbolicc++.cpp"
#endif

#endif
