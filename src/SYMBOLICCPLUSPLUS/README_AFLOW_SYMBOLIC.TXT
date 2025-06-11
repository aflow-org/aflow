****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
SymbolicC++ README for implentation into AFLOW
Describes the changes made to the SYMBOLICCPLUSPLUS/ subdirectory

README written by: David Hicks (DX; email: david.hicks@duke.edu)
****************************************************************************************************

Below is the licensing information for this source code:

    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

****************************************************************************************************

All changes from the original source are comments with "//DX<date>".

Modifications from original source code:

- The directory name was changed from SYMBOLICC++/ to SYMBOLICCPLUSPLUS/ since the former changes the Makefile formatting (i.e., colors for targets, etc.).

- Moved function definitions from the header files into CPP (.cpp) files.
  The original source code was all header (.h) files.
  This prevents inclusion of the code into more than one CPP file (i.e., violates the "One Definition Rule").

- The original source code had a class hierarchy (i.e., symbolic/ subdirectory) and a "redirect" header file.
  For clarity (and to ensure no circular/infinite dependency issues), the redirect was removed and all functions
  were moved to a single directory hierarchy.

- Wrapped all functions in this directory in namespace "symbolic". Otherwise, this can cause namespace conflicts.

- Removed importing of the standard ("std") namespace.  Otherwise, this can cause namespace conflicts.

- Removed directory hierarchy in the SYMBOLICCPLUSPLUS/ subdirectory; otherwise this does not work with the current Makefile scheme.
  Specifically, the original hierarchy had a header redirect which cause issues with the automatic dependency tree in aflow_makefile.cpp.

- Initialization of variables to fix compiler errors/warnings.
