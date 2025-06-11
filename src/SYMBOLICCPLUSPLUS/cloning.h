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


// cloning.h

#ifndef SYMBOLIC_CPLUSPLUS_CLONING
#define SYMBOLIC_CPLUSPLUS_CLONING

#include <typeinfo>

//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

namespace symbolic{ //DX20200625
  class Cloning
  {
    private: int refcount;
             void (*free_p)(Cloning*);
             template <class T> static void free(Cloning*);

    public:  Cloning();
             Cloning(const Cloning&);
             virtual ~Cloning();

             virtual Cloning *clone() const = 0;
             template <class T> static Cloning *clone(const T&);
             static void reference(Cloning*);
             static void unreference(Cloning*);
  };

  class CloningPtr
  {
    protected: Cloning *value;

    public:    CloningPtr();
               CloningPtr(const Cloning&);
               CloningPtr(const CloningPtr&);
               ~CloningPtr();

               CloningPtr &operator=(const Cloning&);
               CloningPtr &operator=(const CloningPtr&);
  };

  template <class T>
    class CastPtr: public CloningPtr
  {
    public:  CastPtr();
             CastPtr(const Cloning&);
             CastPtr(const CloningPtr&);
             ~CastPtr();

             T *operator->() const;
             T &operator*() const;
  };
} //namespace symbolic //DX20200625


#endif
