/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2016, Andrew W. Steiner
  
  This file is part of Bamr.
  
  Bamr is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Bamr is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Bamr. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include "entry.h"

using namespace std;
using namespace o2scl;
using namespace bamr;

entry::entry() {
  np=0;
  ns=0;
}

entry::entry(size_t nps, size_t nso) {
  np=nps;
  params.resize(np);
  ns=nso;
  mass.resize(ns);
  rad.resize(ns);
}

int entry::allocate(size_t nps, size_t nso) {
  np=nps;
  params.resize(np);
  ns=nso;
  mass.resize(ns);
  rad.resize(ns);
  return 0;
}

entry &entry::operator=(const entry &e) {

  if (&e!=this) {
    np=e.np;
    params.resize(np);
    ns=e.ns;
    mass.resize(ns);
    rad.resize(ns);
      
    for(size_t i=0;i<np;i++) {
      params[i]=e.params[i];
    }
    for(size_t i=0;i<ns;i++) {
      mass[i]=e.mass[i];
      rad[i]=e.rad[i];
    }
  }

  return *this;
}
  
entry::entry(const entry &e) {

  if (&e!=this) {
    np=e.np;
    params.resize(np);
    ns=e.ns;
    mass.resize(ns);
    rad.resize(ns);
      
    for(size_t i=0;i<np;i++) {
      params[i]=e.params[i];
    }
    for(size_t i=0;i<ns;i++) {
      mass[i]=e.mass[i];
      rad[i]=e.rad[i];
    }
  }
}

std::ostream &bamr::operator<<(std::ostream &os, entry &e) {

  os << "EOS: ";
  for(size_t k=0;k<e.np-1;k++) {
    os << e.params[k] << " ";
  }
  os << e.params[e.np-1];
  if (e.ns>0) {
    os << " Mass: ";
    for(size_t k=0;k<e.ns;k++) {
      os << e.mass[k] << " ";
    }
    os << "Rad: ";
    for(size_t k=0;k<e.ns-1;k++) {
      os << e.rad[k] << " ";
    }
    os << e.rad[e.ns-1];
  }
  return os;
}

