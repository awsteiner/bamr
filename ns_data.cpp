/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
#include "ns_data.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace bamr;

int ns_data::add_data(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<5) {
    std::cerr << "Not enough arguments given to 'add-data'." << std::endl;
    return o2scl::exc_efailed;
  }
      
  source_names.push_back(sv[1]);
  source_fnames.push_back(sv[2]);
  slice_names.push_back(sv[3]);
  init_mass_fracs.push_back(o2scl::stod(sv[4]));
  if (sv.size()==6) {
    table_names.push_back(sv[5]);
  } else {
    table_names.push_back("");
  }
      
  n_sources++;

  return 0;
}
    
