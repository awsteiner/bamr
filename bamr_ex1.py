#  -------------------------------------------------------------------
#  
#  Copyright (C) 2020, Andrew W. Steiner and Sarah Wellence
#  
#  This file is part of bamr.
#  
#  bamr is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  bamr is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with bamr. If not, see <http://www.gnu.org/licenses/>.
#  
#  -------------------------------------------------------------------

# An example of the use of the python bamr interface

import bamr
import time

bp=bamr.bamr_py(b'twop',1,True)
bp.settings(inc_baryon_mass=True,addl_quants=True,verbose=3)
bp.bamr_init()

start = time.time()
lw=bp.compute_point([1.0,-3.0,0.165,0.644,
                     1.51,0.576,4.6,1.21])
end = time.time()
print('log weight',lw,'time for one point',end-start)

bp.summarize_tables()
