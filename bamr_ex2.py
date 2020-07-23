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

bp=bamr.bamr_py(b'tews_threep_ligo',3)
bp.settings(inc_baryon_mass=True,addl_quants=True,verbose=2,
            norm_max=False,crust_from_L=False,
            compute_cthick=True,apply_intsc=True,
            cached_intsc=True,prior_eta=True)
bp.add_data_alt("6304",
                "data/shb18/6304_H_nopl_syst_wilm.o2",
                "data/shb18/6304_He_nopl_syst_wilm.o2",
                "like",0.7,"rescaled")
bp.add_data_alt("6397",
                "data/shb18/6397_H_syst_wilm.o2",
                "data/shb18/6397_He_syst_wilm3.o2",
                "like",0.7,"rescaled")
bp.add_data_alt("M13",
                "data/shs18/M13_H_rs.o2",
                "data/shs18/M13_He_rs.o2",
                "like",0.7,"rescaled_0")
bp.add_data_alt("M28",
                "data/shb18/M28_H_syst_wilm.o2",
                "data/shb18/M28_He_syst_wilm.o2",
                "like",0.7,"rescaled")
bp.add_data_alt("M30",
                "data/shb18/M30_H_syst_wilm.o2",
                "data/shb18/M30_He_syst_wilm.o2",
                "like",0.7,"rescaled")
bp.add_data_alt("wCen",
                "data/shb18/wCen_H_syst_wilm.o2",
                "data/shb18/wCen_H_syst_wilm.o2",
                "like",0.7,"rescaled")
bp.add_data_alt("X7",
                "data/shb18/X7_H_syst_wilm.o2",
                "data/shb18/X7_He_syst_wilm.o2",
                "like",0.7,"rescaled")
bp.add_data_alt("1810b",
                "data/nks15/1810.o2",
                "data/nks15/1810.o2",
                "weights",0.7,"mcarlo")
bp.add_data_alt("1724b",
                "data/nks15/1724.o2",
                "data/nks15/1724.o2",
                "weights",0.7,"mcarlo")
bp.add_data_alt("1702",
                "data/nat17/1702_D_X_int.o2",
                "data/nat17/1702_D_X_int.o2",
                "avgs",0.7,"hist2_table")
bp.add_data_alt("0030",
                "data/nicer/0030_st_pst.o2",
                "data/nicer/0030_st_pst.o2",
                "prob",0.7,"table3d")
bp.bamr_init()
bp.compute_point([1.322748e+01,4.884697e-01,3.257073e+01,
                  4.455746e+01,4.381964e-01,3.203634e+00,
                  5.517590e+00,6.987111e+00,1.182554e+00,
                  1.197660e+00,2.462268e-01,5.957640e-01,
                  5.420301e-01,4.968797e-01,6.401870e-01,
                  5.838393e-01,7.508562e-01,8.855128e-01,
                  9.716234e-01,5.129025e-01,8.293522e-01,
                  8.446897e-01,4.968797e-01,-1.816410e+00,
                  -1.635335e+00,-1.841114e+00,-4.185872e-01,
                  -1.178018e+00,-1.149620e+00,-1.801794e-01,
                  -4.783507e-01,-8.689520e-01,-8.779179e-01,
                  -1.635335e+00])
bp.summarize_tables()
