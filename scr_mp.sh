#! /usr/bin/bash

./bamr -threads 1 -set prefix out/mp \
	-set max_iters 100000 -set file_update_time 1800 \
	-set verbose 1 -set mcmc_verbose 1 \
	-set min_max_mass 2.0 -set norm_max 0 \
	-set addl_quants 1 -set inc_baryon_mass 1 \
	-set crust_from_L 0 -set compute_cthick 1 \
	-add-data-alt 6304 data/shb18/6304_H_nopl_syst_wilm.o2 \
	data/shb18/6304_He_nopl_syst_wilm.o2 like 0.7 rescaled \
	-add-data-alt 6397 data/shb18/6397_H_syst_wilm.o2 \
	data/shb18/6397_He_syst_wilm3.o2 like 0.7 rescaled \
	-add-data-alt M13 data/shs18/M13_H_rs.o2 \
	data/shs18/M13_He_rs.o2 like 0.7 rescaled_0 \
	-add-data-alt M28 data/shb18/M28_H_syst_wilm.o2 \
	data/shb18/M28_He_syst_wilm.o2 like 0.7 rescaled \
	-add-data-alt M30 data/egz20/M30_echi_H.o2 \
	data/egz20/M30_echi_He.o2 like 0.7 rescaled \
	-add-data-alt wCen data/shb18/wCen_H_syst_wilm.o2 \
	data/shb18/wCen_H_syst_wilm.o2 like 0.7 rescaled \
	-add-data-alt X7 data/shb18/X7_H_syst_wilm.o2 \
	data/shb18/X7_He_syst_wilm.o2 like 0.7 rescaled \
	-add-data-alt 1810b data/nks15/1810.o2 \
	data/nks15/1810.o2 weights 0.7 mcarlo \
	-add-data-alt 1724b data/nks15/1724.o2 \
	data/nks15/1724.o2 weights 0.7 mcarlo \
	-add-data-alt 1702 data/nat17/1702_D_X_int.o2 \
	data/nat17/1702_D_X_int.o2 avgs 0.7 hist2_table \
	-add-data-alt 0030 data/nicer/0030_st_pst.o2 \
	data/nicer/0030_st_pst.o2 prob 0.7 table3d \
	-add-data-alt 0740 data/nicer/J0740_H_MR_t3d.o2 \
	data/nicer/J0740_H_MR_t3d.o2 prob 0.7 table3d \
	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
        -method nsf -model new_poly -set model_dpdm 1 \
	-initial-point-best out/mp_train \
	-mcmc
#> pt_nsf.out

# ./bamr -threads 1 -set prefix out/mp \
# 	-set max_iters 100000 -set file_update_time 1800 \
# 	-set verbose 3 -set mcmc_verbose 3 \
# 	-set min_max_mass 2.0 -set norm_max 0 \
# 	-set addl_quants 1 -set inc_baryon_mass 1 \
# 	-set crust_from_L 0 -set compute_cthick 1 \
# 	-add-data-alt 6304 data/shb18/6304_H_nopl_syst_wilm.o2 \
# 	data/shb18/6304_He_nopl_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt 6397 data/shb18/6397_H_syst_wilm.o2 \
# 	data/shb18/6397_He_syst_wilm3.o2 like 0.7 rescaled \
# 	-add-data-alt M13 data/shs18/M13_H_rs.o2 \
# 	data/shs18/M13_He_rs.o2 like 0.7 rescaled_0 \
# 	-add-data-alt M28 data/shb18/M28_H_syst_wilm.o2 \
# 	data/shb18/M28_He_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt M30 data/egz20/M30_echi_H.o2 \
# 	data/egz20/M30_echi_He.o2 like 0.7 rescaled \
# 	-add-data-alt wCen data/shb18/wCen_H_syst_wilm.o2 \
# 	data/shb18/wCen_H_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt X7 data/shb18/X7_H_syst_wilm.o2 \
# 	data/shb18/X7_He_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt 1810b data/nks15/1810.o2 \
# 	data/nks15/1810.o2 weights 0.7 mcarlo \
# 	-add-data-alt 1724b data/nks15/1724.o2 \
# 	data/nks15/1724.o2 weights 0.7 mcarlo \
# 	-add-data-alt 1702 data/nat17/1702_D_X_int.o2 \
# 	data/nat17/1702_D_X_int.o2 avgs 0.7 hist2_table \
# 	-add-data-alt 0030 data/nicer/0030_st_pst.o2 \
# 	data/nicer/0030_st_pst.o2 prob 0.7 table3d \
# 	-add-data-alt 0740 data/nicer/J0740_H_MR_t3d.o2 \
# 	data/nicer/J0740_H_MR_t3d.o2 prob 0.7 table3d \
# 	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
#         -method kde -model new_poly -set model_dpdm 1 \
# 	-initial-point-best out/mp_train \
# 	-mcmc > pt_kde.out

# ./bamr -threads 1 -set prefix out/mp \
# 	-set max_iters 100000 -set file_update_time 1800 \
# 	-set verbose 1 -set mcmc_verbose 1 \
# 	-set min_max_mass 2.0 -set norm_max 0 \
# 	-set addl_quants 1 -set inc_baryon_mass 1 \
# 	-set crust_from_L 0 -set compute_cthick 1 \
# 	-add-data-alt 6304 data/shb18/6304_H_nopl_syst_wilm.o2 \
# 	data/shb18/6304_He_nopl_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt 6397 data/shb18/6397_H_syst_wilm.o2 \
# 	data/shb18/6397_He_syst_wilm3.o2 like 0.7 rescaled \
# 	-add-data-alt M13 data/shs18/M13_H_rs.o2 \
# 	data/shs18/M13_He_rs.o2 like 0.7 rescaled_0 \
# 	-add-data-alt M28 data/shb18/M28_H_syst_wilm.o2 \
# 	data/shb18/M28_He_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt M30 data/egz20/M30_echi_H.o2 \
# 	data/egz20/M30_echi_He.o2 like 0.7 rescaled \
# 	-add-data-alt wCen data/shb18/wCen_H_syst_wilm.o2 \
# 	data/shb18/wCen_H_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt X7 data/shb18/X7_H_syst_wilm.o2 \
# 	data/shb18/X7_He_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt 1810b data/nks15/1810.o2 \
# 	data/nks15/1810.o2 weights 0.7 mcarlo \
# 	-add-data-alt 1724b data/nks15/1724.o2 \
# 	data/nks15/1724.o2 weights 0.7 mcarlo \
# 	-add-data-alt 1702 data/nat17/1702_D_X_int.o2 \
# 	data/nat17/1702_D_X_int.o2 avgs 0.7 hist2_table \
# 	-add-data-alt 0030 data/nicer/0030_st_pst.o2 \
# 	data/nicer/0030_st_pst.o2 prob 0.7 table3d \
# 	-add-data-alt 0740 data/nicer/J0740_H_MR_t3d.o2 \
# 	data/nicer/J0740_H_MR_t3d.o2 prob 0.7 table3d \
# 	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
#         -method kde_sklearn -model new_poly -set model_dpdm 1 \
# 	-initial-point-best out/mp_train \
# 	-mcmc > pt_kde_sklearn.out

# ./bamr -threads 1 -set prefix out/mp \
# 	-set max_iters 100000 -set file_update_time 1800 \
# 	-set verbose 1 -set mcmc_verbose 1 \
# 	-set min_max_mass 2.0 -set norm_max 0 \
# 	-set addl_quants 1 -set inc_baryon_mass 1 \
# 	-set crust_from_L 0 -set compute_cthick 1 \
# 	-add-data-alt 6304 data/shb18/6304_H_nopl_syst_wilm.o2 \
# 	data/shb18/6304_He_nopl_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt 6397 data/shb18/6397_H_syst_wilm.o2 \
# 	data/shb18/6397_He_syst_wilm3.o2 like 0.7 rescaled \
# 	-add-data-alt M13 data/shs18/M13_H_rs.o2 \
# 	data/shs18/M13_He_rs.o2 like 0.7 rescaled_0 \
# 	-add-data-alt M28 data/shb18/M28_H_syst_wilm.o2 \
# 	data/shb18/M28_He_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt M30 data/egz20/M30_echi_H.o2 \
# 	data/egz20/M30_echi_He.o2 like 0.7 rescaled \
# 	-add-data-alt wCen data/shb18/wCen_H_syst_wilm.o2 \
# 	data/shb18/wCen_H_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt X7 data/shb18/X7_H_syst_wilm.o2 \
# 	data/shb18/X7_He_syst_wilm.o2 like 0.7 rescaled \
# 	-add-data-alt 1810b data/nks15/1810.o2 \
# 	data/nks15/1810.o2 weights 0.7 mcarlo \
# 	-add-data-alt 1724b data/nks15/1724.o2 \
# 	data/nks15/1724.o2 weights 0.7 mcarlo \
# 	-add-data-alt 1702 data/nat17/1702_D_X_int.o2 \
# 	data/nat17/1702_D_X_int.o2 avgs 0.7 hist2_table \
# 	-add-data-alt 0030 data/nicer/0030_st_pst.o2 \
# 	data/nicer/0030_st_pst.o2 prob 0.7 table3d \
# 	-add-data-alt 0740 data/nicer/J0740_H_MR_t3d.o2 \
# 	data/nicer/J0740_H_MR_t3d.o2 prob 0.7 table3d \
# 	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
#         -method gauss -model new_poly -set model_dpdm 1 \
# 	-initial-point-best out/mp_train \
# 	-mcmc
#> pt_gauss.out

