#! /usr/bin/bash
./bamr_nompi -threads 1 -set aff_inv 0 -set couple_threads 0 \
		-set use_population 1 \
		-set min_max_mass 1.98 \
		-set prefix compare_threep_2 \
		-set max_time 86400 \
		-set n_walk 1 -set step_fac 4000.0 \
		-set norm_max 0 -set addl_quants 1 \
		-set inc_baryon_mass 1 \
		-set crust_from_L 0 -set compute_cthick 1 \
		-set file_update_time 10 -set verbose 2 \
		-set mcmc_verbose 2 -add-data-alt 6304 \
		data/shb18/6304_H_nopl_syst_wilm.o2 \
		data/shb18/6304_He_nopl_syst_wilm.o2 \
		like 0.7 rescaled \
		-add-data-alt 6397 \
		data/shb18/6397_H_syst_wilm.o2 \
		data/shb18/6397_He_syst_wilm3.o2 \
		like 0.7 rescaled \
		-add-data-alt M13 \
		data/shs18/M13_H_rs.o2 \
		data/shs18/M13_He_rs.o2 \
		like 0.7 rescaled_0 \
		-add-data-alt M28 \
		data/shb18/M28_H_syst_wilm.o2 \
		data/shb18/M28_He_syst_wilm.o2 \
		like 0.7 rescaled \
		-add-data-alt M30 \
		data/shb18/M30_H_syst_wilm.o2 \
		data/shb18/M30_He_syst_wilm.o2 \
		like 0.7 rescaled \
		-add-data-alt wCen \
		data/shb18/wCen_H_syst_wilm.o2 \
		data/shb18/wCen_H_syst_wilm.o2 \
		like 0.7 rescaled \
		-add-data-alt X7 \
		data/shb18/X7_H_syst_wilm.o2 \
		data/shb18/X7_He_syst_wilm.o2 \
		like 0.7 rescaled \
		-add-data-alt 1810b \
		data/nks15/1810.o2 \
		data/nks15/1810.o2 \
		weights 0.7 mcarlo \
		-add-data-alt 1724b \
		data/nks15/1724.o2 \
		data/nks15/1724.o2 \
		weights 0.7 mcarlo \
		-add-data-alt 1702 \
		data/nat17/1702_D_X_int.o2 \
		data/nat17/1702_D_X_int.o2 \
		avgs 0.7 hist2_table \
		-add-data-alt 0030 \
		data/nicer/0030_st_pst.o2 \
		data/nicer/0030_st_pst.o2 \
		prob 0.7 table3d \
		-set apply_intsc 0 \
		-set cached_intsc 0 \
		-model tews_threep_ligo \
		-initial-point-last guess \
		-set prior_eta 1 \
		-set mmax_deriv 1 \
		-mcmc 