# ----------------------------------------------------------------------
# This section has variables which may need to be modified. The rest
# of the makefile should not need too much editing.
# ----------------------------------------------------------------------

# This variable should include the directories for the O2scl, GSL, and
# HDF5 libraries. By default, this is taken from the enviroment
# variable LDFLAGS.

LDFLAGS = -L/usr/local/lib \
		-L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib

LIB_DIRS = $(LDFLAGS)

# This variable may need to be modified to specify the include
# directories for the GSL, Boost, HDF5, and O2scl header files. By
# default this is taken from the environment variable CXXFLAGS.

CXXFLAGS = -I/usr/lib/x86_64-linux-gnu/openmpi/include \
		-I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
		-I/usr/include/python3.8

INC_DIRS = $(CXXFLAGS)

PYTHON_LDFLAGS = $(SHELL_PYTHON_LDFLAGS) -lpython3.10
PYTHON_CXXFLAGS = $(SHELL_PYTHON_INCLUDES)

# Generic (no MPI necessary) C++ compiler (e.g. g++)
CXX = g++

# C++ compiler (e.g. mpicxx).
MPI_CXX = mpic++

# Set these to be empty if you do not have GNU readline support
READLINE_VAR = -DBAMR_READLINE
READLINE_LIB = -lreadline -lncurses

# Set these to be empty if you do not have FFTW3
FFTW_VAR = -DBAMR_FFTW3
FFTW_LIB = -lfftw3

# Basic compiler flags with and without MPI
COMPILER_FLAGS = -std=c++0x -O3 -Wall -Wno-unused
COMPILER_FLAGS_MPI = -std=c++0x -O3 -Wall -Wno-unused

# ----------------------------------------------------------------------
# UTK makefile
# ----------------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration
LIB_DIRS = $(UTKNA_O2SCL_LIBS)
INC_DIRS = $(UTKNA_O2SCL_INCS)
CXX = $(UTKNA_CXX)
MPI_CXX = $(UTKNA_MPI_CXX)
BAMR_DIR = $(UTKNA_BAMR_DIR)
COMPILER_FLAGS = $(UTKNA_CFLAGS)
COMPILER_FLAGS_MPI = $(UTKNA_MPI_CFLAGS)
COMPILER_FLAGS_OPENMP = $(UTKNA_OPENMP_FLAGS)
PYTHON_INCLUDES = $(UTKNA_PYTHON_INCLUDES)
PYTHON_LDFLAGS = $(UTKNA_PYTHON_LDFLAGS)

endif

# ----------------------------------------------------------------------
# Secondary variables
# ----------------------------------------------------------------------

ALL_FLAGS_MPI = $(COMPILER_FLAGS_MPI) $(INC_DIRS) $(READLINE_VAR) \
	$(FFTW_VAR) -DBAMR_MPI $(COMPILER_FLAGS_OPENMP) $(PYTHON_INCLUDES)

ALL_FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) $(FFTW_VAR) \
	$(PYTHON_INCLUDES)

LIBS = $(LIB_DIRS) $(PYTHON_LDFLAGS) $(READLINE_LIB) $(FFTW_LIB)

# ----------------------------------------------------------------------
# Targets for bamr
# ----------------------------------------------------------------------

bamr: bamr_class.o models.o nstar_cold2.o main.o mcmc_bamr.o ns_data.o \
	pop_data.o ns_pop.o
	$(MPI_CXX) $(ALL_FLAGS_MPI) $(LIB_DIRS) -o bamr \
	bamr_class.o models.o nstar_cold2.o main.o mcmc_bamr.o ns_data.o \
	pop_data.o ns_pop.o $(LIBS)

main.o: main.cpp
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o main.o -c main.cpp

nstar_cold2.o: nstar_cold2.cpp nstar_cold2.h
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o nstar_cold2.o -c nstar_cold2.cpp

ns_data.o: ns_data.cpp ns_data.h
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o ns_data.o -c ns_data.cpp

models.o: models.cpp models.h
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o models.o -c models.cpp

bamr_class.o: bamr_class.cpp bamr_class.h models.o main.o nstar_cold2.o
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o bamr_class.o -c bamr_class.cpp

#bc_wrap.o: bc_wrap.cpp bc_wrap.h models.o main.o nstar_cold2.o
#	$(MPI_CXX) $(ALL_FLAGS_MPI) -o bc_wrap.o -c bc_wrap.cpp

#emulator_bamr.o: emulator_bamr.cpp emulator_bamr.h
#	$(MPI_CXX) $(ALL_FLAGS_MPI) -o emulator_bamr.o -c emulator_bamr.cpp

mcmc_bamr.o: mcmc_bamr.cpp mcmc_bamr.h models.o main.o nstar_cold2.o
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o mcmc_bamr.o -c mcmc_bamr.cpp

ns_pop.o: ns_pop.cpp ns_pop.h
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o ns_pop.o -c ns_pop.cpp

pop_data.o: pop_data.cpp ns_pop.h
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o pop_data.o -c pop_data.cpp

# ----------------------------------------------------------------------
# Help target
# ----------------------------------------------------------------------

help:
	@echo "Main targets:"
	@echo "  bamr"
	@echo "  process"
	@echo "  bamr_nompi"
	@echo "  test_all: runs test_prep, and all MPI tests"
	@echo "  test_all_nompi: runs test_prep, and all _nompi tests"
	@echo "Developer targets:"
	@echo "  doc"
	@echo "  update-tags"

# ----------------------------------------------------------------------
# Targets for bamr_nompi
# ----------------------------------------------------------------------

bamr_nompi: bamr_class_nompi.o models_nompi.o nstar_cold2_nompi.o \
	main_nompi.o mcmc_bamr_nompi.o ns_data_nompi.o pop_data_nompi.o \
	ns_pop_nompi.o
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o bamr_nompi \
	bamr_class_nompi.o models_nompi.o nstar_cold2_nompi.o main_nompi.o \
	mcmc_bamr_nompi.o ns_data_nompi.o pop_data_nompi.o ns_pop_nompi.o \
	$(LIBS) 

main_nompi.o: main.cpp
	$(CXX) $(ALL_FLAGS) -o main_nompi.o -c main.cpp

nstar_cold2_nompi.o: nstar_cold2.cpp nstar_cold2.h
	$(CXX) $(ALL_FLAGS) -o nstar_cold2_nompi.o -c nstar_cold2.cpp

models_nompi.o: models.cpp models.h
	$(CXX) $(ALL_FLAGS) -o models_nompi.o -c models.cpp

ns_data_nompi.o: ns_data.cpp ns_data.h
	$(CXX) $(ALL_FLAGS) -o ns_data_nompi.o -c ns_data.cpp

#emulator_bamr_nompi.o: emulator_bamr.cpp emulator_bamr.h
#	$(CXX) $(ALL_FLAGS) -o emulator_bamr_nompi.o -c emulator_bamr.cpp 

mcmc_bamr_nompi.o: mcmc_bamr.cpp mcmc_bamr.h
	$(CXX) $(ALL_FLAGS) -o mcmc_bamr_nompi.o -c mcmc_bamr.cpp

bamr_class_nompi.o: bamr_class.cpp bamr_class.h
	$(CXX) $(ALL_FLAGS) -o bamr_class_nompi.o -c bamr_class.cpp

ns_pop_nompi.o: ns_pop.cpp ns_pop.h
	$(CXX) $(ALL_FLAGS) -o ns_pop_nompi.o -c ns_pop.cpp

pop_data_nompi.o: pop_data.cpp ns_pop.h
	$(CXX) $(ALL_FLAGS) -o pop_data_nompi.o -c pop_data.cpp

# ----------------------------------------------------------------------
# Targets for process
# ----------------------------------------------------------------------

process.o: process.cpp process.h
	$(CXX) $(ALL_FLAGS) -o process.o -c process.cpp

process_main.o: process_main.cpp
	$(CXX) $(ALL_FLAGS) -o process_main.o -c process_main.cpp

process: process.o process_main.o
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o process process_main.o \
		process.o $(LIBS) 

# ----------------------------------------------------------------------
# Targets for test
# ----------------------------------------------------------------------

test.o: test.cpp
	$(CXX) $(ALL_FLAGS) -o test.o -c test.cpp

test: test.o
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o test \
		test.o $(LIBS) 

# ----------------------------------------------------------------------
# Testing targets
# ----------------------------------------------------------------------

test_all: test_prep test_data test_nodata test_cthick test_fixp \
		test_qt test_qf test_warmup test_ai \
		test_addl test_crustL test_tableseq test_storej bamr \
		test_debug_eos test_debug_star test_rejtab

# test11 test12 

test_all_nompi: test_prep_nompi test_debug_eos_nompi test_debug_star_nompi \
		test_data_nompi test_nodata_nompi test_cthick_nompi \
		test_fixp_nompi test_warmup_nompi test_ai_nompi \
		test_addl_nompi test_crustL_nompi bamr_nompi

# We can't do qt or qf tests until we finish baryon density updates
#
#test_qt_nompi test_qf_nompi 
#test11_nompi test12_nompi 

test_prep: bamr
	-mkdir -p data_temp
	-rm -rf data_temp/*

test_prep_nompi: bamr_nompi
	-mkdir -p data_temp
	-rm -rf data_temp/*

# ----------------------------------------------------------------------
# Individual testing targets with MPI
# ----------------------------------------------------------------------

test_debug_eos: bamr
	mpirun -np 1 bamr -threads 1 -set prefix data_temp/debug_eos \
		-set debug_eos 1 -run default.in -model twop -mcmc \
		> data_temp/debug_eos.scr
	-mv -i debug_eos.o2 data_temp

test_debug_star: bamr
	mpirun -np 1 bamr -threads 1 -set prefix data_temp/debug_star \
		-set debug_star 1 -run default.in -model twop -mcmc \
		> data_temp/debug_star.scr
	-mv -i debug_star.o2 data_temp

test_data: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 300 -set prefix data_temp/twop_data \
		-run default.in -model twop -mcmc \
		> data_temp/twop_data.scr 2> data_temp/twop_data.err

test_data_rp: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 300 -set prefix data_temp/twop_data_rp \
		-run default.in -model twop \
		-read_prev_results "data_temp/twop_data_<rank>_out" -mcmc \
		> data_temp/twop_data_rp.scr 2> data_temp/twop_data_rp.err

test_nodata: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set prefix data_temp/twop_nodata \
		-model twop -mcmc \
		> data_temp/twop_nodata.scr 2> data_temp/twop_nodata.err

test_cthick: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set compute_cthick 1 \
		-set prefix data_temp/twop_cthick -model twop -mcmc \
		> data_temp/twop_cthick.scr 2> data_temp/twop_cthick.err

test_fixp: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set prefix data_temp/fixp_nodata \
		-model fixp -mcmc \
		> data_temp/fixp_nodata.scr 2> data_temp/fixp_nodata.err

test_qt: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 300 -set prefix data_temp/qt_nodata \
		-model qmc_threep -mcmc \
		> data_temp/qt_nodata.scr 2> data_temp/qt_nodata.err

test_qf: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set prefix data_temp/qf_nodata \
		-model qmc_fixp -mcmc \
		> data_temp/qf_nodata.scr 2> data_temp/qf_nodata.err

test_warmup: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set n_warm_up 100 \
		-set prefix data_temp/twop_warmup -model twop -mcmc \
		> data_temp/twop_warmup.scr 2> data_temp/twop_warmup.err

test_ai: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 \
		-set prefix data_temp/twop_ai -set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 10 -mcmc \
		> data_temp/twop_ai.scr 2> data_temp/twop_ai.err

# These are unnecessary now because there is no max_chain_size setting

# test11: bamr
# 	mpirun -np 2 \
# 	bamr -threads 2 -set max_iters 100 -set prefix data_temp/twop_chain \
# 		-set max_chain_size 10 -model twop -mcmc \
# 		> data_temp/twop_chain.scr 2> data_temp/twop_chain.err

# test12: bamr
# 	mpirun -np 2 \
# 	bamr -threads 2 -set max_iters 100 -set prefix data_temp/twop_aic \
# 		-set aff_inv 1 \
# 		-set step_fac 2.0 -model twop -set n_walk 20 \
# 		-set max_chain_size 10 -mcmc \
# 		> data_temp/twop_aic.scr 2> data_temp/twop_aic.err

test_addl: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set addl_quants 1 -set inc_baryon_mass 1 \
		-set prefix data_temp/twop_addl -model twop -mcmc \
		> data_temp/twop_addl.scr 2> data_temp/twop_addl.err

test_crustL: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 100 -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set prefix data_temp/twop_crustL -model twop -mcmc \
		> data_temp/twop_crustL.scr 2> data_temp/twop_crustL.err

test_storej: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 300 \
		-set prefix data_temp/twop_storej -set store_rejects 1 \
		-run default.in -model twop -mcmc \
		> data_temp/twop_storej.scr 2> data_temp/twop_storej.err

test_tableseq: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 300 \
		-set prefix data_temp/twop_tableseq -set table_sequence 0 \
		-run default.in -model twop -mcmc \
		> data_temp/twop_tableseq.scr 2> data_temp/twop_tableseq.err

test_rejtab: bamr
	mpirun -np 2 \
	bamr -threads 2 -set max_iters 300 -set prefix data_temp/twop_rejtab \
		-set table_sequence 0 -set store_rejects 1 \
		-run default.in -model twop -mcmc \
		> data_temp/twop_rejtab.scr 2> data_temp/twop_rejtab.err

# ----------------------------------------------------------------------
# Individual testing targets without MPI
# ----------------------------------------------------------------------

NTEST = 10

test_debug_eos_nompi: bamr_nompi
	./bamr_nompi -set prefix data_temp/debug_eos \
		-set debug_eos 1 -run default.in -model twop -mcmc \
		> data_temp/debug_eos.scr
	-mv -i debug_eos.o2 data_temp

test_debug_star_nompi: bamr_nompi
	./bamr_nompi -set prefix data_temp/debug_star \
		-set debug_star 1 -run default.in -model twop -mcmc \
		> data_temp/debug_star.scr
	-mv -i debug_star.o2 data_temp

test_data_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/twop_data \
		-run default.in -model twop -mcmc \
		> data_temp/twop_data.scr 2> data_temp/twop_data.err

test_nodata_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/twop_nodata \
		-model twop -mcmc \
		> data_temp/twop_nodata.scr 2> data_temp/twop_nodata.err

test_cthick_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set compute_cthick 1 \
		-set prefix data_temp/twop_cthick -model twop -mcmc \
		> data_temp/twop_cthick.scr 2> data_temp/twop_cthick.err

test_fixp_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/fixp_nodata \
		-model fixp -mcmc \
		> data_temp/fixp_nodata.scr 2> data_temp/fixp_nodata.err

test_qt_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/qt_nodata \
		-model qmc_threep -mcmc \
		> data_temp/qt_nodata.scr 2> data_temp/qt_nodata.err

test_qf_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/qf_nodata \
		-model qmc_fixp -mcmc \
		> data_temp/qf_nodata.scr 2> data_temp/qf_nodata.err

test_warmup_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set n_warm_up $(NTEST) \
		-set prefix data_temp/twop_warmup -model twop -mcmc \
		> data_temp/twop_warmup.scr 2> data_temp/twop_warmup.err

test_ai_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/twop_ai \
		-set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 10 -mcmc \
		> data_temp/twop_ai.scr 2> data_temp/twop_ai.err

# These are unnecessary now because there is no max_chain_size setting

# test11_nompi: bamr_nompi
# 	./bamr_nompi -set max_iters 100 -set prefix data_temp/twop_chain \
# 		-set max_chain_size 10 -model twop -mcmc \
# 		> data_temp/twop_chain.scr 2> data_temp/twop_chain.err

# test12_nompi: bamr_nompi
# 	./bamr_nompi -set max_iters $(NTEST) -set prefix data_temp/twop_aic \
# 		-set aff_inv 1 \
# 		-set step_fac 2.0 -model twop -set n_walk 20 \
# 		-set max_chain_size 10 -mcmc \
# 		> data_temp/twop_aic.scr 2> data_temp/twop_aic.err

test_addl_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set addl_quants 1 -set inc_baryon_mass 1 \
		-set prefix data_temp/twop_addl -model twop -mcmc \
		> data_temp/twop_addl.scr 2> data_temp/twop_addl.err

test_crustL_nompi: bamr_nompi
	./bamr_nompi -set max_iters $(NTEST) -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set prefix data_temp/twop_crustL -model twop -mcmc \
		> data_temp/twop_crustL.scr 2> data_temp/twop_crustL.err

# ----------------------------------------------------------------------
# Documentation
# ----------------------------------------------------------------------

update-tags:
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .

doc: empty
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	git rev-parse HEAD | awk \
		'{print "`" $$1 " <http://github.com/awsteiner/bamr/tree/" $$1 ">`_"}' \
		 > sphinx/commit.rst
	cd sphinx/static; cat bib_header.txt > ../bib.rst
	cd sphinx/static; btmanip -parse bamr.bib -rst ../bib_temp.rst
	cd sphinx; cat bib_temp.rst >> bib.rst; rm -f bib_temp.rst
	cd doc; doxygen doxyfile
	cd sphinx; make html

sync-doc:
	rsync -Cavzu sphinx/build/html/* $(STATIC_DOC_DIR)/bamr

test-sync:
	rsync -Cavzun sphinx/build/html/* $(STATIC_DOC_DIR)/bamr

# ----------------------------------------------------------------------
# Misc
# ----------------------------------------------------------------------

empty:

clean:
	rm -f *.o bamr bamr_nompi process *.png

# ----------------------------------------------------------------------
mlines:
	mpirun -np 1 ./bamr -threads 1 -set prefix out/ml2 \
	-set max_time 864300 -set file_update_time 1800 \
	-set verbose 1 -set mcmc_verbose 2 \
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
	-set apply_intsc 0 -set cached_intsc 0 \
	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
	-model new_lines -method hmc -set model_dpdm 1 \
	-set emu_tov 1 \
	-initial-point-last "out/ml_67" \
	-mcmc

nlines:
	mpirun -np 1 ./bamr -threads 1 -set prefix out/nl \
	-set max_time 864300 -set file_update_time 1800 \
	-set verbose 1 -set mcmc_verbose 2 \
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
	-set apply_intsc 0 -set cached_intsc 0 \
	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
	-model new_lines -method hmc -set model_dpdm 0 \
	-initial-point-last "out/nl_56" \
	-mcmc > out/nl.log 2>&1 &

mpoly:
	mpirun -np 1 ./bamr -threads 1 -set prefix out/mp \
	-set max_time 864300 -set file_update_time 1800 \
	-set verbose 1 -set mcmc_verbose 2 \
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
	-set apply_intsc 0 -set cached_intsc 0 \
	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
	-model new_poly -method hmc -set model_dpdm 1 \
	-initial-point-last "out/mp_57" \
	-mcmc > out/mp.log 2>&1 &

npoly:
	mpirun -np 1 ./bamr -threads 1 -set prefix out/np \
	-set max_time 864300 -set file_update_time 1800 \
	-set verbose 1 -set mcmc_verbose 2 \
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
	-set apply_intsc 0 -set cached_intsc 0 \
	-set mmax_deriv 1 -set inc_pop 1 -set inc_ligo 1 \
	-model new_poly -method hmc -set model_dpdm 0 \
	-initial-point-last "out/np_56" \
	-mcmc > out/np.log 2>&1 &

include makefile.aws
