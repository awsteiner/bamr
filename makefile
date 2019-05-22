# ----------------------------------------------------------------------
# This section has variables which may need to be modified. The rest
# of the makefile should not need too much editing.
# ----------------------------------------------------------------------

# This variable should include the directories for the O2scl, GSL, and
# HDF5 libraries.
LIB_DIRS = 

# This variable may need to be modified to specify the include
# directories for the GSL, Boost, HDF5, and O2scl header files. 
INC_DIRS = 

# C++ compiler (e.g. mpicxx)
# MPI_CXX = 

# Generic (no MPI necessary) C++ compiler (e.g. g++)
# CXX = 

# Set these two variables to be empty if you do not have GNU readline
# readline support
READLINE_VAR = -DBAMR_READLINE

READLINE_LIBS = -lreadline -lncurses

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

endif

# ----------------------------------------------------------------------
# Secondary variables
# ----------------------------------------------------------------------

ALL_FLAGS_MPI = $(COMPILER_MPI_FLAGS) $(INC_DIRS) $(READLINE_VAR) \
	-DBAMR_MPI -DO2SCL_MPI -DO2SCL_OPENMP -fopenmp 

ALL_FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) 

LIB = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5_hl -lhdf5 -lgsl -lgslcblas -lm $(READLINE_LIBS)

# ----------------------------------------------------------------------
# Targets for bamr
# ----------------------------------------------------------------------

bamr: bamr_class.o models.o nstar_cold2.o main.o mcmc_bamr.o ns_data.o
	$(MPI_CXX) $(ALL_FLAGS_MPI) $(LIB_DIRS) -o bamr main.o \
		nstar_cold2.o ns_data.o models.o mcmc_bamr.o \
		bamr_class.o $(LIB) 

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

mcmc_bamr.o: mcmc_bamr.cpp mcmc_bamr.h models.o main.o nstar_cold2.o
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o mcmc_bamr.o -c mcmc_bamr.cpp

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

bamr_nompi: bamr_class_nompi.o models_nompi.o \
		nstar_cold2_nompi.o main_nompi.o mcmc_bamr_nompi.o \
		ns_data_nompi.o 
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o bamr_nompi \
		bamr_class_nompi.o models_nompi.o \
		nstar_cold2_nompi.o main_nompi.o mcmc_bamr_nompi.o \
		ns_data_nompi.o $(LIB) 

main_nompi.o: main.cpp
	$(CXX) $(ALL_FLAGS) -o main_nompi.o -c main.cpp

nstar_cold2_nompi.o: nstar_cold2.cpp nstar_cold2.h
	$(CXX) $(ALL_FLAGS) -o nstar_cold2_nompi.o -c nstar_cold2.cpp

models_nompi.o: models.cpp models.h
	$(CXX) $(ALL_FLAGS) -o models_nompi.o -c models.cpp

ns_data_nompi.o: ns_data.cpp ns_data.h
	$(CXX) $(ALL_FLAGS) -o ns_data_nompi.o -c ns_data.cpp

mcmc_bamr_nompi.o: mcmc_bamr.cpp mcmc_bamr.h
	$(CXX) $(ALL_FLAGS) -o mcmc_bamr_nompi.o -c mcmc_bamr.cpp

bamr_class_nompi.o: bamr_class.cpp bamr_class.h
	$(CXX) $(ALL_FLAGS) -o bamr_class_nompi.o -c bamr_class.cpp

# ----------------------------------------------------------------------
# Targets for process
# ----------------------------------------------------------------------

process.o: process.cpp process.h
	$(CXX) $(ALL_FLAGS) -o process.o -c process.cpp

process_main.o: process_main.cpp
	$(CXX) $(ALL_FLAGS) -o process_main.o -c process_main.cpp

process: process.o process_main.o
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o process process_main.o \
		process.o $(LIB) 

# ----------------------------------------------------------------------
# Targets for test
# ----------------------------------------------------------------------

test.o: test.cpp
	$(CXX) $(ALL_FLAGS) -o test.o -c test.cpp

test: test.o
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o test \
		test.o $(LIB) 

# ----------------------------------------------------------------------
# Testing targets
# ----------------------------------------------------------------------

test_all: test_prep test_data test_nodata test_cthick test_fixp \
		test_qt test_qf test_warmup test_ai \
		test_addl test_crustL test_tableseq test_storej bamr \
		test_debug_eos test_debug_star test_rejtab

# test11 test12 

test_all_nompi: test_prep test_debug_eos_nompi test_debug_star_nompi \
		test_data_nompi test_nodata_nompi test_cthick_nompi \
		test_fixp_nompi test_qt_nompi test_qf_nompi \
		test_warmup_nompi test_ai_nompi \
		test_addl_nompi test_crustL_nompi bamr_nompi

#test11_nompi test12_nompi 

test_prep: bamr
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

test_debug_eos_nompi: bamr_nompi
	bamr_nompi -set prefix data_temp/debug_eos \
		-set debug_eos 1 -run default.in -model twop -mcmc \
		> data_temp/debug_eos.scr
	-mv -i debug_eos.o2 data_temp

test_debug_star_nompi: bamr_nompi
	bamr_nompi -set prefix data_temp/debug_star \
		-set debug_star 1 -run default.in -model twop -mcmc \
		> data_temp/debug_star.scr
	-mv -i debug_star.o2 data_temp

test_data_nompi: bamr_nompi
	bamr_nompi -set max_iters 300 -set prefix data_temp/twop_data \
		-run default.in -model twop -mcmc \
		> data_temp/twop_data.scr 2> data_temp/twop_data.err

test_nodata_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_nodata \
		-model twop -mcmc \
		> data_temp/twop_nodata.scr 2> data_temp/twop_nodata.err

test_cthick_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set compute_cthick 1 \
		-set prefix data_temp/twop_cthick -model twop -mcmc \
		> data_temp/twop_cthick.scr 2> data_temp/twop_cthick.err

test_fixp_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set prefix data_temp/fixp_nodata \
		-model fixp -mcmc \
		> data_temp/fixp_nodata.scr 2> data_temp/fixp_nodata.err

test_qt_nompi: bamr_nompi
	bamr_nompi -set max_iters 300 -set prefix data_temp/qt_nodata \
		-model qmc_threep -mcmc \
		> data_temp/qt_nodata.scr 2> data_temp/qt_nodata.err

test_qf_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set prefix data_temp/qf_nodata \
		-model qmc_fixp -mcmc \
		> data_temp/qf_nodata.scr 2> data_temp/qf_nodata.err

test_warmup_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set n_warm_up 100 \
		-set prefix data_temp/twop_warmup -model twop -mcmc \
		> data_temp/twop_warmup.scr 2> data_temp/twop_warmup.err

test_ai_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_ai \
		-set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 10 -mcmc \
		> data_temp/twop_ai.scr 2> data_temp/twop_ai.err

# These are unnecessary now because there is no max_chain_size setting

# test11_nompi: bamr_nompi
# 	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_chain \
# 		-set max_chain_size 10 -model twop -mcmc \
# 		> data_temp/twop_chain.scr 2> data_temp/twop_chain.err

# test12_nompi: bamr_nompi
# 	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_aic \
# 		-set aff_inv 1 \
# 		-set step_fac 2.0 -model twop -set n_walk 20 \
# 		-set max_chain_size 10 -mcmc \
# 		> data_temp/twop_aic.scr 2> data_temp/twop_aic.err

test_addl_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set addl_quants 1 -set inc_baryon_mass 1 \
		-set prefix data_temp/twop_addl -model twop -mcmc \
		> data_temp/twop_addl.scr 2> data_temp/twop_addl.err

test_crustL_nompi: bamr_nompi
	bamr_nompi -set max_iters 100 -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set prefix data_temp/twop_crustL -model twop -mcmc \
		> data_temp/twop_crustL.scr 2> data_temp/twop_crustL.err

# ----------------------------------------------------------------------
# Documentation
# ----------------------------------------------------------------------

update-tags:
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .

doc: empty
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .
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
	rm -f *.o bamr bamr_nompi process

