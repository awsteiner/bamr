# ----------------------------------------------------------------------
# This section has variables which may need to be modified. The rest
# of the makefile should not need too much editing.
# ----------------------------------------------------------------------

# This variable should include the directories for the O2scl, GSL, and
# HDF5 libraries. In my configuration, I use the environment variables
# O2SCL_LIB, but you can just replace this
# entire line with whatever you need.
LIB_DIRS = -L$(O2SCL_LIB) -L$(HDF5_LIB)

# This variable may need to be modified to specify the include
# directories for the GSL, boost, HDF5, and O2scl header files. If
# O2scl was installed with Eigen or armadillo support, those header
# directories may need to be here also.
INC_DIRS = -I$(O2SCL_INC) -I$(HDF5_INC)

# C++ compiler
# MPI_CXX = mpic++ 

# Generic (no MPI necessary) C++ compiler
# CXX = g++

# Comment out these two variables if you do not have GNU readline or
# if O2scl was compiled without readline support
READLINE_VAR = -DO2SCL_READLINE -DBAMR_READLINE

READLINE_LIBS = -lreadline -lncurses

# Basic compiler flags
COMPILER_FLAGS = -std=c++0x -O3

# Specify number of OpenMP threads
THREAD_VAR = -DBAMR_OMP_THREADS=2

# ----------------------------------------------------------------------
# Secondary variables
# ----------------------------------------------------------------------

ALL_FLAGS_MPI = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) \
	-DBAMR_MPI_LOAD -DO2SCL_MPI -DO2SCL_OPENMP -fopenmp \
	$(THREAD_VAR)

ALL_FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) -DBAMR_NO_MPI \
	-DBAMR_OMP_THREADS=1

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

ecsn1:
	bamr_nompi -set prefix ecsn1 -set addl_quants 1 \
		-set inc_baryon_mass 1 -set crust_from_L 1 \
		-set debug_eos 1 -model qmc_threep -mcmc 

ecsn2:
	bamr -set prefix ecsn2 -set addl_quants 1 \
		-set inc_baryon_mass 1 -set crust_from_L 1 \
		-set debug_star 1 -model qmc_threep -mcmc
	acol -read debug_star.o2 -interp gm 1.223 bm \
		-interp gm 1.237 bm \
		-interp gm 1.250 bm 

ecsn3:
	bamr -set prefix ecsn3 -set addl_quants 1 \
		-set inc_baryon_mass 1 -set crust_from_L 1 \
		-set max_iters 10 -model qmc_threep -mcmc
	acol -read ecsn3_0_out -get-row 0

# ----------------------------------------------------------------------
# Testing targets
# ----------------------------------------------------------------------

test-all: test-prep test1 test2 test3 test4 test5 test6 test7 test8 \
		test9 test10 test11 test12 test13 test14

test-all-nompi: test-prep test1-nompi test2-nompi test3-nompi \
		test4-nompi test5-nompi test6-nompi test7-nompi \
		test8-nompi test9-nompi test10-nompi test11-nompi \
		test12-nompi test13-nompi test14-nompi

test-prep: empty
	-mkdir -p data_temp
	-rm -rf data_temp/*

test1: empty
	bamr -set prefix data_temp/debug_eos \
		-set debug_eos 1 -run default.in -model twop -mcmc \
		> data_temp/debug_eos.scr
	-mv -i debug_eos.o2 data_temp

test2: empty
	bamr -set prefix data_temp/debug_star \
		-set debug_star 1 -run default.in -model twop -mcmc \
		> data_temp/debug_star.scr
	-mv -i debug_star.o2 data_temp

test3: empty
	mpirun -np 2 \
	bamr -set max_iters 300 -set prefix data_temp/twop_data \
		-run default.in -model twop -mcmc \
		> data_temp/twop_data.scr 2> data_temp/twop_data.err

test4: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set prefix data_temp/twop_nodata \
		-model twop -mcmc \
		> data_temp/twop_nodata.scr 2> data_temp/twop_nodata.err

test5: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set compute_cthick 1 \
		-set prefix data_temp/twop_cthick -model twop -mcmc \
		> data_temp/twop_cthick.scr 2> data_temp/twop_cthick.err

test6: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set prefix data_temp/fixp_nodata \
		-model fixp -mcmc \
		> data_temp/fixp_nodata.scr 2> data_temp/fixp_nodata.err

test7: empty
	mpirun -np 2 \
	bamr -set max_iters 300 -set prefix data_temp/qt_nodata \
		-model qmc_threep -mcmc \
		> data_temp/qt_nodata.scr 2> data_temp/qt_nodata.err

test8: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set prefix data_temp/qf_nodata \
		-model qmc_fixp -mcmc \
		> data_temp/qf_nodata.scr 2> data_temp/qf_nodata.err

test9: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set n_warm_up 100 \
		-set prefix data_temp/twop_warmup -model twop -mcmc \
		> data_temp/twop_warmup.scr 2> data_temp/twop_warmup.err

test10: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set prefix data_temp/twop_ai -set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 10 -mcmc \
		> data_temp/twop_ai.scr 2> data_temp/twop_ai.err

test11: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set prefix data_temp/twop_chain \
		-set max_chain_size 10 -model twop -mcmc \
		> data_temp/twop_chain.scr 2> data_temp/twop_chain.err

test12: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set prefix data_temp/twop_aic \
		-set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 20 \
		-set max_chain_size 10 -mcmc \
		> data_temp/twop_aic.scr 2> data_temp/twop_aic.err

test13: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set compute_cthick 1 -set crust_from_L 1 \
		-set addl_quants 1 -set inc_baryon_mass 1 \
		-set prefix data_temp/twop_addl -model twop -mcmc \
		> data_temp/twop_addl.scr 2> data_temp/twop_addl.err

test14: empty
	mpirun -np 2 \
	bamr -set max_iters 100 -set compute_cthick 1 -set crust_from_L 1 \
		-set prefix data_temp/twop_crustL -model twop -mcmc \
		> data_temp/twop_crustL.scr 2> data_temp/twop_crustL.err

test1-nompi: empty
	bamr_nompi -set prefix data_temp/debug_eos \
		-set debug_eos 1 -run default.in -model twop -mcmc \
		> data_temp/debug_eos.scr
	-mv -i debug_eos.o2 data_temp

test2-nompi: empty
	bamr_nompi -set prefix data_temp/debug_star \
		-set debug_star 1 -run default.in -model twop -mcmc \
		> data_temp/debug_star.scr
	-mv -i debug_star.o2 data_temp

test3-nompi: empty
	bamr_nompi -set max_iters 300 -set prefix data_temp/twop_data \
		-run default.in -model twop -mcmc \
		> data_temp/twop_data.scr 2> data_temp/twop_data.err

test4-nompi: empty
	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_nodata \
		-model twop -mcmc \
		> data_temp/twop_nodata.scr 2> data_temp/twop_nodata.err

test5-nompi: empty
	bamr_nompi -set max_iters 100 -set compute_cthick 1 \
		-set prefix data_temp/twop_cthick -model twop -mcmc \
		> data_temp/twop_cthick.scr 2> data_temp/twop_cthick.err

test6-nompi: empty
	bamr_nompi -set max_iters 100 -set prefix data_temp/fixp_nodata \
		-model fixp -mcmc \
		> data_temp/fixp_nodata.scr 2> data_temp/fixp_nodata.err

test7-nompi: empty
	bamr_nompi -set max_iters 300 -set prefix data_temp/qt_nodata \
		-model qmc_threep -mcmc \
		> data_temp/qt_nodata.scr 2> data_temp/qt_nodata.err

test8-nompi: empty
	bamr_nompi -set max_iters 100 -set prefix data_temp/qf_nodata \
		-model qmc_fixp -mcmc \
		> data_temp/qf_nodata.scr 2> data_temp/qf_nodata.err

test9-nompi: empty
	bamr_nompi -set max_iters 100 -set n_warm_up 100 \
		-set prefix data_temp/twop_warmup -model twop -mcmc \
		> data_temp/twop_warmup.scr 2> data_temp/twop_warmup.err

test10-nompi: empty
	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_ai \
		-set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 10 -mcmc \
		> data_temp/twop_ai.scr 2> data_temp/twop_ai.err

test11-nompi: empty
	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_chain \
		-set max_chain_size 10 -model twop -mcmc \
		> data_temp/twop_chain.scr 2> data_temp/twop_chain.err

test12-nompi: empty
	bamr_nompi -set max_iters 100 -set prefix data_temp/twop_aic \
		-set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 20 \
		-set max_chain_size 10 -mcmc \
		> data_temp/twop_aic.scr 2> data_temp/twop_aic.err

test13-nompi: empty
	bamr_nompi -set max_iters 100 -set compute_cthick 1 \
		-set crust_from_L 1 \
		-set addl_quants 1 -set inc_baryon_mass 1 \
		-set prefix data_temp/twop_addl -model twop -mcmc \
		> data_temp/twop_addl.scr 2> data_temp/twop_addl.err

test14-nompi: empty
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
	cp -r sphinx/build/html/* $(HOME)/wcs/int4/web/utk/bamr

# ----------------------------------------------------------------------
# Misc
# ----------------------------------------------------------------------

empty: 

clean:
	rm -f *.o bamr bamr_nompi process

