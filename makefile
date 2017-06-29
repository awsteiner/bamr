# ----------------------------------------------------------------------
# This section has variables which may need to be modified. The rest
# of the makefile should not need too much editing.
# ----------------------------------------------------------------------

# This variable should include the directories for the O2scl, GSL, and
# HDF5 libraries. In my configuration, I use the environment variables
# O2SCL_LIB, HDF5_LIB, and GSL_LIB, but you can just replace this
# entire line with whatever you need.
LIB_DIRS = -L$(O2SCL_LIB) -L$(HDF5_LIB) -L$(GSL_LIB)

# This variable should include the parent directories for the GSL,
# boost, HDF5, and O2scl header files. If O2scl was installed with
# Eigen or armadillo support, those header directories may need to be
# here also.
INC_DIRS = -I$(O2SCL_INC) -I$(HDF5_INC) -I$(BOOST_INC) -I$(GSL_INC) \
	-I$(EIGEN_INC) -I$(ARMA_INC)

# C++ compiler
# MPI_CXX = mpic++ 

# Generic (no MPI necessary) C++ compiler
# CXX = g++

# Comment out these two variables if you do not have GNU readline or
# if O2scl was compiled without readline support
READLINE_VAR = -DO2SCL_READLINE -DBAMR_READLINE

READLINE_LIBS = -lreadline -lncurses

# Basic optimization flags
COMPILER_FLAGS = -std=c++0x -O3

# ----------------------------------------------------------------------
# Secondary variables
# ----------------------------------------------------------------------

ALL_FLAGS_MPI = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) -DBAMR_MPI_LOAD

ALL_FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) -DBAMR_NO_MPI

LIB = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5_hl -lhdf5 -lgsl -lgslcblas -lm $(READLINE_LIBS)

# ----------------------------------------------------------------------
# Targets for bamr
# ----------------------------------------------------------------------

bamr: bamr_class.o models.o nstar_cold2.o main.o mcmc_bamr.o
	$(MPI_CXX) $(ALL_FLAGS_MPI) $(LIB_DIRS) -o bamr main.o nstar_cold2.o \
		models.o mcmc_bamr.o bamr_class.o $(LIB) 

main.o: main.cpp
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o main.o -c main.cpp

nstar_cold2.o: nstar_cold2.cpp nstar_cold2.h
	$(MPI_CXX) $(ALL_FLAGS_MPI) -o nstar_cold2.o -c nstar_cold2.cpp

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

bamr_nompi: bamr_nompi.o models_nompi.o \
		nstar_cold2_nompi.o main_nompi.o
	$(CXX) $(ALL_FLAGS) $(LIB_DIRS) -o bamr_nompi main_nompi.o \
		nstar_cold2_nompi.o models_nompi.o \
		bamr_nompi.o $(LIB) 

main_nompi.o: main.cpp
	$(CXX) $(ALL_FLAGS) -o main_nompi.o -c main.cpp

nstar_cold2_nompi.o: nstar_cold2.cpp nstar_cold2.h
	$(CXX) $(ALL_FLAGS) -o nstar_cold2_nompi.o -c nstar_cold2.cpp

models_nompi.o: models.cpp models.h
	$(CXX) $(ALL_FLAGS) -o models_nompi.o -c models.cpp

bamr_nompi.o: bamr.cpp bamr.h models_nompi.o
		nstar_cold2_nompi.o main_nompi.o
	$(CXX) $(ALL_FLAGS) -o bamr_nompi.o -c bamr.cpp

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
# test
# ----------------------------------------------------------------------

test:
	mpirun -np 2 bamr -run default.in -model twop -mcmc

test_all:
	-mkdir -p data_temp
	-rm -rf data_temp/*
	bamr -set prefix data_temp/debug_eos \
		-set debug_eos 1 -run default.in -model twop -mcmc \
		> data_temp/debug_eos.scr
	-mv -i debug_eos.o2 data_temp
	bamr -set prefix data_temp/debug_star \
		-set debug_star 1 -run default.in -model twop -mcmc \
		> data_temp/debug_star.scr
	-mv -i debug_star.o2 data_temp
	bamr -set max_iters 300 -set prefix data_temp/twop_data \
		-run default.in -model twop -mcmc \
		> data_temp/twop_data.scr 2> data_temp/twop_data.err
	bamr -set max_iters 100 -set prefix data_temp/twop_nodata \
		-model twop -mcmc \
		> data_temp/twop_nodata.scr 2> data_temp/twop_nodata.err
	bamr -set max_iters 100 -set compute_cthick 1 -set crust_from_L 1 \
		-set prefix data_temp/twop_crustL -model twop -mcmc \
		> data_temp/twop_crustL.scr 2> data_temp/twop_crustL.err
	bamr -set max_iters 100 -set compute_cthick 1 \
		-set prefix data_temp/twop_cthick -model twop -mcmc \
		> data_temp/twop_cthick.scr 2> data_temp/twop_cthick.err
	bamr -set max_iters 100 -set compute_cthick 1 -set crust_from_L 1 \
		-set addl_quants 1 -set inc_baryon_mass 1 \
		-set prefix data_temp/twop_addl -model twop -mcmc \
		> data_temp/twop_addl.scr 2> data_temp/twop_addl.err
	bamr -set max_iters 100 -set prefix data_temp/fixp_nodata \
		-model fixp -mcmc \
		> data_temp/fixp_nodata.scr 2> data_temp/fixp_nodata.err
	bamr -set max_iters 300 -set prefix data_temp/qt_nodata \
		-model qmc_threep -mcmc \
		> data_temp/qt_nodata.scr 2> data_temp/qt_nodata.err
	bamr -set max_iters 100 -set prefix data_temp/qf_nodata \
		-model qmc_fixp -mcmc \
		> data_temp/qf_nodata.scr 2> data_temp/qf_nodata.err
	bamr -set max_iters 100 -set n_warm_up 100 \
		-set prefix data_temp/twop_warmup -model twop -mcmc \
		> data_temp/twop_warmup.scr 2> data_temp/twop_warmup.err
	bamr -set max_iters 100 -set prefix data_temp/twop_ai -set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 10 -mcmc \
		> data_temp/twop_ai.scr 2> data_temp/twop_ai.err
	bamr -set max_iters 100 -set prefix data_temp/twop_chain \
		-set max_chain_size 10 -model twop -mcmc \
		> data_temp/twop_chain.scr 2> data_temp/twop_chain.err
	bamr -set max_iters 100 -set prefix data_temp/twop_aic -set aff_inv 1 \
		-set step_fac 2.0 -model twop -set n_walk 20 \
		-set max_chain_size 10 -mcmc \
		> data_temp/twop_aic.scr 2> data_temp/twop_aic.err

# ----------------------------------------------------------------------
# Internal 
# ----------------------------------------------------------------------

empty:

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

clean:
	rm -f *.o bamr bamr_nompi process

