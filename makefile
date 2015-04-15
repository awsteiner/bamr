# --------------------------------------------------------
# Variables which must be modified. These default values 
# will probably not work on your system.
# --------------------------------------------------------

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
MPI_CXX = mpic++ 

# Generic (no MPI necessary) C++ compiler
CXX = g++-4.9

# Comment out these two variables if you do not have GNU readline or
# if O2scl was compiled without readline support
READLINE_VAR = -DO2SCL_READLINE -DBAMR_READLINE

READLINE_LIBS = -lreadline -lncurses

# Some basic warning and optimization flags
COMPILER_FLAGS = -Wreturn-type -Wparentheses -Wall -Wno-unused -O3 \
	-std=c++0x -Wno-literal-suffix

# --------------------------------------------------------
# Basic bamr targets
# --------------------------------------------------------

FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) -DBAMR_MPI_LOAD
FLAGS2 = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) -DBAMR_NO_MPI

LIB = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5 -lgsl -lgslcblas -lm $(READLINE_LIBS)

bamr: bamr.o entry.o models.o nstar_cold2.o main.o
	$(MPI_CXX) $(FLAGS) $(LIB_DIRS) -o bamr main.o nstar_cold2.o \
		entry.o models.o bamr.o $(LIB) 

main.o: main.cpp 
	$(MPI_CXX) $(FLAGS) -o main.o -c main.cpp

nstar_cold2.o: nstar_cold2.cpp 
	$(MPI_CXX) $(FLAGS) -o nstar_cold2.o -c nstar_cold2.cpp

models.o: models.cpp 
	$(MPI_CXX) $(FLAGS) -o models.o -c models.cpp

entry.o: entry.cpp 
	$(MPI_CXX) $(FLAGS) -o entry.o -c entry.cpp

bamr.o: bamr.cpp 
	$(MPI_CXX) $(FLAGS) -o bamr.o -c bamr.cpp

bamr_nompi: bamr_nompi.o entry_nompi.o models_nompi.o \
		nstar_cold2_nompi.o main_nompi.o
	$(CXX) $(FLAGS2) $(LIB_DIRS) -o bamr_nompi main_nompi.o \
		nstar_cold2_nompi.o entry_nompi.o models_nompi.o \
		bamr_nompi.o $(LIB) 

main_nompi.o: main.cpp 
	$(CXX) $(FLAGS2) -o main_nompi.o -c main.cpp

nstar_cold2_nompi.o: nstar_cold2.cpp 
	$(CXX) $(FLAGS2) -o nstar_cold2_nompi.o -c nstar_cold2.cpp

models_nompi.o: models.cpp 
	$(CXX) $(FLAGS2) -o models_nompi.o -c models.cpp

entry_nompi.o: entry.cpp 
	$(CXX) $(FLAGS2) -o entry_nompi.o -c entry.cpp

bamr_nompi.o: bamr.cpp 
	$(CXX) $(FLAGS2) -o bamr_nompi.o -c bamr.cpp

test:
	bamr -run default.in -model twop -mcmc run1

process.o: process.cpp 
	$(CXX) $(FLAGS2) -o process.o -c process.cpp

process: process.o
	$(CXX) $(FLAGS2) $(LIB_DIRS) -o process process.o $(LIB) 

# --------------------------------------------------------
# Internal 
# --------------------------------------------------------

empty:

VERSION = 0.2

dist:	
	mkdir -p bamr-$(VERSION)
	mkdir -p bamr-$(VERSION)/doc
	mkdir -p bamr-$(VERSION)/doc/html
	mkdir -p bamr-$(VERSION)/doc/latex
	cp *.cpp *.h makefile bamr-$(VERSION)
	cp doc/latex/refman.pdf bamr-$(VERSION)/doc/latex
	cp -r doc/html/*.html doc/html/*.png doc/html/*.js \
		doc/html/search doc/html/*.css bamr-$(VERSION)/doc/html
	tar cvzf bamr-$(VERSION).tar.gz bamr-$(VERSION)

dist-clean:
	rm -rf bamr-$(VERSION).tar.gz bamr-$(VERSION)

sf-web:
	cd doc/html; rsync -Cavzu * \
		awsteiner,bamr@web.sourceforge.net:htdocs

update-tags:
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .

doc: empty
	cd doc; doxygen doxyfile
	cat doc/doxygen.log
#	cd doc; perl rep.perl < html/search/search.css > temp.css
#	cd doc; mv temp.css html/search/search.css

docp: empty
	cd doc/latex; $(MAKE)

clean:
	rm -f *.o bamr bamr_nompi


