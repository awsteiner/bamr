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
CXX = g++

# Comment out these two variables if you do not have GNU readline or
# if O2scl was compiled without readline support
READLINE_VAR = -DO2SCL_READLINE -DBAMR_READLINE

READLINE_LIBS = -lreadline -lncurses

# Some basic warning and optimization flags
COMPILER_FLAGS = -Wreturn-type -Wparentheses -Wall -Wno-unused -O3 \
	-DBAMR_MPI_LOAD -DO2SCL_NO_TR1_MEMORY -std=c++0x

# --------------------------------------------------------
# Basic bamr targets
# --------------------------------------------------------

FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR)

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

test:
	bamr -run default.in -model twop -mcmc run1

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

update_tag:
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
	rm -f *.o bamr


