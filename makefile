# --------------------------------------------------------
# Variables which must be modified. These default values 
# will probably not work on your system.
# --------------------------------------------------------

# This variable should include the directories for the o2scl, gsl, and
# hdf libraries. 
LIB_DIRS = -L$(HOME)/install/lib -L$(HOME)/install/o2scl-0.915/lib

# This variable should include the parent directories for the 
# gsl and o2scl header files as well as the directories for the
# hdf5 include files
INC_DIRS = -I$(HOME)/install/include -I$(HOME)/install/include/hdf5 \
	-I$(HOME)/install/o2scl-0.915/include \
	-I$(HOME)/pkgs/Eigen-3.1.3 -I$(HOME)/install/arma/include \
	-I$(HOME)/pkgs/boost_1_53_0

# Path to MPI C++ compiler
MPI_CXX = $(HOME)/install/bin/mpic++ -std=c++0x

# Path to generic (no MPI necessary) C++ compiler
CXX = g++ -std=c++0x

# Comment out these two variables if you do not have GNU readline or
# if O2scl was compiled without readline support
READLINE_VAR = -DO2SCL_READLINE -DBAMR_READLINE

READLINE_LIBS = -lreadline -lncurses

# Some basic warning and optimization flags
COMPILER_FLAGS = -Wreturn-type -Wparentheses -Wall -Wno-unused -O3 \
	-DBAMR_MPI_LOAD

# The root include directory (only necessary if you're using the 
# plotting code, otherwise can be blank)
ROOT_INC = -I$(HOME)/root/include -I/usr/include/root

# The root library directory (only necessary if you're using the 
# plotting code, otherwise can be blank)
ROOT_LIB = -L$(HOME)/root/lib \
	-lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d \
	-lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics \
	-lMathCore -lThread -lz -pthread -ldl

# The location of graph.cpp (only necessary if you're using the 
# plotting code, otherwise can be blank)
GRAPH_CPP_DIR = $(HOME)/o2scl/src/other

# --------------------------------------------------------
# Basic bamr targets
# --------------------------------------------------------

FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR)

LIB = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5 -lgsl -lgslcblas -lm $(READLINE_LIBS)

bamr: bamr.o entry.o models.o misc.o main.o
	$(MPI_CXX) $(FLAGS) $(LIB_DIRS) -o bamr main.o misc.o \
		entry.o models.o bamr.o $(LIB) 

main.o: main.cpp 
	$(MPI_CXX) $(FLAGS) -o main.o -c main.cpp

misc.o: misc.cpp 
	$(MPI_CXX) $(FLAGS) -o misc.o -c misc.cpp

models.o: models.cpp 
	$(MPI_CXX) $(FLAGS) -o models.o -c models.cpp

entry.o: entry.cpp 
	$(MPI_CXX) $(FLAGS) -o entry.o -c entry.cpp

bamr.o: bamr.cpp 
	$(MPI_CXX) $(FLAGS) -o bamr.o -c bamr.cpp

test:
	bamr -run default.in -model twop -mcmc run1 &

# --------------------------------------------------------
# Plotting targets
# --------------------------------------------------------

$(GRAPH_CPP_DIR)/graph.o: $(GRAPH_CPP_DIR)/graph.cpp
	cd $(GRAPH_CPP_DIR); $(CXX) $(FLAGS) $(ROOT_INC) -c graph.cpp

plot.o: plot.cpp 
	$(CXX) $(FLAGS) $(ROOT_INC) -c plot.cpp

plot: plot.o $(GRAPH_CPP_DIR)/graph.o
	$(CXX) $(FLAGS) -o plot plot.o $(GRAPH_CPP_DIR)/graph.o\
		$(ROOT_LIB) $(LIB_DIRS) $(LIB) 

plot2d.o: plot2d.cpp 
	$(CXX) $(FLAGS) $(ROOT_INC) -c plot2d.cpp

plot2d: plot2d.o $(GRAPH_CPP_DIR)/graph.o
	$(CXX) $(FLAGS) -o plot2d plot2d.o $(GRAPH_CPP_DIR)/graph.o\
		$(ROOT_LIB) $(LIB_DIRS) $(LIB) 

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
	rm -f *.o bamr plot plot2d


