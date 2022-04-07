
LCXX = $(CXX)

# ----------------------------------------------------------------------
# UTK generic makefile
# ----------------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration
LIBS = $(UTKNA_O2SCL_LIBS)
LCXX = $(UTKNA_CXX)
LMPI_CXX = $(UTKNA_MPI_CXX)
BAMR_DIR = $(UTKNA_BAMR_DIR)
LCFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) \
        -I$(BAMR_DIR)
LMPI_CFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_MPI_CFLAGS) \
        -I$(BAMR_DIR) -DO2SCL_MPI -DBAMR_MPI

endif

all: 
	$(LCXX) $(LCFLAGS) -o like likelihood.cpp $(LIBS)

like:
	$(LCXX) $(LCFLAGS) -o like likelihood.cpp $(LIBS)

test:
	$(LCXX) $(LCFLAGS) -o test test.cpp likelihood.cpp mass_data.cpp $(LIBS)

plot:	
	gnuplot "likelihood.gnu"

clean:	
	rm -f like test *.cpp~

