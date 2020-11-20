LIB_DIRS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
        -L$/usr/local/lib -lo2scl_eos \

INC_DIRS = -I$(EIGEN_INC) -I/usr/local/o2scl \
     		-I/usr/include/hdf5/openmpi \
         
CXX = c++ 
MPI_CXX = mpic++

BAMR_DIR = ./

COMPILER_FLAGS = $(UTKNA_CFLAGS)
COMPILER_FLAGS_MPI = $(UTKNA_MPI_CFLAGS)

PYTHON_INC_DIR = /usr/include/python3.8
PYTHON_LDFLAGS = -lpython3.8

