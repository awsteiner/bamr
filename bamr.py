#!/usr/bin/env python3

import ctypes
import platform

# Load ctypes library
if platform.system()=='Darwin':
    bamr_lib=ctypes.CDLL('libbamr.dylib')
else:
    bamr_lib=ctypes.CDLL('./libbamr.so')

# Type definitions
void_pp=ctypes.POINTER(ctypes.c_void_p)
double_ptr=ctypes.POINTER(ctypes.c_double)
double_ptr_ptr=ctypes.POINTER(double_ptr)
int_ptr=ctypes.POINTER(ctypes.c_int)

# Create the bamr pointers
bamr_lib.create_pointers.argtypes=[void_pp,void_pp]
bamr_class_ptr=ctypes.c_void_p()
model_data_ptr=ctypes.c_void_p()
bamr_ptr=bamr_lib.create_pointers(bamr_class_ptr,model_data_ptr)

# Compute one point
cp_fun=bamr_lib.py_compute_point
cp_fun.argtypes=[ctypes.c_void_p,ctypes.c_void_p,
                 ctypes.c_int,double_ptr]
point=(ctypes.c_double * 8)()
point[0]=1.0
point[1]=-3.0
point[2]=0.165
point[3]=0.644
point[4]=1.51
point[5]=0.576
point[6]=4.6
point[7]=1.21
cp_fun(bamr_class_ptr,model_data_ptr,8,point)

# Get a column from the table
gmc_fun=bamr_lib.get_mvsr_column
col_name=ctypes.c_char_p(b'gm')
nrows=ctypes.c_int(0)
ptr=double_ptr()
gmc_fun.argtypes=[ctypes.c_void_p,ctypes.c_char_p,
                 int_ptr,double_ptr_ptr]
gmc_fun.restype=ctypes.c_int
gmc_ret=gmc_fun(model_data_ptr,col_name,ctypes.byref(nrows),
                ctypes.byref(ptr))

for i in range(0,nrows.value):
    print('%d %7.6e' % (i,ptr[i]))

# Free the memory from the bamr pointers
bamr_lib.destroy_pointers.argtypes=[ctypes.c_void_p,ctypes.c_void_p]
bamr_lib.destroy_pointers(bamr_class_ptr,model_data_ptr)

