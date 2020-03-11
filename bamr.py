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

# Create the bamr pointers
bamr_lib.create_pointers.argtypes=[void_pp,void_pp]
bamr_class_ptr=ctypes.c_void_p()
model_data_ptr=ctypes.c_void_p()
bamr_ptr=bamr_lib.create_pointers(bamr_class_ptr,model_data_ptr)

# Compute one point
cp_fun=bamr_lib.py_compute_point
cp_fun.argtypes=[ctypes.c_void_p,ctypes.c_void_p,
                 ctypes.c_int,ctypes.POINTER(ctypes.c_double)]
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

# Free the memory from the bamr pointers
bamr_lib.destroy_pointers.argtypes=[ctypes.c_void_p,ctypes.c_void_p]
bamr_lib.destroy_pointers(bamr_class_ptr,model_data_ptr)

