#!/usr/bin/env python3

import ctypes
import platform

if platform.system()=='Darwin':
    bamr_lib=ctypes.CDLL('libbamr.dylib')
else:
    bamr_lib=ctypes.CDLL('libbamr.so')

void_pp=ctypes.POINTER(ctypes.c_void_p)

bamr_lib.create_pointers.argtypes=[void_pp,void_pp]
bamr_class_ptr=ctypes.c_void_p()
model_data_ptr=ctypes.c_void_p()
bamr_ptr=bamr_lib.create_pointers(bamr_class_ptr,model_data_ptr)

bamr_lib.destroy_pointers.argtypes=[ctypes.c_void_p,ctypes.c_void_p]
bamr_lib.destroy_pointers(bamr_class_ptr,model_data_ptr)

