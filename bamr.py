#!/usr/bin/env python3

import ctypes
import platform

if platform.system()=='Darwin':
    bamr_lib=ctypes.CDLL('libbamr.dylib')
else:
    bamr_lib=ctypes.CDLL('libbamr.so')

bamr_lib.create_bamr_class.restype=ctypes.c_void_p
bamr_ptr=bamr_lib.create_bamr_class()

bamr_lib.destroy_bamr_class.argtypes=[ctypes.c_void_p]
bamr_lib.destroy_bamr_class(bamr_ptr)

