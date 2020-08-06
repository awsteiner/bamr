#  -------------------------------------------------------------------
#  
#  Copyright (C) 2020, Andrew W. Steiner and Sarah Wellence
#  
#  This file is part of bamr.
#  
#  bamr is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  bamr is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with bamr. If not, see <http://www.gnu.org/licenses/>.
#  
#  -------------------------------------------------------------------

import ctypes
import platform
import math
import numpy as np

def force_bytes(obj):
    """
    This function returns the bytes object corresponding to ``obj``
    in case it is a string using UTF-8. 
    """
    if isinstance(obj,np.bytes_)==False and isinstance(obj,bytes)==False:
        return bytes(obj,'utf-8')
    return obj

class bamr_py:
    """
    Wrapper class for accessing and using bamr

    The use of this class follows a particular order. The settings()
    function must be called first, then optionally add_data(), then
    bamr_init(), then compute_point().

    bp=bamr.bamr_py(b'twop')
    bp.settings()
    bp.add_data(...)
    bp.bamr_init()
    bp.compute_point()

    This wrapper does not currently support any parallelization.
    """

    bamr_lib=0
    """
    bamr library object (set in __init())
    """

    bamr_class_ptr=0
    """ 
    Pointer to the bamr_class object (set in __init())
    """
    
    model_data_ptr=0
    """
    Pointer to the model_data object (set in __init())
    """
    
    ns_data_ptr=0
    """
    Pointer to the ns_data object (set in __init())
    """
    
    settings_ptr=0
    """
    Pointer to the settings object (set in __init())
    """

    n_params=0
    """ 
    Number of model parameters
    """
    
    param_names=[]
    """
    Names of model parameters
    """
    
    param_units=[]
    """
    Units of model parameters
    """
    
    init_called=False
    """ 
    True only if bamr_init() has been called
    """
    
    settings_called=False
    """ 
    True only if settings() has been called
    """

    apply_intsc=False
    """
    Desc
    """
    
    def __free_pointers(self):
        """
        Private method to free bamr pointers
        """
        self.bamr_lib.destroy_pointers.argtypes=[ctypes.c_void_p,
                                                 ctypes.c_void_p]
        self.bamr_lib.destroy_pointers(self.bamr_class_ptr,
                                       self.model_data_ptr)
        self.init_called=False
        self.settings_called=False
        self.bamr_class_ptr=0
        self.model_data_ptr=0
        # End of __free_pointers()
        return

    def __make_pointers(self,model=b'twop',data_dir=b'data',verbose=1):
        """
        Desc
        """

        if force_bytes(model)==b'twop':
            self.n_params=8
            self.param_names=["comp","kprime","esym","gamma",
                              "trans1","exp1","trans2","exp2"];
            self.param_units=["1/fm","1/fm","1/fm","","1/fm^4","",
                              "1/fm^4",""];
        elif force_bytes(model)==b'tews_threep_ligo':
            self.n_params=12
            self.param_names=["a","alpha","param_S","param_L",
                              "index1","trans1","index2","trans2",
                              "index3","M_chirp_det","eta","z_cdf"]
            self.param_units=["MeV","","MeV","MeV",
                              "","1/fm^4","","1/fm^4",
                              "","Msun","",""]
        elif force_bytes(model)==b'tews_fixp_ligo':
            self.n_params=11
            self.param_names=["a","alpha","param_S","param_L",
                              "pres1","pres2","pres3","pres4",
                              "M_chirp_det","eta","z_cdf"]
            self.param_units=["MeV","","MeV","MeV",
                              "","1/fm^4","","1/fm^4",
                              "","Msun","",""]
        else:
            print('Unknown model type',model)
            quit()

        if len(self.param_names)!=self.n_params:
            print('Problem in parameter names.')
            quit()
        if len(self.param_units)!=self.n_params:
            print('Problem in parameter units.')
            quit()
            
        # Type definitions
        void_pp=ctypes.POINTER(ctypes.c_void_p)
        double_ptr=ctypes.POINTER(ctypes.c_double)
        double_ptr_ptr=ctypes.POINTER(double_ptr)
        int_ptr=ctypes.POINTER(ctypes.c_int)

        # Create the bamr pointers
        self.bamr_lib.create_pointers.argtypes=[ctypes.c_char_p,void_pp,
                                                void_pp,void_pp,void_pp,
                                                ctypes.c_char_p,
                                                ctypes.c_int]
        self.bamr_class_ptr=ctypes.c_void_p()
        self.model_data_ptr=ctypes.c_void_p()
        self.ns_data_ptr=ctypes.c_void_p()
        self.settings_ptr=ctypes.c_void_p()
        model_c=ctypes.c_char_p(force_bytes(model))
        data_dir_c=ctypes.c_char_p(force_bytes(data_dir))
        if verbose>2:
            print('Going to create_pointers() function.')
        bamr_ptr=self.bamr_lib.create_pointers(model_c,
                                               self.bamr_class_ptr,
                                               self.model_data_ptr,
                                               self.ns_data_ptr,
                                               self.settings_ptr,
                                               data_dir_c,
                                               verbose)

        if verbose>2:
            print('Done with create_pointers() function.')

        # End of __make_pointers()
        return
        
    def __init__(self,model=b'twop',data_dir=b'data',verbose=1,openmp=False):
        """ 
        Load bamr_lib using ctypes and create the associated bamr
        pointers given the specified model
        """

        # Load ctypes library
        if platform.system()=='Darwin':
            self.bamr_lib=ctypes.CDLL('libbamr.dylib',
                                      mode=ctypes.RTLD_GLOBAL)
        else:
            if openmp:
                print('Loading openmp')
                self.openmp_lib=ctypes.CDLL('/usr/lib/gcc/x86_64-'+
                                            'linux-gnu/9/libgomp.so',
                                            mode=ctypes.RTLD_GLOBAL)
            self.bamr_lib=ctypes.CDLL('libbamr.so',mode=ctypes.RTLD_GLOBAL)

        if verbose>2:
            print('Loaded libbamr.')
            
        self.__make_pointers(model,data_dir,verbose)

        if verbose>2:
            print('Done in __init__().')

        # End of __init__()
        return

    def bamr_init(self):
        """
        """
        
        if self.bamr_class_ptr==0:
            print('Bamr class pointer is 0.')
            quit()
        if self.settings_called==False:
            print('The settings() function has not been called.')
            quit()
        
        double_ptr=ctypes.POINTER(ctypes.c_double)
        double_ptr_ptr=ctypes.POINTER(double_ptr)
        int_ptr=ctypes.POINTER(ctypes.c_int)
        int_ptr_ptr=ctypes.POINTER(int_ptr)
        char_ptr=ctypes.POINTER(ctypes.c_char)
        char_ptr_ptr=ctypes.POINTER(char_ptr)
        
        self.bamr_lib.init.argtypes=[ctypes.c_void_p,ctypes.c_void_p,
                                     ctypes.c_void_p,ctypes.c_void_p,
                                     int_ptr,int_ptr_ptr,char_ptr_ptr,
                                     int_ptr_ptr,char_ptr_ptr,
                                     double_ptr_ptr,double_ptr_ptr]
        npar=ctypes.c_int(0)
        name_counts=int_ptr()
        name_c=char_ptr()
        unit_counts=int_ptr()
        unit_c=char_ptr()
        low=double_ptr()
        high=double_ptr()
        
        self.bamr_lib.init.restype=ctypes.c_int
        iret=self.bamr_lib.init(self.bamr_class_ptr,
                                self.model_data_ptr,
                                self.ns_data_ptr,
                                self.settings_ptr,
                                ctypes.byref(npar),
                                ctypes.byref(name_counts),
                                ctypes.byref(name_c),
                                ctypes.byref(unit_counts),
                                ctypes.byref(unit_c),
                                ctypes.byref(low),
                                ctypes.byref(high))

        npar=npar.value
        nix=0
        uix=0
        names=[]
        units=[]
        for i in range(0,npar):
            tname=b''
            for j in range(0,name_counts[i]):
                tname=tname+name_c[nix]
                nix=nix+1
            tname=tname.decode('utf-8')
            names.append(tname)
            tunit=b''
            for j in range(0,unit_counts[i]):
                tunit=tunit+unit_c[uix]
                uix=uix+1
            tunit=tunit.decode('utf-8')
            units.append(tunit)

        low=[low[i] for i in range(0,npar)]
        high=[high[i] for i in range(0,npar)]
            
        if iret!=0:
            print('Function bamr_init() failed.')
            quit()
        self.init_called=True
        return iret,npar,names,units,low,high

    def settings(self,inc_baryon_mass=False,addl_quants=False,verbose=0,
                 norm_max=False,crust_from_L=False,
                 compute_cthick=False,apply_intsc=True,
                 cached_intsc=True,prior_eta=True,data_dir='data'):
                 
        """
        Apply various settings (must be called before init())
        """
        
        if self.bamr_class_ptr==0:
            print('Bamr class pointer is 0.')
            quit()
        if self.bamr_lib==0:
            print('Bamr lib object is 0.')
            quit()
            
        self.bamr_lib.set_parameter.argtypes=[ctypes.c_void_p,
                                              ctypes.c_void_p,
                                              ctypes.c_char_p,
                                              ctypes.c_double]
        self.bamr_lib.set_parameter_string.argtypes=[ctypes.c_void_p,
                                                     ctypes.c_void_p,
                                                     ctypes.c_char_p,
                                                     ctypes.c_char_p]
        if inc_baryon_mass:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'inc_baryon_mass'),
                                        1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'inc_baryon_mass'),
                                        0.0)
        if addl_quants:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'addl_quants'),1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'addl_quants'),0.0)
        if norm_max:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'norm_max'),
                                        1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'norm_max'),
                                        0.0)
        if crust_from_L:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'crust_from_L'),
                                        1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'crust_from_L'),
                                        0.0)
        if compute_cthick:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'compute_cthick'),
                                        1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'compute_cthick'),
                                        0.0)
        if apply_intsc:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'apply_intsc'),
                                        1.0)
            self.apply_intsc=True
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'apply_intsc'),
                                        0.0)
            self.apply_intsc=False
        if cached_intsc:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'cached_intsc'),
                                        1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'cached_intsc'),
                                        0.0)
        if prior_eta:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'prior_eta'),
                                        1.0)
        else:
            self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                        ctypes.c_char_p(b'prior_eta'),
                                        0.0)
        self.bamr_lib.set_parameter(self.bamr_class_ptr,self.settings_ptr,
                                    ctypes.c_char_p(b'verbose'),
                                    float(verbose))
        if data_dir:
            fun=self.bamr_lib.set_parameter_string
            print('Setting data_dir to ',data_dir)
            fun(self.bamr_class_ptr,
                self.settings_ptr,
                ctypes.c_char_p(b'data_dir'),
                ctypes.c_char_p(force_bytes(data_dir)))
            print('Done setting data_dir to ',data_dir)
        self.settings_called=True
        return
    
    def change_model(self,model=b'twop',data_dir=b'data',verbose=1):
        """
        Change to a different model 
        """
        self.__free_pointers()
        self.__make_pointers(model,data_dir,verbose)
        return

    def compute_point(self,params,verbose=1):
        """
        Compute one point
        """
        
        if self.init_called==False:
            print('Function bamr_init() was not called.')
            quit()
        if self.bamr_class_ptr==0:
            print('Bamr class pointer is 0.')
            quit()
        double_ptr=ctypes.POINTER(ctypes.c_double)
        log_wgt=ctypes.c_double(0.0)
        cp_fun=self.bamr_lib.compute_point
        cp_fun.argtypes=[ctypes.c_void_p,ctypes.c_void_p,
                         ctypes.c_int,double_ptr,double_ptr]
        cp_fun.restype=ctypes.c_int
        point=(ctypes.c_double * self.n_params)()
        for i in range(0,self.n_params):
            if verbose>1:
                print(i,params[i])
            point[i]=params[i]

        if verbose>2:
            print('Going to compute_point().')
        iret=cp_fun(self.bamr_class_ptr,self.model_data_ptr,
                    self.n_params,point,ctypes.byref(log_wgt))
        if verbose>2:
            print('iret,log_wgt',iret,log_wgt.value)
        log_wgt=log_wgt.value
        if iret!=0:
            low_wgt=-800
        if verbose>2:
            print('Done in compute_point().')

        # End of compute_point()
        return log_wgt

    def summarize_tables(self):
        self.bamr_lib.summarize_tables.argtypes=[ctypes.c_void_p]
        print(' ')
        self.bamr_lib.summarize_tables(self.model_data_ptr)
        return

    def get_mvsr_constant(self,name):
        """
        Get a constant from the list associated with the M-R table
        """
        self.bamr_lib.get_mvsr_constant.argtypes=[ctypes.c_void_p,
                                                  ctypes.c_char_p]
        self.bamr_lib.get_mvsr_constant.restype=ctypes.c_double
        con_name=ctypes.c_char_p(force_bytes(name))
        value=self.bamr_lib.get_mvsr_constant(self.model_data_ptr,
                                              con_name)
        return value

    def get_mvsr_column(self,name):
        """
        Setup to call the get_mvsr_column function
        """
        gmc_fun=self.bamr_lib.get_mvsr_column
        nrows=ctypes.c_int(0)
        gmc_fun.argtypes=[ctypes.c_void_p,ctypes.c_char_p,
                          int_ptr,double_ptr_ptr]
        gmc_fun.restype=ctypes.c_int

        # Get the gravitational masses from the table
        col_name=ctypes.c_char_p(force_bytes(name))
        ptr_gm=double_ptr()
        gmc_ret=gmc_fun(model_data_ptr,col_name,ctypes.byref(nrows),
                        ctypes.byref(ptr_gm))
        return gmc_ret

    def get_eos_column(self,name):
        """
        Setup to call the get_eos_column function
        """
        gec_fun=self.bamr_lib.get_eos_column
        nrows=ctypes.c_int(0)
        gec_fun.argtypes=[ctypes.c_void_p,ctypes.c_char_p,
                           int_ptr,double_ptr_ptr]
        gec_fun.restype=ctypes.c_int
        
        # Get the energy densities from the table
        col_name=ctypes.c_char_p(force_bytes(name))
        ptr_ed=double_ptr()
        gec_ret=gec_fun(model_data_ptr,col_name,ctypes.byref(nrows),
                          ctypes.byref(ptr_ed))
        return gec_ret

    def get_grid_column(self,name):
        """
        Setup to call the get_grid_column function
        """
        ggc_fun=self.bamr_lib.get_grid_column
        nrows=ctypes.c_int(0)
        ggc_fun.argtypes=[ctypes.c_void_p,ctypes.c_char_p,
                           int_ptr,double_ptr_ptr]
        ggc_fun.restype=ctypes.c_int
        
        col_name=ctypes.c_char_p(force_bytes(name))
        ptr_I_bar=double_ptr()
        
        ggc_ret=ggc_fun(model_data_ptr,col_name,ctypes.byref(nrows),
                          ctypes.byref(ptr_I_bar))
        return ggc_ret

    def add_data(self,name,fname,fname_alt,slice_name,mass_frac,
                     table):
        """
        Desc
        """
        if self.init_called==True:
            print('Cannot add data after bamr_init() was called.')
            quit()
            
        ad_fun=self.bamr_lib.add_data
        ad_fun.argtypes=[ctypes.c_void_p,ctypes.c_char_p,
                          ctypes.c_char_p,
                          ctypes.c_char_p,ctypes.c_double,
                          ctypes.c_char_p]
        name2=ctypes.c_char_p(force_bytes(name))
        fname2=ctypes.c_char_p(force_bytes(fname))
        slice_name2=ctypes.c_char_p(force_bytes(slice_name))
        table2=ctypes.c_char_p(force_bytes(table))

        if self.apply_intsc:
            self.n_params=self.n_params+2
        else:
            self.n_params=self.n_params+1
        
        ad_fun(self.ns_data_ptr,name2,fname2,slice_name2,
                mass_frac,table2)
        return

    def add_data_alt(self,name,fname,fname_alt,slice_name,mass_frac,
                     table):
        """
        Desc
        """
        if self.init_called==True:
            print('Cannot add data after bamr_init() was called.')
            quit()
            
        ada_fun=self.bamr_lib.add_data_alt
        ada_fun.argtypes=[ctypes.c_void_p,ctypes.c_char_p,
                          ctypes.c_char_p,ctypes.c_char_p,
                          ctypes.c_char_p,ctypes.c_double,
                          ctypes.c_char_p]
        name2=ctypes.c_char_p(force_bytes(name))
        fname2=ctypes.c_char_p(force_bytes(fname))
        fname_alt2=ctypes.c_char_p(force_bytes(fname_alt))
        slice_name2=ctypes.c_char_p(force_bytes(slice_name))
        table2=ctypes.c_char_p(force_bytes(table))

        if self.apply_intsc:
            self.n_params=self.n_params+2
        else:
            self.n_params=self.n_params+1
        
        ada_fun(self.ns_data_ptr,name2,fname2,fname_alt2,slice_name2,
                mass_frac,table2)
        return

    def test_twop(self,verbose=1):
        """
        """
        (iret,lw)=self.compute_point([1.0,-3.0,0.165,0.644,
                                      1.51,0.576,4.6,1.21],
                                     verbose)
        print('log_wgt: %7.6e' % lw)
        print('M_max: %7.6e' % self.get_mvsr_constant('M_max'))
        return
    

