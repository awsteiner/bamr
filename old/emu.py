#!/usr/bin/python3
import math
import h5py
import os
import re
import sys
import getopt
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF

# normalization function
def norm(data):
    d_mean = np.mean(data)
    d_std = np.std(data)
    d_norm = np.array([])
    for i in range(0, len(data)):
        temp = (data[i]-d_mean)/d_std
        d_norm = np.append(d_norm, temp)
    return d_norm, d_mean, d_std

# Inverse normalization
def inv_norm(data, d_mean, d_std):
    inv_data = []
    for i in range(0, len(data)):
        temp = data[i] * d_std + d_mean
        inv_data.append(temp)
    return inv_data

# GPR Class
class modGpr:
    
    def __init__(self):
        
        # Class variable for GPR
        self.gpr = 0

        # Class variable for parameter names
        self.params = 0

        # Class variable for parameter values
        self.init_vals = 0

        # Class variable for target columns
        self.target_cols = ['log_wgt', 'R_43', 'M_max', 'e_max']

        # additional target column for speed of sound
        for i in range(0, 100):
            self.target_cols.append('cs2_'+str(i))

        # Class variable for parameter columns std
        self.param_std_train = {}

        # Class variable for parameter columns mean
        self.param_mean_train = {}

        # Class variable for target columns mean
        self.target_mean_train = {}

        # Class variable for target columns std
        self.target_std_train = {}

        # Check the source file numbers
        self.sources = 0

    # PyObject_CallObject does not update the instance.
    # Issue with calling different methods
    '''
    def show(self, param_vals):
        # copy parameter initial values
        init_vals = param_vals
        print("PyMOdule : Parameter values as arguments : ", init_vals)
        print("PyModule : Trained model values : ", self.mod_vals)
        # print("PyModule : Trained models :", self.gpr)
    '''

    def modTrain(self, hdf_file, param_name, param_vals, n_sources):
        """
        Train the emulator based on file named `hdf_file` using
        parameters from the list `param_name` or evaluate the emulator
        based on the values in `param_vals` presuming the number of
        sources is `n_sources`.
        """
        
        # Update class parameter names
        self.params = param_name

        # Update Parameter values 
        self.init_vals = param_vals

        # Addtinal source check
        self.sources = n_sources

        # Check if training is done once
        if(self.gpr == 0):
            if(self.sources > 0):
                print("Pymodule : adding atmosphere columns")
                for i in range(0, self.sources):
                    temp_atm = "alt_"+str(i)
                    self.params.append(temp_atm)

            print("PyModule : Parameter names :", self.params)
            print("PyModule : Parameter names :", self.target_cols)

            # Read train file
            train_file = h5py.File(hdf_file, 'r')

            # Read the training columns
            X_train = np.array([[]])
            for i in range(0, len(self.params)):
                if(i == 0):
                    temp = np.array(
                        train_file.file["markov_chain_0/data/"+
                                        self.params[i]])
                    X_train = np.array([temp])
                else:
                    temp = np.array(
                        train_file.file["markov_chain_0/data/"+
                                        self.params[i]])
                    X_train = np.vstack([X_train, [temp]])

            # Read the target columns
            Y_train = np.array([[]])
            for i in range(0, len(self.target_cols)):
                if(i == 0):
                    temp = np.array(
                        train_file.file["markov_chain_0/data/"+
                                        self.target_cols[i]])
                    Y_train = np.array([temp])
                else:
                    temp = np.array(
                        train_file.file["markov_chain_0/data/"+
                                        self.target_cols[i]])
                    Y_train = np.vstack([Y_train, [temp]])

            train_file.close()

            # Standardize the parameters
            for i in range(0, len(X_train)):
                temp, temp_mean, temp_std = norm(X_train[i])
                self.param_mean_train.update({self.params[i]: temp_mean})
                self.param_std_train.update({self.params[i]: temp_std})
                X_train[i] = temp

            X_train = X_train.transpose()

            # Standardize the targets
            for i in range(0, len(Y_train)):
                temp, temp_mean, temp_std = norm(Y_train[i])
                self.target_mean_train.update({self.target_cols[i]: temp_mean})
                self.target_std_train.update({self.target_cols[i]: temp_std})
                Y_train[i] = temp
            Y_train = Y_train.transpose()

            print("PyModule : Training array : ", X_train.shape)
            print("PyModule : Target array : ", Y_train.shape)

            print("PyModule : Training GPR model.")
            kernel = 1.0 * RBF(1)
            self.gpr = GPR(kernel=kernel, random_state=0).fit(X_train, Y_train)

            print("PyModule : Training done : ", self.gpr)

            # delete training arrays
            del X_train
            del Y_train
            
            return 0

        # If trained model is available, predictions from given
        # parameter values
        else:
            
            #print("PyModule : GPR model already exist.")
            # normalize initial parameter values
            
            norm_init_vals = np.array([])

            ##################################################################
            # parameter values from bint class does not
            # contain "alt" values.
            # Need to include those values for other models except nodata
            #################################################################

            for i in range(0, len(self.init_vals)):
                temp_mean = self.param_mean_train[self.params[i]]
                temp_std = self.param_std_train[self.params[i]]
                norm_init_vals = np.append(norm_init_vals,
                                           (self.init_vals[i]-temp_mean)/
                                           temp_std)

            norm_init_vals = norm_init_vals.reshape(1, len(self.params))
            #print("PyModule : normalized parameter values: ", norm_init_vals)

            # prediction from given parameter values
            predicted = self.gpr.predict(norm_init_vals, return_std=False)

            re_predicted = []
            for i in range(0, len(self.target_cols)):
                re_temp = (predicted[0][i] *
                           self.target_std_train[self.target_cols[i]] +
                           self.target_mean_train[self.target_cols[i]])
                re_predicted.append(re_temp)
                
                #print("Pymodule : predicted M_max :", re_predicted[2])
                
            predicted_log_wgt = (predicted[0][0] *
                                 self.target_std_train['log_wgt'] +
                                 self.target_mean_train['log_wgt'])
            #print(predicted)

            return re_predicted
