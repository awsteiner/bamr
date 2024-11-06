import numpy
from o2sclpy.utils import string_to_dict2

class interpm_tf_dnn:
    """
    Interpolate one or many multimensional data sets using 
    a deep neural network from TensorFlow
    """

    def __init__(self):
        
        self.verbose=0
        self.kernel=0
        self.outformat='numpy'
        self.dnn=0
        self.SS1=0
        self.SS2=0
        self.transform_in=0
        self.transform_out=0
        
        return
    
    def set_data(self,in_data,out_data,outformat='numpy',verbose=0,
                 transform_in='none',transform_out='none',
                 test_size=0.0,evaluate=False,
                 hlayers=[8,8],activations=['relu','relu'],
                 out_act='linear',batch_size=32,epochs=100,
                 loss='mean_squared_error',monitor='val_loss',
                 ls_factor=0.5,ls_patience=10,ls_min_lr=1.0e-6,
                 es_min_delta=1.0e-4,es_patience=100):
        """
        Set the input and output data to train the interpolator

        some activation functions: 
        relu [0,\\infty]
        sigmoid [0,1]
        tanh [-1,1]

        transformations:
        quantile transforms to [0,1]
        MinMaxScaler transforms to [a,b]
        """

        from sklearn.model_selection import train_test_split
        import tensorflow as tf
        
        if verbose>0:
            print('interpm_tf_dnn::set_data():')
            print('  outformat:',outformat)
            print('  in_data shape:',numpy.shape(in_data))
            print('  out_data shape:',numpy.shape(out_data))
            print('  batch_size:',batch_size)
            print('  layers:',hlayers)
            print('  activation functions:',activations)
            print('  transform_in:',transform_in)
            print('  transform_out:',transform_out)
            print('  epochs:',epochs)
            print('  test_size:',test_size)

        self.outformat=outformat
        self.verbose=verbose
        self.transform_in=transform_in
        self.transform_out=transform_out

        # ----------------------------------------------------------
        # Handle the data transformations
        
        from sklearn.preprocessing import QuantileTransformer
        from sklearn.preprocessing import MinMaxScaler
        from sklearn.preprocessing import StandardScaler
        
        if self.transform_in=='minmax_0':
            self.SS1=MinMaxScaler(feature_range=(0,1))
            i_nonzeros=in_data!=0
            in_data_trans=numpy.copy(in_data)
            in_data_trans[i_nonzeros]=self.SS1.fit_transform(
                in_data[i_nonzeros].reshape(-1,1)).flatten()
        elif self.transform_in=='minmax_1':
            self.SS1=MinMaxScaler(feature_range=(0,1))        
            in_data_trans=self.SS1.fit_transform(in_data)
        elif self.transform_in=='minmax_pm':
            self.SS1=MinMaxScaler(feature_range=(-1,1))
            in_data_trans=self.SS1.fit_transform(in_data)
        elif self.transform_in=='quant':
            self.SS1=QuantileTransformer(n_quantiles=in_data.shape[0])
            in_data_trans=self.SS1.fit_transform(in_data)
        elif self.transform_in=='standard':
            self.SS1=StandardScaler()
            in_data_trans=self.SS1.fit_transform(in_data)
        else:
            in_data_trans=in_data
        
        if self.transform_out=='minmax_0':
            try:
                i_nonzeros=out_data!=0
                out_data_trans=numpy.copy(out_data)
                out_data_trans[i_nonzeros]=self.SS1.fit_transform(
                    out_data[i_nonzeros].reshape(-1,1)).flatten()
            except:
                out_data=numpy.array(out_data).reshape(-1,1)
                out_data_trans=self.SS1.transform(out_data)
        elif self.transform_out=='minmax_1':
            out_data_trans=self.SS1.fit_transform(out_data)
        elif self.transform_out=='minmax_pm':
            self.SS2=MinMaxScaler(feature_range=(-1,1))
            out_data_trans=self.SS2.fit_transform(out_data)
        elif self.transform_out=='quant':
            self.SS2=QuantileTransformer(n_quantiles=out_data.shape[0])
            out_data_trans=self.SS2.fit_transform(out_data)
        elif self.transform_out=='standard':
            self.SS2=StandardScaler()
            out_data_trans=self.SS2.fit_transform(out_data)
        else:
            out_data_trans=out_data

        if self.verbose>0:
            try:
                minv=out_data_trans[0,0]
                maxv=out_data_trans[0,0]
                minv_old=out_data[0,0]
                maxv_old=out_data[0,0]
            except Exception as e:
                print('Exception in interpm_tf_dnn::set_data()',
                      'at min,max().',e)
            
            for j in range(0,numpy.shape(out_data)[0]):
                if out_data[j,0]<minv_old:
                    minv_old=out_data[j,0]
                if out_data[j,0]>maxv_old:
                    maxv_old=out_data[j,0]
                if out_data_trans[j,0]<minv:
                    minv=out_data_trans[j,0]
                if out_data_trans[j,0]>maxv:
                    maxv=out_data_trans[j,0]

            print('min,max before transformation: %7.6e %7.6e' %
                  (minv_old,maxv_old))
            print('min,max after transformation : %7.6e %7.6e' %
                  (minv,maxv))
            
        if test_size>0.0:
            try:
                x_train,x_test,y_train,y_test=train_test_split(
                    in_data_trans,out_data_trans,test_size=test_size)
            except Exception as e:
                print('Exception in interpm_tf_dnn::set_data()',
                      'at test_train_split().',e)
        else:
            x_train=in_data_trans
            y_train=out_data_trans

        nd_in=numpy.shape(in_data)[1]
        nd_out=numpy.shape(out_data)[1]
        
        if self.verbose>0:
            print('nd_in,nd_out:',nd_in,nd_out)
            print('  Training DNN model.')
            
        try:
            nl=len(hlayers)
            na=len(activations)
            inp=tf.keras.Input(shape=(nd_in,))
            layers=[inp,tf.keras.layers.Dense(hlayers[0],
                                              activation=activations[0])]
            if self.verbose>0:
                print('Layer: dense',hlayers[0],nd_in,activations[0])
            if nl>0:
                for i in range(1,nl):
                    act=activations[i%na]
                    layers.append(tf.keras.layers.Dense(hlayers[i],
                                                        activation=act))
                    if self.verbose>0:
                        print('Layer: dense',hlayers[i],act)
            layers.append(tf.keras.layers.Dense(nd_out,
                                                activation=out_act))
            if self.verbose>0:
                print('Layer: dense',nd_out,out_act)
            model=tf.keras.Sequential(layers)
        except Exception as e:
            print('Exception in interpm_tf_dnn::set_data()',
                  'at model definition.',e)
        
        if self.verbose>0:
            print('summary:',model.summary())

        try:
            import keras
            from keras.callbacks import EarlyStopping # type: ignore
            from keras.callbacks import ReduceLROnPlateau # type: ignore

            class loss_history(keras.callbacks.Callback):
                def on_train_begin(self,logs={}):
                    self.losses=[]
                def on_batch_end(self,batch,logs={}):
                    self.losses.append(logs.get('loss'))
            history=loss_history()

            early_stopping=EarlyStopping(monitor=monitor,
                                         min_delta=es_min_delta,
                                         patience=es_patience,
                                         restore_best_weights=True,
                                         mode='min')
            lr_schedule=ReduceLROnPlateau(monitor=monitor,
                                          factor=ls_factor,
                                          patience=ls_patience,
                                          min_lr=ls_min_lr)

            model.compile(loss=loss,optimizer='adam')
                
            if test_size>0.0:
                # Fit the model to training data
                model.fit(x_train,y_train,batch_size=batch_size,
                          epochs=epochs,validation_data=(x_test,y_test),
                          verbose=self.verbose,
                          callbacks=[history,early_stopping,lr_schedule])
                          
            else:
                # Fit the model to training data
                model.fit(x_train,y_train,batch_size=batch_size,
                          epochs=epochs,verbose=self.verbose,
                          callbacks=[history,early_stopping,lr_schedule])
                
        except Exception as e:
            print('Exception in interpm_tf_dnn::set_data()',
                  'at model fitting.',e)
            raise

        if evaluate==True:
            # Return loss value and metrics
            if self.verbose>0:
                print('  Training done.')
            print('  Test Score: [loss, accuracy]:',
                  model.evaluate(x_test,y_test,verbose=self.verbose))
            
        self.dnn=model

        return
    
    def set_data_str(self,in_data,out_data,options):
        """
        Set the input and output data to train the interpolator,
        using a string to specify the keyword arguments.
        """

        try:
            dct=string_to_dict2(options,list_of_ints=['verbose',
                                                      'batch_size',
                                                      'epochs'],
                                list_of_floats=['test_size'],
                                list_of_bools=['evaluate'])
            if "hlayers" in dct:
                htemp=dct["hlayers"]
                htemp=htemp[1:-1]
                htemp=htemp.split(',')
                htemp2=[]
                for i in range(0,len(htemp)):
                    htemp2.append(int(htemp[i]))
                dct["hlayers"]=htemp2
            if "activations" in dct:
                atemp=dct["activations"]
                atemp=atemp[1:-1]
                atemp=atemp.split(',')
                atemp2=[]
                for i in range(0,len(atemp)):
                    atemp2.append(atemp[i])
                dct["activations"]=atemp2
            print('String:',options,'Dictionary:',dct)

            self.set_data(in_data,out_data,**dct)
        except Exception as e:
            print('Calls in interpm_tf_dnn::set_data_str() failed.',e)

        return
    
    def eval(self,v):
        """
        Evaluate the NN at point ``v``.
        """

        if self.transform_in!='none':
            v_trans=0
            try:
                v_trans=self.SS1.transform(v.reshape(1,-1))
            except Exception as e:
                print('Exception at input transformation in interpm_tf_dnn:',
                      e)
        else:
            v_trans=v.reshape(1,-1)

        try:
            # We don't want output at every point, even if verbose is
            # 1, so we use self.verbose-1 here for the argument to
            # the predict function.
            if self.verbose>1:
                pred=self.dnn.predict(v_trans, verbose=self.verbose-1)
            else:
                pred=self.dnn.predict(v_trans,verbose=0)
        except Exception as e:
            print('Exception 4 in interpm_tf_dnn:',e)
            
        if self.transform_out!='none':
            try:
                pred_trans=self.SS2.inverse_transform(pred)
            except Exception as e:
                print('Exception 5 in interpm_tf_dnn:',e)
        elif self.transform_in=='minmax_0':
            try:
                i_nonzeros=pred!=0
                pred_trans=numpy.copy(pred)
                pred_trans[i_nonzeros]=self.SS1.inverse_transform(
                    pred[i_nonzeros].reshape(-1,1)).flatten()
            except Exception as e:
                print('Exception 6 in interpm_tf_dnn:',e)
        else:
            pred_trans=pred
    
        if self.outformat=='list':
            return pred_trans.tolist()

        if pred_trans.ndim==1:
            
            if self.verbose>1:
                print('interpm_tf_dnn::eval():',
                      'type(pred_trans),pred_trans:',
                      type(pred_trans),pred_trans)
            # The output from tf.keras is float32, so we have to convert to
            # float64 
            n_out=numpy.shape(pred_trans[0])[0]
            out_double=numpy.zeros((n_out))
            for i in range(0,n_out):
                out_double[i]=pred_trans[i]
                    
            return numpy.ascontiguousarray(out_double)
        
        if self.verbose>1:
            print('interpm_tf_dnn::eval():',
                  'type(pred_trans[0]),pred_trans[0]:',
                  type(pred_trans[0]),pred_trans[0])

        # The output from tf.keras is float32, so we have to convert to
        # float64 
        n_out=numpy.shape(pred_trans[0])[0]
        out_double=numpy.zeros((n_out))
        for i in range(0,n_out):
            out_double[i]=pred_trans[0][i]
            
        return numpy.ascontiguousarray(out_double)