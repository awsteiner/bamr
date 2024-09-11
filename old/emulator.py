import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import keras_tuner as kt
from keras_tuner import HyperModel
import tensorflow as tf
from tensorflow import keras
from keras import layers, callbacks
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler



class hyper_search(HyperModel):

  def __init__(self, x_shape, y_shape):
    self.x_shape = x_shape
    self.y_shape = y_shape


  def build(self, hp):
    model = keras.Sequential()
    nx, ny = self.x_shape, self.y_shape
    for i in range(hp.Int('n_layers', min_value=1, max_value=6)):
      model.add(
        layers.Dense(
          units=hp.Int(f"units_{i}", min_value=nx, max_value=10*nx, step=nx),
          activation='relu'
        )
      )
    model.add(layers.Dense(ny, activation='linear'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    return model



class neural_network:

  def __init__(self, filename, x_list, y_list, project):
    self.filename = filename
    self.x_list = x_list
    self.y_list = y_list
    self.project = project
  

  def read_data(self):
    file = h5py.File(self.filename, 'r')
    mchain = file['markov_chain_0']
    nlines = mchain['nlines'][0]
    dtable = mchain['data']
    x = np.zeros((nlines, len(self.x_list)))
    y = np.zeros((nlines, len(self.y_list)))
    for i, params in enumerate(self.x_list):
      x[:,i] = np.array(dtable[self.x_list[i]])
    for i, quants in enumerate(self.y_list):
      y[:,i] = np.array(dtable[self.y_list[i]])
    return x, y


  def scale_data(self, x, y, is_normed=False, x_raw=[], y_raw=[]):
    xy = np.concatenate((x,y), axis=1)
    scaler = MinMaxScaler()
    if (not is_normed):
      xy_scaled = scaler.fit_transform(xy)
    else:
      xy_raw = np.concatenate((x_raw,y_raw), axis=1)
      scaler.fit(xy_raw)
      xy_scaled = scaler.inverse_transform(xy)
    x_scaled = xy_scaled[:,:x.shape[1]]
    y_scaled = xy_scaled[:,x.shape[1]:]
    return x_scaled, y_scaled
    

  def process_data(self):
    x, y = self.read_data()
    xs, ys = self.scale_data(x, y)
    return xs, ys


  def build_model(self, arch, act_func, loss, optimizer):
    model = keras.Sequential()
    model.add(layers.Dense(arch[1], input_shape=(arch[0],), \
                           activation=act_func[0]))
    for i, nn in enumerate(arch):
      if (i>1 and i<len(arch)-1):
        model.add(layers.Dense(arch[i], activation=act_func[0]))
    model.add(layers.Dense(arch[len(arch)-1], activation=act_func[1]))
    model.compile(loss=loss, optimizer=optimizer)
    model.summary()
    return model


  def load_model(self, is_tuned=False):
    if (not is_tuned):
      model_path = './emulator/checkpoints/'
      if not os.path.exists(model_path):
        print("Loading failed: Trained model not found!")
        exit(1)
    else:
      model_path = './emulator/model_'+self.project+'.h5'
      if not os.path.exists(model_path):
        print("Loading failed: Tuned model not found!")
        exit(1)
    model = tf.keras.models.load_model(model_path)
    model.summary()
    return model


  def set_callbacks(self, model_ckpt=False, backup_restore=False, \
                    early_stop=False, var_lr=False, reduce_lr=False, \
                      term_nan=False, verb=0):
    call_backs = []
    
    if (model_ckpt):
      ckpt_path = './emulator/checkpoints/'
      if not os.path.exists(ckpt_path):
        os.makedirs(ckpt_path)
      save_wgts = callbacks.ModelCheckpoint(filepath=ckpt_path, verbose=verb, \
                                            save_best_only=True)
      call_backs.append(save_wgts)

    if (backup_restore):
      backup_dir='./emulator/backup/'
      if not os.path.exists(backup_dir):
        os.makedirs(backup_dir)
      backup_wgts = callbacks.BackupAndRestore(backup_dir, save_freq='epoch', \
                                           delete_checkpoint=True)
      call_backs.append(backup_wgts)

    if (early_stop):
      stop_early = callbacks.EarlyStopping(monitor='loss', min_delta=1.0e-6, \
                                           patience=5, verbose=verb)
      call_backs.append(stop_early)

    if (term_nan):
      term_on_nan = callbacks.TerminateOnNaN()
      call_backs.append(term_on_nan)

    return call_backs
  

  def train(self):
    x, y = self.process_data()
    x_tr, x_tv, y_tr, y_tv = \
      train_test_split(x, y, test_size=0.2, shuffle=True, random_state=42)
    x_ts, x_vl, y_ts, y_vl = \
      train_test_split(x_tv, y_tv, test_size=0.25, shuffle=True, random_state=42)
    call_backs = self.set_callbacks(model_ckpt=True, early_stop=True, \
                                   term_nan=True, verb=0)
    model = self.load_model(is_tuned=True)
    train = model.fit(x=x_tr, y=y_tr, batch_size=512, validation_data=(x_ts,y_ts), \
                    epochs=5000, callbacks=call_backs, verbose=2)
    loss = model.evaluate(x_vl, y_vl, verbose=0)
    print("Loss = {:.4e}".format(loss))
    return train
  

  def predict(self, x):
    model_path = './emulator/checkpoints/'
    model = tf.keras.models.load_model(model_path)
    y = model(x)
    return y
  

  def end_session(self):
    from numba import cuda
    print('Releasing VRAM...')
    cuda.select_device(0)
    cuda.close()
    print('GPU session ended.')


  def emulate(self):
    pass


  def set_tuner(self, hypermodel):
    tuner = kt.GridSearch(
      hypermodel, 
      objective='val_loss', 
      max_trials=1000, 
      seed=42,
      hyperparameters=None,
      tune_new_entries=True,
      allow_new_entries=True,
      max_retries_per_trial=0,
      max_consecutive_failed_trials=3,
      executions_per_trial=1, 
      directory='./emulator/grid_search/',
      project_name=self.project, 
      overwrite=True
    )
    tuner.search_space_summary()
    return tuner


  def search(self):
    x, y = self.process_data()
    x_tr, x_ts, y_tr, y_ts = train_test_split(x, y, test_size=0.2, \
                                              shuffle=True, random_state=42)
    x_shape, y_shape = len(self.x_list), len(self.y_list)
    hypermodel = hyper_search(x_shape, y_shape)
    tuner = self.set_tuner(hypermodel)
    call_backs = self.set_callbacks(early_stop=True)
    tuner.search(x_tr, y_tr, batch_size=512, epochs=1000, validation_data=(x_ts, y_ts), \
                callbacks=call_backs, verbose=2)
    tuner.results_summary(num_trials=5)
    #best_hp = tuner.get_best_hyperparameters()[0]
    #best_model = hypermodel.build(best_hp)
    models = tuner.get_best_models(num_models=1)
    best_model = models[0]
    best_model.build(input_shape=(None, x_shape))
    best_model.summary()
    best_model.save('./emulator/model_'+self.project+'.h5')


if __name__ == '__main__':
  file_name='out/nl_20'
  x_list=['a','alpha','param_S','param_L','csq1','trans1','csq2','trans2','csq3', \
        'M_chirp_det','q','z_cdf','m1_gw19','mf_6304','mf_6397','mf_M13','mf_M28', \
        'mf_M30','mf_wCen','mf_X7','mf_1810b','mf_1724b','mf_1702','mf_0030', \
        'mf_0740','mean_NS','log10_width_NS','skewness_NS','mean_WD','log10_width_WD', \
        'skewness_WD','mean_LMS','log10_width_LMS','skewness_LMS','M_J0453p','M_J0453c', \
        'M_J1906p','M_J1906c','M_B1534p','M_B1534c','M_B1913p','M_B1913c','M_B2127p', \
        'M_B2127c','M_J0737A','M_J0737B','M_J1756p','M_J1756c','M_J1807p','M_J1807c', \
        'M_J1518p','M_J1518c','M_J1811p','M_J1811c','M_J1829p','M_J1829c','M_J2045', \
        'M_J2053','M_J1713','M_B1855','M_J0751','M_J1141','M_J1738','M_J1614','M_J0348', \
        'M_J2222','M_J2234','M_J1949','M_J1012','M_J0437','M_J1909','M_J1802','M_J1911', \
        'M_J2043','M_J0337','M_J1946','M_J1918','M_J1600','M_J0621','M_B2303','M_J0024', \
        'M_J0514','M_B1516','M_J1748I','M_J1748J','M_B1802','M_B1911','M_J0740','M_CygX2', \
        'M_XTEJ2123','M_4U1822','M_HerX1','M_2S0921']
  y_list=['log_wgt']
  project = 'nl'
  print(x_list)
  em = neural_network(file_name, x_list, y_list, project)
  em.search()
