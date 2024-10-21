import h5py
import numpy as np
import o2sclpy as o2py
from sklearn.model_selection import train_test_split

def read_table(tab, list):
    data = []
    for i in list:
        data.append(tab[i][0:tab.get_nlines()])
    data = np.array(data).T
    return data

def get_mvsr(eosp, fname):
    radii = [f"R_{i}" for i in range(100)]
    o2set=o2py.lib_settings_class()
    
    hf = o2py.hdf_file()
    hf.open(fname)
    tab = o2py.table()
    name = b''
    o2py.hdf_input_table(hf, tab, name)
    hf.close()

    x = read_table(tab, eosp)
    y = read_table(tab, radii)
    x_tr, x_ts, y_tr, y_ts = train_test_split(x, y, test_size=0.2, random_state=42)
    model = o2py.interpm_sklearn_dtr()
    model.set_data(x_tr, y_tr)
    return model
    