import h5py
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeRegressor 
from sklearn.model_selection import GridSearchCV
plt.rc('font', size=12.0)

dir = '/home/anik/bamr/out/aff_inv/'
model, label = ['mp'], ['MP']
x_name = ['a', 'alpha', 'param_S', 'param_L', '', 'trans1', '',  'trans2', '']

if model[0] == 'mp' or model[0] == 'np':
    x_name[4], x_name[6], x_name[8] = 'exp1', 'exp2', 'exp3'
else:
    x_name[4], x_name[6], x_name[8] = 'csq1', 'csq2', 'csq3'

x_name = x_name + ['M_chirp_det', 'q']
y_name = [f"R_{i}" for i in range(100)]
y_name = y_name + ['M_max', 'dpdM', 'I1', 'I2']

mchain = h5py.File(dir + model[0] + '_all', 'r')['markov_chain_0']
nrow, data = mchain['nlines'][0], mchain['data']
x_ncol, y_ncol = len(x_name), len(y_name)
x, y = np.zeros((x_ncol, nrow)), np.zeros((y_ncol, nrow))
for i in range(x_ncol):
    x[i] = data[x_name[i]]
for i in range(y_ncol):
    y[i] = data[y_name[i]]
x, y = x.T, y.T
x_tr, x_ts, y_tr, y_ts = train_test_split(x, y, test_size=0.1, random_state=42)

DTR=DecisionTreeRegressor(criterion='absolute_error', random_state=42)
scoring='neg_mean_squared_error'
cv, n_jobs=10, -1
hp_dtr = {
    'max_depth': [18, 20, 22],   # NP: 10, MP: 22
    'min_samples_leaf': [2],     # NP: 1,  MP: 2
    #'min_samples_split': [2],   # NP: 2,  MP: 2
}
gs_dtr = GridSearchCV(estimator=DTR, param_grid=hp_dtr,
                       scoring=scoring, cv=cv, n_jobs=n_jobs,
                       return_train_score=True, verbose=2)
gs_dtr.fit(x_tr, y_tr)
best_dtr = gs_dtr.best_estimator_
best_hps = gs_dtr.best_params_
print(f"Best HPs for {label[0]}: {best_hps}")