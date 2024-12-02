import h5py
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score

models = ['np']
labels = ['NL']
x_name = ['a', 'alpha', 'param_S', 'param_L', '', 'trans1', '', 
           'trans2', '']
for k in range(len(models)):
    if models[k] == 'mp' or models[k] == 'np':
        x_name[4], x_name[6], x_name[8] = 'exp1', 'exp2', 'exp3'
    else:
        x_name[4], x_name[6], x_name[8] = 'csq1', 'csq2', 'csq3'
x_name = x_name + ['M_chirp_det', 'q']
y_name = [f"R_{i}" for i in range(100)]
y_name = y_name + ['M_max', 'I1', 'I2']

dir = '/home/anik/bamr/out/aff_inv/'
mchain = h5py.File(dir + models[0] + '_all', 'r')['markov_chain_0']
nrow, data = mchain['nlines'][0], mchain['data']

x_ncol, y_ncol = len(x_name), len(y_name)
x, y = np.zeros((x_ncol, nrow)), np.zeros((y_ncol, nrow))
for i in range(x_ncol):
    x[i] = data[x_name[i]]
for i in range(y_ncol):
    y[i] = data[y_name[i]]
x, y = x.T, y.T

x_tr, x_ts, y_tr, y_ts = train_test_split(x, y, test_size=0.2, random_state=42)

DTR=DecisionTreeRegressor(criterion='absolute_error', random_state=42)
RFR=RandomForestRegressor(criterion='absolute_error', random_state=42)

scoring='neg_mean_squared_error'
cv, n_jobs=5, -1

hp_dtr = {
    'max_depth': [10, 15, 20, None],
    'min_samples_leaf': [1, 2],
    'min_samples_split': [2, 3],
}

gs_dtr = GridSearchCV(estimator=DTR, param_grid=hp_dtr,
                       scoring=scoring, cv=cv, n_jobs=n_jobs,
                       return_train_score=True, verbose=2)

gs_dtr.fit(x_tr, y_tr)

best_dtr = gs_dtr.best_estimator_
best_hps = gs_dtr.best_params_
print(f"Best HPs for DTR-1: {best_hps}")

y_pr = best_dtr.predict(x_ts)

r2 = r2_score(y_ts, y_pr)
print(f"R2 score: {r2}")