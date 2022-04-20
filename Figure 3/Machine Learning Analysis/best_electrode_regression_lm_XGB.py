# -*- coding: utf-8 -*-
"""
Created on Mon May 18 17:55:08 2020

@author: Pilin
"""

#%%
import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import train_test_split
import itertools
import mne
from mne.decoding import CSP
from mne import Epochs
from mne.decoding import SPoC
mne.set_log_level(verbose='warning') #to avoid info at terminal
import pickle
import sys
import IO
import os
import multiprocessing

#import tensorflow
#import tensorflow as tf
#import keras
#from keras.layers import BatchNormalization
#from keras.models import Sequential
#from keras.layers import Dense, Dropout
#from keras.wrappers.scikit_learn import KerasClassifier, KerasRegressor
#from keras.optimizers import Adam

import tensorflow
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dense, Activation, Permute, Dropout
from tensorflow.keras.layers import Conv2D, MaxPooling2D, AveragePooling2D
from tensorflow.keras.layers import SeparableConv2D, DepthwiseConv2D
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.layers import SpatialDropout2D
from tensorflow.keras.regularizers import l1_l2
from tensorflow.keras.layers import Input, Flatten
from tensorflow.keras.constraints import max_norm
from tensorflow.keras import backend as K
from tensorflow.keras.optimizers import Adam
import tensorflow.keras as keras
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras import backend as K
from tensorflow.keras.models import load_model

#from tensorflow.python.keras import backend as K
from sklearn.model_selection import StratifiedKFold

from sklearn.linear_model import Ridge
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LassoCV
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import ElasticNet
#from keras.callbacks import EarlyStopping
#from keras.callbacks import ModelCheckpoint
#from keras.models import load_model
#from tensorflow.keras.layers import Dense, Activation, Permute, Dropout

from scipy import stats
from collections import OrderedDict
from sklearn.model_selection import KFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn import metrics
from bayes_opt import BayesianOptimization
from skopt.space import Real, Integer, Categorical
from skopt.utils import use_named_args
from skopt import gp_minimize
import gc
from sklearn.preprocessing import StandardScaler

import xgboost as xgb
from xgboost import XGBRegressor

VICTORIA = False
WRITE_OUT_CH_IND = False
USED_MODEL = 2 # 0 - Enet, 1 - XGB, 2 - NN
settings = {}
VERBOSE_ALL = 0

if VICTORIA is True:
    # insert at 1, 0 is the script path (or '' in REPL)
    sys.path.insert(1, '/home/victoria/icn/icn_m1')
    settings['BIDS_path'] = "//mnt/Datos/BML_CNCRS/Data_BIDS_new/"
    settings['out_path'] = "/mnt/Datos/BML_CNCRS/Data_processed/Derivatives/Raw_pipeline/"
    if USED_MODEL==0 :
           settings['out_path_process'] = "/mnt/Datos/BML_CNCRS/Spoc/ECoG_STN/LM_Out/"
    if USED_MODEL==1 :
           settings['out_path_process'] = "/mnt/Datos/BML_CNCRS/Spoc/ECoG_STN/XGB_Out/"
else:
    settings['BIDS_path'] = "C:\\Users\\ICN_admin\\Dropbox (Brain Modulation Lab)\\Shared Lab Folders\\CRCNS\\MOVEMENT DATA\\"
    settings['out_path'] = "C:\\Users\\ICN_admin\\Dropbox (Brain Modulation Lab)\\Shared Lab Folders\\CRCNS\\MOVEMENT DATA\\derivatives\\Int_old_grid\\"
    settings['out_path_process'] = "C:\\Users\\ICN_admin\\Dropbox (Brain Modulation Lab)\\Shared Lab Folders\\CRCNS\MOVEMENT DATA\\ECoG_STN\\NN_Out\\"

settings['frequencyranges']=[[4, 8], [8, 12], [13, 20], [20, 35], [13, 35], [60, 80], [90, 200], [60, 200]]
settings['seglengths']=[1, 2, 2, 3, 3, 3, 10, 10, 10]
settings['num_patients']=['000', '001', '004', '005', '006', '007', '008', '009', '010', '013', '014']
settings['BIDS_path']=settings['BIDS_path'].replace("\\", "/")
settings['out_path']=settings['out_path'].replace("\\", "/")



#%%
space_LM = [Real(0, 1, "uniform", name='alpha'),
           Real(0, 1, "uniform", name='l1_ratio')]

space_XGB  = [Integer(1, 100, name='max_depth'),
          Real(10**-5, 10**0, "log-uniform", name='learning_rate'),
          Real(10**0, 10**1, "uniform", name="gamma")]

space_NN = [Real(low=1e-4, high=1e-2, prior='log-uniform', name='learning_rate'),
              Integer(low=1, high=3, name='num_dense_layers'),
              Integer(low=1, high=10, prior='uniform', name='num_input_nodes'),
              Integer(low=1, high=10, name='num_dense_nodes'),
              Categorical(categories=['relu', 'tanh'], name='activation')]

def optimize_enet(x,y):

    @use_named_args(space_LM)
    def objective(**params):
        reg=ElasticNet(max_iter=1000, normalize=False)
        reg.set_params(**params)
        cval = cross_val_score(reg, x, y, scoring='r2', cv=3)
        cval[np.where(cval < 0)[0]] = 0

        return -cval.mean()

    res_gp = gp_minimize(objective, space_LM, n_calls=20, random_state=0)
    return res_gp

def create_model_NN():
        """
        Create NN tensorflow with different numbers of hidden layers / hidden units
        """

        #start the model making process and create our first layer
        model = tensorflow.keras.Sequential()
        model.add(Dense(num_input_nodes, input_shape=(40,), activation=activation))

        #create a loop making a new dense layer for the amount passed to this model.
        #naming the layers helps avoid tensorflow error deep in the stack trace.
        for i in range(num_dense_layers):
            name = 'layer_dense_{0}'.format(i+1)
            model.add(BatchNormalization())
            model.add(Dropout(0.2))
            model.add(Dense(num_dense_nodes,
                     activation=activation,
                            name=name
                     ))
        #add our classification layer.
        model.add(BatchNormalization())
        model.add(Dropout(0.2))
        model.add(Dense(1,activation='linear'))

        #setup our optimizer and compile
        adam = Adam(lr=learning_rate)
        model.compile(optimizer=adam, loss='mean_squared_error',
                     metrics=['mse'])
        return model

def optimize_nn(x,y):


    @use_named_args(space_NN)
    def objective(**params):
        print(params)
        global learning_rate
        learning_rate=params["learning_rate"]
        global num_dense_layers
        num_dense_layers=params["num_dense_layers"]
        global num_input_nodes
        num_input_nodes=params["num_input_nodes"]
        global num_dense_nodes
        num_dense_nodes=params["num_dense_nodes"]
        global activation
        activation=params["activation"]

        cv = KFold(n_splits=3, shuffle=True)
        cv_res = []
        with tf.device(tf.DeviceSpec(device_type="")):
            for train_index, test_index in cv.split(x):
                X_train, X_test=x[train_index, :], x[test_index, :]
                y_train, y_test=y[train_index], y[test_index]
                X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, train_size=0.8,shuffle=False)
                #model = KerasRegressor(build_fn=create_model_NN, epochs=1000, batch_size=500, verbose=2)
                es = EarlyStopping(monitor='val_mse', mode='min', verbose=VERBOSE_ALL, patience=10)
                mc = ModelCheckpoint('best_model.h5', monitor='val_mse', mode='min', verbose=VERBOSE_ALL, save_best_only=True)
                model = create_model_NN()
                model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=1000, batch_size=500, verbose=VERBOSE_ALL, callbacks=[mc,es])
                model = load_model('best_model.h5')
                sc = metrics.r2_score(model.predict(X_test), y_test)
                if sc < 0: sc = 0
                cv_res.append(sc)

        #cv_res = cross_val_score(model, x, y, cv=3, n_jobs=59, scoring="r2")
        #cv_res[np.where(cv_res < 0)[0]] = 0
        return -np.mean(cv_res)

    res_gp = gp_minimize(objective, space_NN, n_calls=20, random_state=0)
    return res_gp

def optimize_xgb(x,y):

    def evalerror(preds, dtrain):
        """
        Custom defined r^2 evaluation function
        """
        labels = dtrain.get_label()
        # return a pair metric_name, result. The metric name must not contain a
        # colon (:) or a space since preds are margin(before logistic
        # transformation, cutoff at 0)

        r2 = metrics.r2_score(labels, preds)

        if r2 < 0:
            r2 = 0

        return 'r2', r2

    @use_named_args(space_XGB)
    def objective(**params):
        print(params)

        params_ = {'max_depth': int(params["max_depth"]),
             'gamma': params['gamma'],
             #'n_estimators': int(params["n_estimators"]),
             'learning_rate': params["learning_rate"],
             'subsample': 0.8,
             'eta': 0.1,
             'disable_default_eval_metric' : 1}
             # 'scale_pos_weight ' : 1}
             #'nthread':59}
             #'tree_method' : 'gpu_hist'}
             #'gpu_id' : 1}

        cv_result = xgb.cv(params_, xgb.DMatrix(x, label=y), num_boost_round=30, feval=evalerror, nfold=3)
        return -cv_result['test-r2-mean'].iloc[-1]

    res_gp = gp_minimize(objective, space_XGB, n_calls=20, random_state=0)
    return res_gp


def get_int_runs(subject_id, subfolder):
    """

    :param patient_idx:
    :return: list with all run files for the given patient
    """
    os.listdir(settings['out_path'])

    if 'right' in str(subfolder):
        list_subject = [i for i in os.listdir(settings['out_path']) if i.startswith('sub_'+subject_id+'_sess_right') and i.endswith('.p')]
    else:
        list_subject = [i for i in os.listdir(settings['out_path']) if i.startswith('sub_'+subject_id+'_sess_left') and i.endswith('.p')]

    return list_subject


def append_time_dim(arr, y_, time_stamps):
    """
    apply added time dimension for the data array and label given time_stamps (with downsample_rate=100) in 100ms / need to check with 1375Hz
    """
    time_arr = np.zeros([arr.shape[0]-time_stamps, int(time_stamps*arr.shape[1])])
    for time_idx, time_ in enumerate(np.arange(time_stamps, arr.shape[0])):
        for time_point in range(time_stamps):
            time_arr[time_idx, time_point*arr.shape[1]:(time_point+1)*arr.shape[1]] = arr[time_-time_point,:]
    return time_arr, y_[time_stamps:]
#%%
cv = KFold(n_splits=3, shuffle=False)
laterality=[("CON"), ("IPS")]
signal=["STN", "ECOG"]

#%%cross-val within subject
#for sub_idx in np.arange(0, len(settings['num_patients']), 1):

def run_patient(sub_idx):
    print("IN there")

    for signal_idx, signal_ in enumerate(signal):
        subject_path=settings['BIDS_path'] + 'sub-' + settings['num_patients'][sub_idx]
        subfolder=IO.get_subfolders(subject_path)

        for sess_idx in range(len(subfolder)):
            if os.path.exists(os.path.join(settings['out_path_process'],
                        settings['num_patients'][sub_idx]+'BestChpredictions_'+signal_+'-'+
                            str(subfolder[sess_idx])+'.npy')) is True:
                continue
            X=[]
            Y_con=[]
            Y_ips=[]
            list_subject=get_int_runs(settings['num_patients'][sub_idx], subfolder[sess_idx])
            list_subject=sorted(list_subject)
            if signal_=="ECOG":
                if sub_idx==4 and sess_idx==0: #for sake of comparison with spoc
                    list_subject.pop(0)
                if sub_idx==4 and sess_idx==1:
                    list_subject.pop(2)

            print('RUNNIN SUBJECT_'+ settings['num_patients'][sub_idx]+ '_SESS_'+ str(subfolder[sess_idx]) + '_SIGNAL_' + signal_)

            for run_idx in range(len(list_subject)):
                with open(settings['out_path']+ '/'+ list_subject[run_idx], 'rb') as handle:
                    run_ = pickle.load(handle)

                #concatenate features
                #get cortex data only
                if signal_=="ECOG":
                    ind_cortex=run_['used_channels']['cortex']
                    rf=run_['rf_data_median']
                    x=rf[:,ind_cortex,:]
                    x=np.clip(x, -2,2) # this should have been implemented in the pipeline
                    y=run_['label_baseline_corrected']
                    con_true=run_['label_con_true']
                    y_con=np.squeeze(y[con_true==True])
                    y_ips=np.squeeze(y[con_true==False])
                    X.append(x)
                    Y_con.append(y_con)
                    Y_ips.append(y_ips)
                else:
                    ind_subcortex=run_['used_channels']['subcortex']
                    if ind_subcortex is not None:

                        rf=run_['rf_data_median']
                        x=rf[:,ind_subcortex,:]
                        x=np.clip(x, -2,2)

                        y=run_['label_baseline_corrected']
                        con_true=run_['label_con_true']
                        y_con=np.squeeze(y[con_true==True])
                        y_ips=np.squeeze(y[con_true==False])

                        X.append(x)
                        Y_con.append(y_con)
                        Y_ips.append(y_ips)

            gc.collect() # free unreferenced memory

            X=np.concatenate(X, axis=0)
            Y_con=np.concatenate(Y_con, axis=0)
            Y_ips=np.concatenate(Y_ips, axis=0)

            Yp_tr= OrderedDict() # Y_predict_train
            sc_tr= OrderedDict() # score_train
            Yp_te= OrderedDict()
            sc_te= OrderedDict()
            Yt_tr= OrderedDict()
            Yt_te= OrderedDict()
            Model_coef= OrderedDict()
            Model_hyperarams= OrderedDict()

            for laterality_idx, laterality_ in enumerate(laterality):
                print("training %s" %laterality_)
                sc_tr[laterality_] = []
                sc_te[laterality_] = []
                Yp_tr[laterality_] = []
                Yp_te[laterality_] = []
                Yt_tr[laterality_] = []
                Yt_te[laterality_] = []
                Model_coef[laterality_] = []
                Model_hyperarams[laterality_] = []


                if laterality_=="CON":
                    label=Y_con
                else:
                    label=Y_ips

                Score_tr=np.empty(X.shape[1], dtype=object)
                Score_te=np.empty(X.shape[1], dtype=object)
                Label_te=np.empty(X.shape[1], dtype=object)
                Label_tr=np.empty(X.shape[1], dtype=object)
                Labelpre_te=np.empty(X.shape[1], dtype=object)
                Labelpre_tr=np.empty(X.shape[1], dtype=object)
                COEF_ = np.empty(X.shape[1],dtype=object)
                Hyperarapms= np.empty(X.shape[1],dtype=object)


                #for each electrode
                for ch_idx in range(X.shape[1]):
                    print("running channel "+str(ch_idx))
                    Ypre_te= []
                    Ypre_tr= []
                    score_tr= []
                    Ypre_te= []
                    score_te= []
                    label_test=[]
                    label_train=[]
                    coords = []
                    coef_ = []
                    hyp_=[]
                    for train_index, test_index in cv.split(X):
                        Xtr, Xte=X[train_index,ch_idx,:], X[test_index,ch_idx,:]
                        Ytr, Yte=label[train_index], label[test_index]
                        label_test.append(Yte)
                        label_train.append(Ytr)
                        dat_tr,label_tr = append_time_dim(Xtr, Ytr, time_stamps=5)
                        dat_te,label_te = append_time_dim(Xte, Yte, time_stamps=5)

                        if USED_MODEL == 0: # Enet
                            optimizer=optimize_enet(x=dat_tr,y=label_tr)
                            model=ElasticNet(alpha=optimizer['x'][0],
                                               l1_ratio=optimizer['x'][1],
                                               max_iter=1000,
                                               normalize=False)
                        elif USED_MODEL == 1: # XGB
                            optimizer=optimize_xgb(x=dat_tr, y=label_tr)
                            model=XGBRegressor(max_depth=optimizer['x'][0],
                                               learning_rate=optimizer['x'][1],
                                               gamma=optimizer['x'][2], subsample= 0.8,
                                               eta= 0.1,
                                               disable_default_eval_metric= 1)

                        elif USED_MODEL == 2:
                            optimizer=optimize_nn(x=dat_tr, y=label_tr)
                            global learning_rate
                            learning_rate=optimizer['x'][0]
                            global num_dense_layers
                            num_dense_layers=optimizer['x'][1]
                            global num_input_nodes
                            num_input_nodes=optimizer['x'][2]
                            global num_dense_nodes
                            num_dense_nodes=optimizer['x'][3]
                            global activation
                            activation=optimizer['x'][4]

                            model = create_model_NN()
                        else:
                            break
                            print("ARCHITECTURE IS NOT DEFINED")

                        if USED_MODEL == 2:

                            es = EarlyStopping(monitor='val_mse', mode='min', verbose=VERBOSE_ALL, patience=10)
                            mc = ModelCheckpoint('best_model.h5', monitor='val_mse', mode='min', verbose=VERBOSE_ALL, save_best_only=True)
                            X_train, X_val, y_train, y_val = train_test_split(dat_tr, label_tr, train_size=0.8,shuffle=True)
                            model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=1000, batch_size=500, verbose=VERBOSE_ALL, callbacks=[mc,es])
                            r2_tr = metrics.r2_score(model.predict(X_train), y_train)
                            if r2_tr < 0: r2_tr = 0
                            r2_te = metrics.r2_score(model.predict(dat_te), label_te)
                            print("channel: "+str(ch_idx)+" r2 test: "+str(r2_te))
                            if r2_te < 0: r2_te = 0

                        else:
                            model.fit(dat_tr, label_tr)
                            r2_tr=model.score(dat_tr, label_tr)
                            if r2_tr < 0: r2_tr = 0
                            r2_te=model.score(dat_te, label_te)
                            if r2_te < 0: r2_te = 0

                        score_tr.append(r2_tr)
                        score_te.append(r2_te)
                        #with tf.device(tf.DeviceSpec(device_type="CPU", device_index=0)):
                        Ypre_te.append(model.predict(dat_te) if USED_MODEL != 2 else model.predict(dat_te)[:,0])
                        Ypre_tr.append(model.predict(dat_tr) if USED_MODEL != 2 else model.predict(dat_tr)[:,0])
                        if USED_MODEL == 0: coef_.append(model.coef_)
                        hyp_.append(optimizer['x'])

                    Score_tr[ch_idx]=np.mean(score_tr)
                    Score_te[ch_idx]=np.mean(score_te)
                    Label_te[ch_idx]=label_test
                    Label_tr[ch_idx]=label_train
                    Labelpre_te[ch_idx]=Ypre_te
                    Labelpre_tr[ch_idx]=Ypre_tr
                    COEF_[ch_idx]=coef_
                    Hyperarapms[ch_idx]=hyp_

                sc_tr[laterality_] = Score_tr
                sc_te[laterality_] = Score_te
                Yp_tr[laterality_] = Labelpre_te
                Yp_te[laterality_] = Labelpre_tr
                Yt_tr[laterality_] = Label_te
                Yt_te[laterality_] = Label_tr
                Model_coef[laterality_] = COEF_
                Model_hyperarams[laterality_]=Hyperarapms

            predict_ = {
                "y_pred_test": Yp_te,
                "y_test": Yt_te,
                "y_pred_train": Yp_tr,
                "y_train": Yt_tr,
                "score_tr": sc_tr,
                "score_te": sc_te,
                "coord_patient" : run_["coord_patient"],
                "coef" :Model_coef,
                "model_hyperparams": Model_hyperarams


            }
            out_path_file = os.path.join(settings['out_path_process']+ \
                settings['num_patients'][sub_idx]+'BestChpredictions_'+signal_+'-'+ str(subfolder[sess_idx])+'.npy')
            np.save(out_path_file, predict_)

if __name__ == '__main__':
    #for sub_idx in np.arange(0, len(settings['num_patients']), 1):
    #my_devices = tf.config.experimental.list_physical_devices(device_type='CPU')
    #tf.config.experimental.set_visible_devices(devices= my_devices, device_type='CPU')
    #pool = multiprocessing.Pool()
    #pool.map(run_patient, np.arange(0, len(settings['num_patients']), 1))
    sess = tf.Session(config=tf.ConfigProto(
      allow_soft_placement=True, log_device_placement=True))
    for sub in np.arange(0, len(settings['num_patients']), 1):
        run_patient(sub)
