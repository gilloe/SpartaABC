# -*- coding: utf-8 -*-
"""
Created on Tue March  3 17:34:38 2020

@author: gillo
"""

#TODO: ridge from feature importance, CV score, V
#TODO: explain std that it is on nn, mean vs ridge V
#TODO: add click, check if file size is large enough
#TODO: based on agg_results_v2, but needs a polish V
#TODO: copy MSF functions to avoid problems V

# pwd = "D:/university/projects/abc/code/abc_nn/abc_nn/pipeline/"
import sys
# sys.path.insert(1, pwd)
import logging
# logging.basicConfig(filename=pwd+'log/infer_abc_params_single_folder.log',level=logging.INFO,
#                     format='%(asctime)s %(levelname)-8s %(message)s',
#                     datefmt='%Y-%m-%d %H:%M:%S') 
#logging.warning('start')

# import model_selection_functions as msf
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
import os
# import time
#from multiprocessing.dummy import Pool as ThreadPool 
#from multiprocessing import Pool 
from scipy.stats import pearsonr
#from scipy.special import zeta
from sklearn import linear_model 
from sklearn import model_selection
# from sklearn import metrics
import matplotlib.pyplot as plt

from sklearn.utils.multiclass import unique_labels
from keras.layers import Input, Dense
from keras.models import Model
#from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.cross_decomposition import PLSRegression

import random
import click

#%% MSF functions




### Load lib data
def load_lib_data(path,lib='data_name',models_list=['ideq','iddif'] ,size_th=1E6,rel_path=''):
# Written for ideq and iddfif
#models_list = ['ideq','iddif']
#lib = 'EMGT00050000000002' #EMGT00050000000151
#i=0
#file_name = lib + '/SpartaABC_EMGT00050000000151_model_' + models_list[i] + '.posterior_params'
    for i in range(len(models_list)):
        if rel_path == '':
            file_name = path + 'SpartaABC_'+ lib + '_' + models_list[i] + '.posterior_params'
        else:
            file_name = path + rel_path + lib + '/SpartaABC_'+ lib + '_' + models_list[i] + '.posterior_params'
        if os.stat(file_name).st_size<size_th:
            logging.warning(f'{file_name} too small.')
            return None,None
        with open(file_name) as f:
            file_rows_num = sum(1 for line in f)
        df_tmp = pd.read_csv(file_name, delimiter='\t',
                             skiprows=[i for i in range(1,4)],nrows=(file_rows_num-11))
        df_tmp['model_id'] = i
        df_tmp['model_name'] = models_list[i]
        if (i==0):
            df = df_tmp[(df_tmp.IR - df_tmp.DR) == 0]
            n_drop = len(df_tmp[(df_tmp.IR - df_tmp.DR) != 0])
        else:
            df_tmp = df_tmp.iloc[n_drop:]
            df = pd.concat((df,df_tmp))
    #    del df_tmp
        df_meta = pd.read_csv(file_name, delimiter='\t',
                              nrows=2)
#df = df.dropna(axis='columns')
    df = df.reset_index(drop=True)
    return df, df_meta

### Bayes factor calculations
    
def sort_df_by_dist(df,dist_col_name='DISTANCE'):
    return df.sort_values(by=dist_col_name,ascending=True)

def calc_bayes_factor(df_sort,num_top=100,
                      dist_col_name='DISTANCE',mode_id_col_name='model_id'):
    # For two models comparison, assuming model id=0/1
    num_ones = np.sum(df_sort.head(num_top)[mode_id_col_name])
    num_zeros = num_top - num_ones
    bayes_factor = num_ones / num_zeros  
    bayes_prop = num_ones / num_top
    return bayes_factor, bayes_prop

def plot_bayes_factor_vs_n(df_sort,b_num_vec,
                      dist_col_name='DISTANCE',mode_id_col_name='model_id'):
    #b_num_vec = np.array([1,3,5,7,10,20,30,50,70,100,300,500,700,1E3,3E3,7E3,1E4,3E4,len(df_bayes)]).astype(int)
    b_frac_vec = np.zeros_like(b_num_vec).astype(float)
    b_fact_vec = np.zeros_like(b_num_vec).astype(float)
    
    for ind,b_num in enumerate(b_num_vec):
        b_fact_vec[ind], b_frac_vec[ind]  = calc_bayes_factor(df_sort,num_top=b_num,
                      dist_col_name=dist_col_name,mode_id_col_name=mode_id_col_name)
    
    plt.figure()
    plt.semilogx(b_num_vec,b_frac_vec,'bx-') 
    plt.semilogx(b_num_vec,0.5* np.ones_like(b_num_vec).astype(float),'r--')   
    plt.xlabel('Number of nn to use',fontsize=14)
    plt.ylabel('Fraction of IR!=DR model in top nn',fontsize=14)
    plt.ylim([0,1])
    plt.show()
    
    plt.figure()
    plt.semilogx(b_num_vec,b_fact_vec,'x-')  
    plt.xlabel('Number of nn to use',fontsize=14)
    plt.ylabel('Bayes factor - n_difID/n_eqID',fontsize=14)
    plt.show()
    
    return

def abc_param_estimation(df_sort,num_top=100,models_list=['ideq','iddif'],
                         mode_name_col_name='model_name'):
    # Assumes for now 2 models
    df_model1_res = df_sort[df_sort[mode_name_col_name] == models_list[0]].head(num_top).mean()
    df_model2_res = df_sort[df_sort[mode_name_col_name] == models_list[1]].head(num_top).mean()
    df_model_combo_res = df_sort.head(num_top).mean()
    return df_model1_res, df_model2_res, df_model_combo_res

### https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html#sphx-glr-auto-examples-model-selection-plot-confusion-matrix-py

def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax

### NN classification functions

def data_preperation(df, df_meta=None, verbose=1, test_fraction=0.25,
                     val_fraction=0.05,y_col_name='model_id',
                     model_name_col='model_name',
                     y_col_reg_list = ['RL', 'AIR', 'ADR', 'IR', 'DR'],
                     num_first_col_to_exclude=7,
                     models_list=['ideq','iddif'],seed_num=3,
                     large_meta_flag=False):
#    y_col_name = 'model_id'
    exclude_col_list = list(df.columns[:num_first_col_to_exclude])
    exclude_col_list.append(y_col_name)
    x_col_list = [x for x in df.columns if x not in exclude_col_list and model_name_col not in x]
    #X_train_orig,X_test_orig,Y_train_orig,Y_test_orig = train_test_split(df[x_col_list],
    #                                            df[y_col_name],
    #                                            test_size=0.3,random_state=3)
    
    df_len = len(df)
#    test_fraction = 0.25
#    val_fraction = 0.05
    train_fraction = 1 - test_fraction - val_fraction
    indices = list(range(df_len))
    random.Random(seed_num).shuffle(indices)
    indices_test = indices[:int(test_fraction*df_len)]
    indices_val = indices[int(test_fraction*df_len):int((test_fraction+val_fraction)*df_len)]
    indices_train = indices[int((test_fraction+val_fraction)*df_len):]
    X_train_orig = df.iloc[indices_train][x_col_list]
    X_val_orig = df.iloc[indices_val][x_col_list]
    X_test_orig = df.iloc[indices_test][x_col_list]
    if df_meta is not None:
        if large_meta_flag:
            X_data_orig = df_meta[x_col_list]
        else:
            X_data_orig = df_meta.iloc[0][x_col_list]
    Y_train_orig = df.iloc[indices_train][y_col_name]
    Y_val_orig = df.iloc[indices_val][y_col_name]
    Y_test_orig = df.iloc[indices_test][y_col_name]
    Y_train_orig_reg = df.iloc[indices_train][y_col_reg_list]
    Y_val_orig_reg = df.iloc[indices_val][y_col_reg_list]
    Y_test_orig_reg = df.iloc[indices_test][y_col_reg_list]

    
    
    
    if verbose:
        print( len(Y_train_orig), len(Y_test_orig))
    X_train_orig = X_train_orig.values
    X_test_orig = X_test_orig.values
    X_val_orig = X_val_orig.values
    if df_meta is not None:
        if large_meta_flag:
            X_data_orig = X_data_orig.values    
        else:
            X_data_orig = np.array(X_data_orig.values).astype(float).reshape(1,-1)
    Y_train_orig = Y_train_orig.values
    Y_test_orig = Y_test_orig.values
    Y_val_orig = Y_val_orig.values
    Y_train_orig_reg = Y_train_orig_reg.values
    Y_test_orig_reg = Y_test_orig_reg.values
    Y_val_orig_reg = Y_val_orig_reg.values
    
    #
    X_train_mean = np.average(X_train_orig,axis=0)
    X_train_std  = np.std(X_train_orig,axis=0)
    X_train = (X_train_orig - X_train_mean)/X_train_std
    X_test = (X_test_orig - X_train_mean)/X_train_std
    X_val = (X_val_orig - X_train_mean)/X_train_std
    if df_meta is not None:
        X_data = (X_data_orig - X_train_mean)/X_train_std
    else: 
        X_data = None
    
    # Reshape
    Y_train_mean = 0 #np.average(Y_train_orig,axis=0) #classification
    Y_train_std  = 1 #np.std(Y_train_orig,axis=0) #classification
    Y_train_reg_mean = np.average(Y_train_orig_reg,axis=0).reshape(1,-1)
    Y_train_reg_std  = np.std(Y_train_orig_reg,axis=0).reshape(1,-1)
    Y_train = (Y_train_orig - Y_train_mean)/Y_train_std
    Y_train = Y_train.reshape(-1,1)
    Y_test = (Y_test_orig - Y_train_mean)/Y_train_std
    Y_test = Y_test.reshape(-1,1)
    Y_val = (Y_val_orig - Y_train_mean)/Y_train_std
    Y_val = Y_val.reshape(-1,1)
    Y_train_reg = (Y_train_orig_reg - Y_train_reg_mean)/Y_train_reg_std
    Y_test_reg = (Y_test_orig_reg - Y_train_reg_mean)/Y_train_reg_std
    Y_val_reg = (Y_val_orig_reg - Y_train_reg_mean)/Y_train_reg_std

    # Converting to one hot encoding
    b = np.zeros((len(Y_test),len(models_list)))
    b[np.arange(len(Y_test)), Y_test.flatten().astype(int)] = 1
    Y_test = b
    b = np.zeros((len(Y_train),len(models_list)))
    b[np.arange(len(Y_train)), Y_train.flatten().astype(int)] = 1
    Y_train = b
    b = np.zeros((len(Y_val),len(models_list)))
    b[np.arange(len(Y_val)), Y_val.flatten().astype(int)] = 1
    Y_val = b
    #Y_train = Y_train_orig.T
    #Y_test = Y_test_orig.T
    if verbose:
        print ("number of training examples = " + str(X_train.shape[0]))
        print ("number of test examples = " + str(X_test.shape[0]))
        print ("X_train shape: " + str(X_train.shape))
        print ("Y_train shape: " + str(Y_train.shape))
        print ("X_test shape: " + str(X_test.shape))
        print ("Y_test shape: " + str(Y_test.shape))
        
    return X_train, Y_train, X_val, Y_val, X_test, Y_test, X_data, Y_train_reg, Y_test_reg, Y_val_reg, Y_train_reg_mean, Y_train_reg_std

def NNModel_class(input_shape,num_classes=2):
    X_input = Input(input_shape)  # Define the input placeholder as a tensor with shape input_shape. Think of this as your input image!
    X = Dense(input_shape[0], activation='relu', name='fc1')(X_input)
    X = Dense(10, activation='relu', name='fc2')(X)
    X = Dense(5, activation='relu', name='fc3')(X)
#    X = Dense(1, activation='sigmoid', name='fc_out')(X)
    X = Dense(num_classes, activation='sigmoid', name='fc_out')(X)
    model = Model(inputs = X_input, outputs = X, name='ClassModelSelection')  # Create model. This creates your Keras model instance, you'll use this instance to train/test the model.
    return model

def NNModel_reg(input_shape,num_param=5):
    X_input = Input(input_shape)  # Define the input placeholder as a tensor with shape input_shape. Think of this as your input image!
    X = Dense(input_shape[0], activation='relu', name='fc1')(X_input)
    X = Dense(12, activation='relu', name='fc2')(X)
    X = Dense(8, activation='relu', name='fc3')(X)
#    X = Dense(1, activation='sigmoid', name='fc_out')(X)
    X = Dense(num_param, activation='linear', name='fc_out')(X)
    model = Model(inputs = X_input, outputs = X, name='RegModel')  # Create model. This creates your Keras model instance, you'll use this instance to train/test the model.
    return model

def res_vec_to_metrics(Y_pred,reg_pred = 1E-3):
    Y_pred_class = np.argmax(Y_pred,axis=1)
    Y_pred_arg = np.max(Y_pred,axis=1) 
    Y_pred_fac = Y_pred_arg/(np.min(Y_pred,axis=1)+reg_pred)
    return Y_pred_class, Y_pred_arg, Y_pred_fac


#%% Calculating ABC mean, ABC ridge, nn, len_stats (TODO)

def calc_abc_mean_stats(df_sort,
                        b_num_top=100,bayes_combo_flag=False):
    """
    Inputs:
        df_sort - df sorted by distance
    Output:
        res_dict including calculated metrics
    """
    res_dict={}
    
    # df_model1_res, df_model2_res, df_model_combo_res = msf.abc_param_estimation(df_sort=df_sort)
    df_model1_res, df_model2_res, df_model_combo_res = abc_param_estimation(df_sort=df_sort)
    param_list = ['RL', 'AIR', 'ADR', 'IR', 'DR']
    eq_model_std = df_sort[df_sort['model_id'] == 0].head(b_num_top)[param_list].std()
    dif_model_std = df_sort[df_sort['model_id'] == 1].head(b_num_top)[param_list].std()
    # if bayes_combo_flag: #combining eq and dif
    df_bayes_res = df_model_combo_res
    b_std = df_sort.head(b_num_top)[param_list].std()
    # elif bayes_prop<0.5: #eq
    #     df_bayes_res = df_model1_res 
    #     b_std = df_sort[df_sort['model_id'] == 0].head(b_num_top)[param_list].std()
    # else: #dif
    #     df_bayes_res = df_model2_res
    #     b_std = df_sort[df_sort['model_id'] == 1].head(b_num_top)[param_list].std()
    # m stands for calculated by mean abc
    res_dict = {    'm_combo_RL':[df_bayes_res.RL], 'm_combo_AIR':[df_bayes_res.AIR],
                    'm_combo_ADR':[df_bayes_res.ADR],'m_combo_IR':[df_bayes_res.IR],
                    'm_combo_DR':[df_bayes_res.DR],
                    'm_combo_RL_std':[b_std['RL']],'m_combo_AIR_std':[b_std['AIR']],
                    'm_combo_ADR_std':[b_std['ADR']],'m_combo_IR_std':[b_std['IR']],
                    'm_combo_DR_std':[b_std['DR']],
                    'm_eq_RL':[df_model1_res.RL], 'm_eq_AIR':[df_model1_res.AIR],
                    'm_eq_ADR':[df_model1_res.ADR],'m_eq_IR':[df_model1_res.IR],
                    'm_eq_DR':[df_model1_res.DR],
                    'm_dif_RL':[df_model2_res.RL], 'm_dif_AIR':[df_model2_res.AIR],
                    'm_dif_ADR':[df_model2_res.ADR],'m_dif_IR':[df_model2_res.IR],
                    'm_dif_DR':[df_model2_res.DR],
                    'm_eq_RL_std':[eq_model_std['RL']],'m_eq_AIR_std':[eq_model_std['AIR']],
                    'm_eq_ADR_std':[eq_model_std['ADR']],'m_eq_IR_std':[eq_model_std['IR']],
                    'm_eq_DR_std':[eq_model_std['DR']],
                    'm_dif_RL_std':[dif_model_std['RL']],'m_dif_AIR_std':[dif_model_std['AIR']],
                    'm_dif_ADR_std':[dif_model_std['ADR']],'m_dif_IR_std':[dif_model_std['IR']],
                    'm_dif_DR_std':[dif_model_std['DR']],}
    # print('bla')
    # print(res_dict)
    return res_dict

#def calc_abc_ridge_stats(df_sort,df_meta,
#                        b_num_top=100):
#    param_list = ['RL', 'AIR', 'ADR', 'IR', 'DR']
#    df_sort_eq_model = df_sort[df_sort['model_id'] == 0].head(b_num_top)
#    df_sort_dif_model = df_sort[df_sort['model_id'] == 1].head(b_num_top)
#    df_sort_combo_model = df_sort.head(b_num_top)
    
#    # df_eq
#    res_eq = ridge_for_abc(df_sort=df_sort_eq_model,df_meta=df_meta,num_bayes=b_num_top,
#                  params=param_list,
#                  cols_to_exclude=['DISTANCE','model_id','model_name'],
#                  norm_flag=True,prefix='r_eq_')
#    res_dif = ridge_for_abc(df_sort=df_sort_dif_model,df_meta=df_meta,num_bayes=b_num_top,
#                  params=param_list,
#                  cols_to_exclude=['DISTANCE','model_id','model_name'],
#                  norm_flag=True,prefix='r_dif_')
#    res_combo = ridge_for_abc(df_sort=df_sort_combo_model,df_meta=df_meta,num_bayes=b_num_top,
#                  params=param_list,
#                  cols_to_exclude=['DISTANCE','model_id','model_name'],
#                  norm_flag=True,prefix='r_combo_')
#    res_out = res_eq
#    res_out.update(res_dif)
#    res_out.update(res_combo)
#    return res_out

def calc_abc_ridge_stats(df_sort,df_meta,
                        b_num_top=100):
    param_list = ['RL', 'AIR', 'ADR', 'IR', 'DR']
    df_sort_eq_model = df_sort[df_sort['model_id'] == 0].head(b_num_top)
    df_sort_dif_model = df_sort[df_sort['model_id'] == 1].head(b_num_top)
    df_sort_combo_model = df_sort.head(b_num_top)
    
    # df_eq
    res_eq = ridge_for_abc(df_sort=df_sort_eq_model,df_meta=df_meta,num_bayes=b_num_top,
                  params=param_list,
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='r_eq_')
    res_eq_l = lasso_for_abc(df_sort=df_sort_eq_model,df_meta=df_meta,num_bayes=b_num_top,
                  params=param_list,
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='l_eq_')
    res_eq_p = pls_for_abc(df_sort=df_sort_eq_model, df_meta=df_meta, num_bayes=b_num_top,
                             params=param_list,
                             cols_to_exclude=['DISTANCE', 'model_id', 'model_name'],
                             norm_flag=True, prefix='p_eq_')
    res_dif = ridge_for_abc(df_sort=df_sort_dif_model,df_meta=df_meta,num_bayes=b_num_top,
                  params=param_list,
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='r_dif_')
    res_dif_l = lasso_for_abc(df_sort=df_sort_dif_model,df_meta=df_meta,num_bayes=b_num_top,
                  params=param_list,
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='l_dif_')
    res_dif_p = pls_for_abc(df_sort=df_sort_eq_model, df_meta=df_meta, num_bayes=b_num_top,
                             params=param_list,
                             cols_to_exclude=['DISTANCE', 'model_id', 'model_name'],
                             norm_flag=True, prefix='p_dif_')
    res_combo = ridge_for_abc(df_sort=df_sort_combo_model,df_meta=df_meta,num_bayes=b_num_top,
                  params=param_list,
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='r_combo_')
    res_combo_l = lasso_for_abc(df_sort=df_sort_combo_model,df_meta=df_meta,num_bayes=b_num_top,
                  params=param_list,
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='l_combo_')
    res_out = res_eq
    res_out.update(res_dif)
    res_out.update(res_combo)
    res_out.update(res_eq_l)
    res_out.update(res_dif_l)
    res_out.update(res_combo_l)
    res_out.update(res_eq_p)
    res_out.update(res_dif_p)
    return res_out    

    
def ridge_for_abc(df_sort,df_meta,num_bayes=100,
                  params=['RL','AIR','ADR','IR','DR'],
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='r_'):
    # df_ridge = df_sort.head(num_bayes_vec[0])
    # params = ['RL','AIR','ADR','IR','DR']
    res_out = {}
    df_ridge = df_sort.head(num_bayes)
    y = np.float64(df_ridge[params].values)
    cols_to_exclude = ['DISTANCE','model_id','model_name']
    x_cols = [col for col in df_ridge.columns if (col not in params) and (col not in cols_to_exclude)]
    x_norm = df_meta.iloc[1][x_cols].values
    if norm_flag:
        x = np.float64(df_ridge[x_cols].values*x_norm)
    else:
        x = np.float64(df_ridge[x_cols].values)
    x_to_pred = np.float64(df_meta.iloc[0][x_cols].values*x_norm)
    # y_hat_list = []
    for ind,p in enumerate(params):
        # reg = linear_model.RidgeCV(alphas=np.logspace(-4,4,20), cv=10,
        #                             fit_intercept=True, normalize=False,store_cv_values=True)      
        reg = linear_model.Ridge(normalize=False,fit_intercept=True)
        parameters = {'alpha':np.logspace(-4,4,20)}
        # clf = model_selection.GridSearchCV(estimator = reg,
        #                                    param_grid = parameters, cv=10,
        #                                    scoring = 'neg_root_mean_squared_error')
        # clf.fit(x,y[:,ind])
        # # reg.fit(x,y[:,ind])
        # cv_rmse = np.min(-clf.cv_results_['mean_test_score'])
        clf = model_selection.GridSearchCV(estimator = reg,
                                   param_grid = parameters, cv=10,
                                   scoring = 'neg_mean_squared_error')
        clf.fit(x,y[:,ind])
        # reg.fit(x,y[:,ind])
        cv_rmse = np.min(np.sqrt(-clf.cv_results_['mean_test_score']))
        train_rmse = np.std(clf.predict(x) - y[:,ind])
        best_lambda = clf.best_params_['alpha']
        res = clf.predict(x_to_pred.reshape(1,-1))
        res_out[prefix+p] = [res[0]]
        res_out[prefix+p+'_cv_rmse'] = [cv_rmse]
        res_out[prefix+p+'_train_rmse'] = [train_rmse]
        res_out[prefix+p+'_best_lambda'] = [best_lambda]
        # plt.figure()
        # plt.plot(y[:,ind],clf.predict(x),'x')
        # plt.show()
        # print(cv_rmse,train_rmse, best_lambda)
        # cv_mse = np.mean(reg.cv_values_, axis=0)
        
    return res_out

def lasso_for_abc(df_sort,df_meta,num_bayes=100,
                  params=['RL','AIR','ADR','IR','DR'],
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='l_'):
    # df_ridge = df_sort.head(num_bayes_vec[0])
    # params = ['RL','AIR','ADR','IR','DR']
    res_out = {}
    df_lasso = df_sort.head(num_bayes)
    y = np.float64(df_lasso[params].values)
    cols_to_exclude = ['DISTANCE','model_id','model_name']
    x_cols = [col for col in df_lasso.columns if (col not in params) and (col not in cols_to_exclude)]
    x_norm = df_meta.iloc[1][x_cols].values
    if norm_flag:
        x = np.float64(df_lasso[x_cols].values*x_norm)
    else:
        x = np.float64(df_lasso[x_cols].values)
    x_to_pred = np.float64(df_meta.iloc[0][x_cols].values*x_norm)
    # y_hat_list = []
    for ind,p in enumerate(params):
        # reg = linear_model.RidgeCV(alphas=np.logspace(-4,4,20), cv=10,
        #                             fit_intercept=True, normalize=False,store_cv_values=True)      
        reg = linear_model.Lasso(normalize=False,fit_intercept=True)
        parameters = {'alpha':np.logspace(-9,2,20)}
        # clf = model_selection.GridSearchCV(estimator = reg,
        #                                    param_grid = parameters, cv=10,
        #                                    scoring = 'neg_root_mean_squared_error')
        # clf.fit(x,y[:,ind])
        # # reg.fit(x,y[:,ind])
        # cv_rmse = np.min(-clf.cv_results_['mean_test_score'])
        clf = model_selection.GridSearchCV(estimator = reg,
                                   param_grid = parameters, cv=10,
                                   scoring = 'neg_mean_squared_error')
        clf.fit(x,y[:,ind])
        # reg.fit(x,y[:,ind])
        cv_rmse = np.min(np.sqrt(-clf.cv_results_['mean_test_score']))
        train_rmse = np.std(clf.predict(x) - y[:,ind])
        best_lambda = clf.best_params_['alpha']
        res = clf.predict(x_to_pred.reshape(1,-1))
        res_out[prefix+p] = [res[0]]
        res_out[prefix+p+'_cv_rmse'] = [cv_rmse]
        res_out[prefix+p+'_train_rmse'] = [train_rmse]
        res_out[prefix+p+'_best_lambda'] = [best_lambda]
        # plt.figure()
        # plt.plot(y[:,ind],clf.predict(x),'x')
        # plt.show()
        # print(cv_rmse,train_rmse, best_lambda)
        # cv_mse = np.mean(reg.cv_values_, axis=0)
        
    return res_out

def pls_for_abc(df_sort,df_meta,num_bayes=100,
                  params=['RL','AIR','ADR','IR','DR'],
                  cols_to_exclude=['DISTANCE','model_id','model_name'],
                  norm_flag=True,prefix='p_'): #TODO https://www.statology.org/partial-least-squares-in-python/
    # df_ridge = df_sort.head(num_bayes_vec[0])
    # params = ['RL','AIR','ADR','IR','DR']
    max_num_pls_comps = 10
    res_out = {}
    df_lasso = df_sort.head(num_bayes)
    y = np.float64(df_lasso[params].values)
    cols_to_exclude = ['DISTANCE','model_id','model_name']
    x_cols = [col for col in df_lasso.columns if (col not in params) and (col not in cols_to_exclude)]
    x_norm = df_meta.iloc[1][x_cols].values
    if norm_flag:
        x = np.float64(df_lasso[x_cols].values*x_norm)
    else:
        x = np.float64(df_lasso[x_cols].values)
    x_to_pred = np.float64(df_meta.iloc[0][x_cols].values*x_norm)

    cv = model_selection.RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
    mse = []
    n = len(x)

    # Calculate MSE with only the intercept
    # score = -1 * model_selection.cross_val_score(skl   PLSRegression(n_components=1),
    #                                              np.ones((n, 1)), y, cv=cv, scoring='neg_mean_squared_error').mean()
    # mse.append(score)

    # Calculate MSE using cross-validation, adding one component at a time
    for i in np.arange(1, max_num_pls_comps):
        pls = PLSRegression(n_components=i)
        score = -1 * model_selection.cross_val_score(pls, x, y, cv=cv,
                                                     scoring='neg_mean_squared_error').mean()
        mse.append(score)
    chosen_num_comp = sorted(np.arange(1,max_num_pls_comps), key=lambda t: mse[t-1])[0]
    pls = PLSRegression(n_components=chosen_num_comp)
    pls.fit(x, y)
    y_train_hat = pls.predict(x)
    mse_train = ((y_train_hat - y)**2).sum(axis=0)**0.5
    y_hat = pls.predict(x_to_pred.reshape(1, -1)).flatten()
    # y_hat_list = []
    for ind, p in enumerate(params):
        # reg = linear_model.RidgeCV(alphas=np.logspace(-4,4,20), cv=10,

        res_out[prefix+p] = [y_hat[ind]]
        # res_out[prefix+p+'_cv_rmse'] = [cv_rmse]
        res_out[prefix+p+'_train_rmse'] = [mse_train[ind]]
        res_out[prefix+p+'_num_comp'] = [chosen_num_comp]
    return res_out

def nn_class_and_reg(df,df_meta,models_list=['ideq','iddif'],
               num_epochs = 50, batch_size=4096,verbose=0):
        # X_train, Y_train, X_val, Y_val, X_test, Y_test, X_data, Y_train_reg, Y_test_reg, Y_val_reg, Y_train_reg_mean, Y_train_reg_std = msf.data_preperation(df=df,
        #                                                                           df_meta=df_meta,
        #                                                                          verbose=verbose)
        X_train, Y_train, X_val, Y_val, X_test, Y_test, X_data, Y_train_reg, Y_test_reg, Y_val_reg, Y_train_reg_mean, Y_train_reg_std = data_preperation(df=df,
                                                                          df_meta=df_meta,
                                                                         verbose=verbose)
        # ms_model = msf.NNModel_class(X_train.shape[1:],num_classes=len(models_list))
        ms_model = NNModel_class(X_train.shape[1:],num_classes=len(models_list))
        #ms_model.compile('adam', 'categorical_crossentropy', metrics=['accuracy'])
        ms_model.compile('SGD', 'categorical_crossentropy', metrics=['accuracy']) #Changed optimizer, as adam was not stable.
        ms_model.fit(X_train, Y_train, epochs=num_epochs, batch_size=batch_size,
                      validation_data=(X_val, Y_val),verbose=verbose)
        Y_pred_train = ms_model.predict(X_train) #predict on train set
        Y_pred_test = ms_model.predict(X_test) # Predict test set
        Y_pred_data = ms_model.predict(X_data)
        nn_c_class = 'eq' if Y_pred_data[0][0]>Y_pred_data[0][1] else 'dif'
        nn_c_test_max_out_mean = Y_pred_test.max(axis=1).mean()
        nn_c_test_max_out_std = Y_pred_test.max(axis=1).std()
        # Y_pred_class_train, Y_pred_arg_train, Y_pred_fac_train = msf.res_vec_to_metrics(Y_pred=Y_pred_train)
        # Y_pred_class_test, Y_pred_arg_test, Y_pred_fac_test = msf.res_vec_to_metrics(Y_pred=Y_pred_test)
        # Y_pred_class_data, Y_pred_arg_data, Y_pred_fac_data = msf.res_vec_to_metrics(Y_pred=Y_pred_data)
        # Y_real_class_test, _, _ = msf.res_vec_to_metrics(Y_pred=Y_test)
        Y_pred_class_train, Y_pred_arg_train, Y_pred_fac_train = res_vec_to_metrics(Y_pred=Y_pred_train)
        Y_pred_class_test, Y_pred_arg_test, Y_pred_fac_test = res_vec_to_metrics(Y_pred=Y_pred_test)
        Y_pred_class_data, Y_pred_arg_data, Y_pred_fac_data = res_vec_to_metrics(Y_pred=Y_pred_data)
        Y_real_class_test, _, _ = res_vec_to_metrics(Y_pred=Y_test)
        cm = confusion_matrix(y_true=Y_real_class_test,
                              y_pred=Y_pred_class_test)/(len(Y_real_class_test)/2)
        if verbose:
            # msf.plot_confusion_matrix(y_true=Y_real_class_test,y_pred=Y_pred_class_test, 
            #                   classes=np.array(models_list), normalize=True)
            plot_confusion_matrix(y_true=Y_real_class_test,y_pred=Y_pred_class_test, 
                              classes=np.array(models_list), normalize=True)
            
        # reg_model = msf.NNModel_reg(X_train.shape[1:],num_param=Y_train_reg.shape[1])
        reg_model = NNModel_reg(X_train.shape[1:],num_param=Y_train_reg.shape[1])
        #reg_model.compile('adam', 'mean_squared_error', metrics=['accuracy'])
        reg_model.compile('SGD', 'mean_squared_error', metrics=['accuracy']) #Changed optimizer, as adam was not stable.
        reg_model.fit(X_train, Y_train_reg, epochs=num_epochs, batch_size=batch_size,
                      validation_data=(X_val, Y_val_reg),verbose=verbose)
    #    Y_pred_train_reg = reg_model.predict(X_train)*Y_train_reg_std + Y_train_reg_mean #predict on train set
        Y_pred_test_reg = reg_model.predict(X_test)*Y_train_reg_std + Y_train_reg_mean # Predict test set
        Y_pred_data_reg = reg_model.predict(X_data)*Y_train_reg_std + Y_train_reg_mean
        nn_std_err = np.std(Y_pred_test_reg - (Y_test_reg*Y_train_reg_std + Y_train_reg_mean),axis=0)
        # Gathering results
        res_dict = {'nn_c_out_eq':[Y_pred_data[0][0]],
                    'nn_c_out_dif':[Y_pred_data[0][1]], 
                    'nn_c_class':[nn_c_class],
                    'nn_c_test_max_out_mean':[nn_c_test_max_out_mean],
                    'nn_c_test_max_out_std':[nn_c_test_max_out_std],
                    'nn_c_cm00':[cm[0][0]],
                    'nn_c_cm01':[cm[0][1]], 'nn_c_cm10':[cm[1][0]],'nn_c_cm11':[cm[1][1]],
                    'nn_r_corr_RL':[pearsonr(Y_test_reg[:,0], Y_pred_test_reg[:,0])[0]],
                    'nn_r_corr_AIR':[pearsonr(Y_test_reg[:,1], Y_pred_test_reg[:,1])[0]],
                    'nn_r_corr_ADR':[pearsonr(Y_test_reg[:,2], Y_pred_test_reg[:,2])[0]],
                    'nn_r_corr_IR':[pearsonr(Y_test_reg[:,3], Y_pred_test_reg[:,3])[0]],
                    'nn_r_corr_DR':[pearsonr(Y_test_reg[:,4], Y_pred_test_reg[:,4])[0]],
                    'nn_r_RL':[Y_pred_data_reg[0,0]],'nn_r_AIR':[Y_pred_data_reg[0,1]],
                    'nn_r_ADR':[Y_pred_data_reg[0,2]],'nn_r_IR':[Y_pred_data_reg[0,3]],
                    'nn_r_DR':[Y_pred_data_reg[0,4]],
                    'nn_r_RL_test_stderr':[nn_std_err[0]],'nn_r_AIR_test_stderr':[nn_std_err[1]],
                    'nn_r_ADR_test_stderr':[nn_std_err[2]],'nn_r_IR_test_stderr':[nn_std_err[3]],
                    'nn_r_DR_test_stderr':[nn_std_err[4]],}
        return res_dict


#%% Aggregating all the functions to create the data on a single lib
# results_path = 'results2/eggnog/' # relative to pwd
# csv_out_path = 'results_small/infer_abc_params/eggnog/'
# lib = 'Artropoda_ENOG41034WH' #'Artropoda_ENOG41034WH_dummy'#'PF00747'   
# verbose = 1
# models_list=['ideq','iddif']
# bayes_combo_flag = False
# num_epochs = 50
# batch_size=4096
# b_num_top=100
        
# relative to pwd
# csv_out_path = results_path = pipeline_path+'results/example4/'
# lib = 'data_name' #'Artropoda_ENOG41034WH_dummy'#'PF00747'   
# verbose = 1
# models_list=['ideq','iddif']
# bayes_combo_flag = False
# num_epochs = 50
# batch_size=4096
# b_num_top=100

# calc_stats(lib=lib,path='' ,csv_out_path=csv_out_path,
#            verbose = 1 , models_list=['ideq','iddif'],
#            num_epochs = 50, batch_size=4096,b_num_top=100)

def calc_stats(csv_out_path='results_small/infer_abc_params/',lib='data_name',path='',     
               verbose = 0, models_list=['ideq','iddif'],
               num_epochs = 50, batch_size=4096,b_num_top=100):
    """
    Inputs: 
        lib - lib path # need to modify
        models_list - currently supported 'ideq','iddif'
    """
    logging.basicConfig(filename=csv_out_path+'infer_abc_params_single_folder.log',level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S') 
    res_dict = {}
    res_dict['dir_name'] = [lib]
    print(lib) # For debugging purpose
    logging.info('started calc stats '+lib)
    # Add a check if file size is big enough, add custom made path
    # df, df_meta = msf.load_lib_data(lib=lib) # maybe switch to other function
    # df_sort = msf.sort_df_by_dist(df=df)
    # df, df_meta = load_lib_data(lib=lib,pwd=pwd, rel_path=path) # maybe switch to other function
    df, df_meta = load_lib_data(path=csv_out_path,lib=lib,models_list=['ideq','iddif'] ,size_th=1E6,rel_path='')
    df = df[(df.DISTANCE != 100000)] #High distance that in new version is a sign of too long run
    if df is None:
        logging.warning(f'skipping {lib} due to small size of file.')
        return
    logging.info('ran load_lib_data '+lib)
    df_sort = sort_df_by_dist(df=df)
    # calculating bayes factor
    # bayes_factor, bayes_prop = msf.calc_bayes_factor(df_sort=df_sort,num_top=b_num_top)
    bayes_factor, bayes_prop = calc_bayes_factor(df_sort=df_sort,num_top=b_num_top)
    logging.info('ran calc_bayes_factor '+lib)
    res_dict['abc_num_top'] = [b_num_top]
    res_dict['bayes_factor'] = [bayes_factor]
    res_dict['bayes_prop'] = [bayes_prop]
    res_dict['bayes_class'] = ['eq'] if bayes_prop<0.5 else ['dif']
    logging.info('done bayes class '+lib)    
    
    #ABC mean
    tmp_res_dict = calc_abc_mean_stats(df_sort=df_sort, b_num_top=b_num_top) 
    res_dict.update(tmp_res_dict)
    logging.info('done ABC mean '+lib)
    #ABC ridge
    tmp_res_dict = calc_abc_ridge_stats(df_sort=df_sort,df_meta=df_meta,b_num_top=b_num_top)
    res_dict.update(tmp_res_dict)
    logging.info('done ABC ridge '+lib)
    #nn
    tmp_res_dict = nn_class_and_reg(df=df,df_meta=df_meta,
                                    models_list=models_list,
                                    num_epochs = num_epochs,
                                    batch_size=batch_size,verbose=verbose)
    res_dict.update(tmp_res_dict)
    logging.info('done nn class. and reg. '+lib)
    
    #
   
    df_res = pd.DataFrame(res_dict)
    df_res = df_res.reindex(sorted(df_res.columns), axis=1)
    df_res.to_csv(csv_out_path+lib+'_res.csv')
    logging.info(f'Done {lib}. Res file written to {csv_out_path+lib+"_res.csv"}')
    return

#%% click command

@click.command()
@click.option('--lib', default='', help='Data directory name.')
@click.option('--pwd', default='/sparta/', help='pwd')
@click.option('--path', default='' ,help='Where the data directoy relative to pwd')
@click.option('--csv_out_path', default='results_small/infer_abc_params/', help='csv out path relative to pwd')
@click.option('--verbose', default=0 ,help='Verbosity (0/1)')
@click.option('--ne', default=50, help='NN num epochs')
@click.option('--bs', default=4096 ,help='NN batch size')
@click.option('--bn', default=100 ,help='epsilon num to use for ABC')
def calc_stats_click(lib,pwd,path,csv_out_path,verbose,
               ne, bs,bn):
        calc_stats(lib=str(lib), pwd=str(pwd),path=str(path) ,
                   csv_out_path=str(csv_out_path),
                   verbose=int(verbose),models_list=['ideq','iddif'],
                   num_epochs=int(ne), batch_size=int(bs),b_num_top=int(bn))
# bn = 100
# bs = 4096 
# ne = 50       
# path = 'results2/eggnog/' # relative to pwd
# csv_out_path = 'results_small/infer_abc_params/eggnog/'
# lib = 'Artropoda_ENOG41034WH' #'Artropoda_ENOG41034WH_dummy'#'PF00747'   
# verbose = 1

# models_list=['ideq','iddif']
# bayes_combo_flag = False
# num_epochs = 50
# batch_size=4096
# b_num_top=100
        
# Command example
# python infer_abc_params_single_folder.py --lib Artropoda_ENOG41034WH --path results2/eggnog/ --csv_out_path results_small/infer_abc_params/eggnog/ --verbose 1
#%% main
        
if __name__ == '__main__':
    calc_stats_click()

