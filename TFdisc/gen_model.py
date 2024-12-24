import os
import numpy as np
import pandas as pd
import pickle
import time
import joblib
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
from multiprocessing.pool import Pool
from functools import partial
from itertools import repeat
import warnings
warnings.filterwarnings("ignore")
<<<<<<< HEAD
from .pre_model import fix_TF_matrix_predict,fix_HVG_matrix_predict
=======
from pre_model import fix_TF_matrix_predict,fix_HVG_matrix_predict
>>>>>>> 144b3617b0f9ac7af893c5570c0093b8ddc5b51a


def combine_predict(imp_data,TF_list,HVG_list,grn,krr_premodel,rf_premodel,
                    TF=None,krr_time=5,rf_time=0,
                    core=10,matrix_err = 10000,min_matrix_err = 0.01) : 
    """
    Simulation prediction with a fixed number of iterations
    
    Arguments:
    # imp_data : imputation data which represents normalized data(imputation).  Rows represent cells, columns represent genes
    # TF_list : List of TFs that have been screened
    # grn : The dictionary form represents top TF that regulates target gene
    # krr_premodel : the folder location of pre-trained gene prediction models
    # rf_premodel : the folder location of pre-trained gene prediction models
    # TF : the simulate knockout TF,if TF=None it represents self-iteration
    # krr_time : the number of KRR prediction iterations
    # rf_time : the number of RF prediction iterations
    # core : The core of parallel computing
    # matrix_err : the initial value of matrix mse
    # min_matrix_err : the min value of matrix mse
    
    Returns: expression matrix
    
    """
    if(krr_time>0):
        TFdata = fix_TF_matrix_predict(imp_data,TF_list,grn,krr_premodel,save=None,
                      countlist=None,countnumber=krr_time,TF=None,core=core,
                      matrix_err = matrix_err,min_matrix_err = min_matrix_err,
                      err_save = None ,TF_save=None)
    if(rf_time>0):
        TFdata = fix_TF_matrix_predict(TFdata,TF_list,grn,rf_premodel,save=None,
                      countlist=None,countnumber=rf_time,TF=None,core=core,
                      matrix_err = matrix_err,min_matrix_err =min_matrix_err,
                      err_save = None ,TF_save=None)
    data = fix_HVG_matrix_predict(imp_data,TFdata,HVG_list,grn,krr_premodel,core=10)
    return data
    
    