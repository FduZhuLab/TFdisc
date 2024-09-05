import os
import numpy as np
import pandas as pd
import pickle
import time
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error,r2_score
import joblib
from sklearn.ensemble import RandomForestRegressor
from sklearn.kernel_ridge import KernelRidge
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")


def TF_krr(imp_data,GRN,gene,save_path,krr_test_size=0.1,score=False):
    
    """
    Train prediction model(KRR)
    
    Arguments:
    imp_data: a pandas DataFrame of gene expression data(imputation data), row(cell), col(gene).
    GRN: Gene regulatory links
    gene: Gene for predictive model construction
    save_path: The path of model. 
    krr_test_size: The proportion of test data sets
    score: Whether to save the model score
    Returns: Prediction krr model of gene
    
    """
    X, y = imp_data[GRN[gene]].values, imp_data[gene].values
    train_data, test_data, train_targets, test_targets = train_test_split(X, y, test_size=krr_test_size,random_state=42)
    param_grid = {"alpha": [1e0, 1e-1, 1e-2],
             "gamma":[1e0, 1e-1,0.02, 1e-2]}
    kr = GridSearchCV(KernelRidge(kernel='rbf'), param_grid=param_grid,n_jobs=-1)
    kr.fit(train_data, train_targets)
    joblib.dump(kr, save_path)
    rf_y_predict = kr.predict(test_data)
    rf_y_predict_train = kr.predict(train_data)
    if score:
        traing_score.append(kr.score(train_data, train_targets))
        testing_score.append(kr.score(test_data, test_targets))
        mse.append(mean_squared_error(test_targets, rf_y_predict))
        mae.append(mean_absolute_error(test_targets, rf_y_predict))
        mse_train.append(mean_squared_error(train_targets, rf_y_predict_train))
        mae_train.append(mean_absolute_error(train_targets, rf_y_predict_train))

def TF_rf(imp_data,GRN,gene,save_path,rf_test_size=0.1,score=False):
    
    """
    Train prediction model(KRR)
    
    Arguments:
    imp_data: a pandas DataFrame of gene expression data(imputation data), row(cell), col(gene).
    GRN: Gene regulatory links
    gene: Gene 
    save_path: The path of model. 
    rf_test_size: the proportion of test data sets
    score: Whether to save the model score
    
    Returns: Prediction krr model of gene
    
    """
    X, y = imp_data[GRN[gene]].values, imp_data[gene].values
    train_data, test_data, train_targets, test_targets = train_test_split(X, y, test_size=rf_test_size,random_state=42)
    rf1 = RandomForestRegressor(n_estimators=1000,bootstrap=True, random_state=42,oob_score=True,n_jobs=-1)
    rf1.fit(train_data, train_targets)
    importance1 = rf1.feature_importances_
    imp_save1 = np.argsort(importance1)[::-1][:len(GRN[gene])]
    raw_TF_list_save1 = GRN[gene]
    train_TF_list_save1 = []
    for i in range(len(imp_save1.tolist())):
        train_TF_list_save1.append(raw_TF_list_save1[imp_save1[i]])
    imp_save1.shape = (len(GRN[gene]),1)
    df = pd.DataFrame(importance1[imp_save1], index = train_TF_list_save1, columns =["importance"])
    #df.to_pickle(new_importance_path)
    joblib.dump(rf1, save_path)
    rf_y_predict = rf1.predict(test_data)
    rf_y_predict_train = rf1.predict(train_data)
    if score:
        traing_score.append(rf1.score(train_data, train_targets))
        testing_score.append(rf1.score(test_data, test_targets))
        oob_score.append(rf1.oob_score_)
        mse.append(mean_squared_error(test_targets, rf_y_predict))
        mae.append(mean_absolute_error(test_targets, rf_y_predict))
        mse_train.append(mean_squared_error(train_targets, rf_y_predict_train))
        mae_train.append(mean_absolute_error(train_targets, rf_y_predict_train))


def TF_model(imp_data,
             gene_name,
             GRN,
             save = None,
             method = "krr",
             verbose = True,
             test_size=0.1,
             model_score=False):
    """
    Train prediction model
    
    Arguments:
    imp_data: a pandas DataFrame of gene expression data(imputation data), row(cell), col(gene).
    gene_name: list of gene. If None or 'all', the list of gene names in expression_data will be used. 
    GRN: Gene regulatory links.
    save: The path of model.
    method: Method of prediction. "krr" or "rf".
    verbose: print info.
    test_size: the proportion of test data sets
    model_score: Whether to save the model score
    
    Returns: Prediction model of every gene
    
    """
    traing_score = []
    testing_score = []
    mse = []
    mae = []
    mse_train=[]
    mae_train=[]
    oob_score = []
    if isinstance(imp_data, pd.DataFrame):
        gene = list(imp_data.columns)
    if gene_name is None:
        gene_name = gene
    elif gene_name == 'all':
        gene_name = gene
    else:
        if len(gene_name) == 0:
            raise ValueError('gene names is empty')
        if not set(gene).intersection(set(gene_name)):
            raise ValueError('Intersection of gene and imp_data is empty.')
        if  set(gene).intersection(set(gene_name)):
            gene_name = list(set(gene) & set(gene_name))
    if save is None:
        if method == "krr":
            save = "./krr/"
            if not os.path.exists(save):
                os.makedirs(save)
        if method == "rf":
            save = "./rf/"
            if not os.path.exists(save):
                os.makedirs(save)
    else:
        if method == "krr":
            save = save+"/krr/"
            if not os.path.exists(save):
                os.makedirs(save)
        if method == "rf":
            save = save+"/rf/"
            if not os.path.exists(save):
                os.makedirs(save)
    if method == "krr":
        start = time.time()
        for i in tqdm(gene_name, desc="Task Progress", unit="items", ncols=80):
            path = save + "krr_" + i + ".m"
            TF_krr(imp_data,GRN,i,path,test_size,model_score)
        if model_score:
            df = pd.DataFrame({'traing_score': np.array(traing_score),
                               'testing_score': np.array(testing_score),
                               'mse': np.array(mse),
                               'mae': np.array(mae),
                               'mse_train': np.array(mse_train),
                               'mae_train': np.array(mae_train),
                               'gene': np.array(gene_name)})
            df.to_csv(save+"score.csv")
        if verbose:
                print("running time =",time.time()-start)
    if method == "rf":
        start = time.time()
        for i in tqdm(gene_name, desc="Task Progress", unit="items", ncols=80):
            start = time.time()
            path = save + "rf_" + i + ".m"
            TF_rf(imp_data,GRN,i,path,test_size,model_score)
        if model_score:
            df = pd.DataFrame({'traing_score': np.array(traing_score),
                               'testing_score': np.array(testing_score),
                               'mse': np.array(mse),
                               'mae': np.array(mae),
                               'mse_train': np.array(mse_train),
                               'mae_train': np.array(mae_train),
                               'oob_score': np.array(oob_score),
                               'gene': np.array(gene_name)})
            df.to_csv(save+"score.csv")
        if verbose:
                print("running time =",time.time()-start)
    
