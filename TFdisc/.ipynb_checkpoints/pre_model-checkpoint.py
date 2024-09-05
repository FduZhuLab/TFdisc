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
from multiprocessing.pool import Pool
from functools import partial
from itertools import repeat
import warnings
warnings.filterwarnings("ignore")

def gene_predict(gene,gene_premodel,model,TFdata,grn) :
    
    """
    single gene predict
    
    Arguments:
    # gene : gene name.
    # gene_premodel : the folder location of pre-trained gene prediction models.
    # TFdata_pd : TF data.
    # method: Method of prediction. "krr" or "rf".
    # grn : The dictionary form represents top TF that regulates target gene.
    
    Returns: Prediction values of gene
    
    """
    name = gene_premodel+method+"_"+gene+".m"
    model = joblib.load(name).set_params(n_jobs = 1)
    return model.predict(TFdata[grn[gene]].values)


def fix_TF_matrix_predict(imp_data,TF_list,grn,gene_premodel,save,
                      countlist,countnumber,TF=None,core=10,
                      matrix_err = 10000,min_matrix_err = 0.01,
                      err_save = None ,TF_save=None) : 
    """
    Simulation prediction with a fixed number of iterations
    
    Arguments:
    # imp_data : imputation data which represents normalized data(imputation).  Rows represent cells, columns represent genes
    # TF_list : List of TFs that have been screened
    # TF : the simulate knockout TF,if TF=None it represents self-iteration
    # grn : The dictionary form represents top TF that regulates target gene
    # countlist : The list of iterations number of the generated matrix
    # countnumber : the max number of iterations
    # gene_premodel : the folder location of pre-trained gene prediction models
    # core : The core of parallel computing
    # matrix_err : the initial value of matrix mse
    # min_matrix_err : the min value of matrix mse
    # save : the folder location of save TF knock-out simulation matrix
    # err_save : the folder location of err
    # TF_save : the folder location of TF
    
    Returns: TF expression matrix
    
    """
    start = time.time()
    TFdata_pd = imp_data[TF_list].copy()
    count = 0
    gene_errlist = []
    tmp_TF_list=[]
    matrix_errlist = []
    matrix_err_percent_list = []
    if TF is not None and TF in TF_list:
        TFdata_pd[TF]=0 
    while count< countnumber :
        pool = Pool(processes=core)
        a = list(pool.starmap(gene_predict, zip(TF_list, repeat(gene_premodel),
                                                repeat(TFdata_pd), repeat(grn)) ))
        pool.close()
        pool.join()
        inputdata_pd = pd.DataFrame(np.array(a).T, index=TFdata_pd.index, columns=TFdata_pd.columns)
        
        if (TF is not None) and (TF in TF_list):
            gene_err = np.sqrt(np.sum(inputdata_pd[TF].values**2))
            gene_errlist.append(gene_err)
            tmp_TF_list.append(inputdata_pd[TF])
            inputdata_pd[TF] = 0
        inputdata_pd[inputdata_pd<=0.0001]=0
        matrix_err = np.sqrt(np.sum((TFdata_pd.values - inputdata_pd.values)**2))
        matrix_err_percent = np.sqrt(np.sum((TFdata_pd.values - inputdata_pd.values)**2))/np.sqrt(np.sum((TFdata_pd.values)**2))
        matrix_errlist.append(matrix_err)
        matrix_err_percent_list.append(matrix_err_percent)
        TFdata_pd = inputdata_pd.copy()
        count = count + 1
        
        if (TF is not None) and (TF in TF_list):
            print("counts =",count,"    ","simulate_TF =",TF,"    ","matrix_err =",matrix_err,
                  "    ","matrix_err_percent =",matrix_err_percent, "    ",
                  "gene_err =",gene_err,"time =",time.time()-start)
            if (count in countlist) or (count >= countnumber):
                TFdata_pd.to_pickle(save+TF+"_TFdata"+str(count)+".pkl")
        else :
            print("counts =",count,"    ","simulate_TF =","cycle","    ","matrix_err =",matrix_err,
                  "    ","matrix_err_percent =",matrix_err_percent, "    ",
                  "time =",time.time()-start)
            if (count in countlist) or (count >= countnumber):
                TFdata_pd.to_pickle(save+"cycle_TFdata"+str(count)+".pkl")
                
    if (err_save is not None) and (TF is not None and TF in TF_list):
        pd.DataFrame({"gene_err":gene_errlist,
              "matrix_err":matrix_errlist,
              "matrix_err_percent":matrix_err_percent_list}).to_pickle(save+"err/"+TF+".pkl")
    if (TF_save is not None)  and (TF is not None and TF in TF_list): 
        pd.DataFrame(np.array(tmp_TF_list).T, index=TFdata_pd.index).to_pickle(save+"gene/"+TF+".pkl")
        
def fix_HVG_matrix_predict(imp_data,HVG_list,grn,gene_premodel,save,countnumber,TF=None,core=10) : 
    """
    Simulation prediction with a fixed number of iterations
    
    Arguments:
    # imp_data : imputation data which represents normalized data.  Rows represent cells, columns represent genes
    # HVG_list : List of HVG that have been screened
    # TF : the simulate knockout TF ,if TF=None it represents self-iteration
    # grn : The dictionary form represents top TF that regulates target gene
    # countnumber : the number of iterations
    # gene_premodel : the folder location of pre-trained gene prediction models
    # core : The core of parallel computing
    # save : the folder location of save TF knock-out simulation matrix
    
    Returns: HVG expression matrix
    
    """

    start = time.time()
    HVGdata_pd = imp_data[HVG_list].copy()
    if TF is not None:
        TFdata_pd = pd.read_pickle(save+TF+"_TFdata"+str(countnumber)+".pkl")
    else :
        TFdata_pd = pd.read_pickle(save+"cycle_TFdata"+str(countnumber)+".pkl")
    pool = Pool(processes=core)
    a = list(pool.starmap(gene_predict, zip(HVG_list, repeat(gene_premodel), repeat(TFdata_pd), repeat(grn)) ))
    pool.close()
    pool.join()
    inputdata_pd = pd.DataFrame(np.array(a).T, index=HVGdata_pd.index, columns=HVGdata_pd.columns)
    inputdata_pd[inputdata_pd<=0.0001]=0
    data_pd = pd.concat([TFdata_pd,inputdata_pd],axis=1)
    if TF is not None:
        data_pd.to_pickle(save+TF+str(countnumber)+".pkl")
    else :
        data_pd.to_pickle(save+"cycle"+str(countnumber)+".pkl")
    print("simulate_TF =",TF,"time =",time.time()-start)
    
    
    
def TF_matrix_predict(imp_data,TF_list,TF,grn,countnumber,
                      gene_premodel,save,pattern,core=10,
                      matrix_err = 10000,min_matrix_err = 0.01,
                      err_save = None ,TF_save=None) : 
    """
    Simulation prediction with a flexible number of iterations
    
    Arguments:
    
    # imp_data : imputation data which represents normalized data.  Rows represent cells, columns represent genes
    # TF_list : List of TFs that have been screened
    # TF : the simulate knockout TF
    # grn : The dictionary form represents top TF that regulates target gene
    # countnumber : the max number of iterations
    # gene_premodel : the folder location of pre-trained gene prediction models
    # pattern : "transform" Determine the number of iterations based on matrix err percent, "perturb" gene_err
    # core : The core of parallel computing
    # matrix_err : the initial value of matrix mse
    # min_matrix_err : the min value of matrix mse
    # save : the folder location of save TF knock-out simulation matrix
    # err_save : the folder location of err
    # TF_save : the folder location of TF
    
    Returns: TF expression matrix
    
    """
    start = time.time()
    TFdata_pd = imp_data[TF_list].copy()
    count = 0
    gene_errlist = []
    tmp_TF_list=[]
    matrix_errlist = []
    matrix_err_percent_list = []
    TFdata_pd[TF]=0 
    while count< countnumber and matrix_err > min_matrix_err:
        pool = Pool(processes=core)
        a = list(pool.starmap(gene_predict, zip(TF_list, repeat(gene_premodel), repeat(TFdata_pd), repeat(grn)) ))
        pool.close()
        pool.join()
        inputdata_pd = pd.DataFrame(np.array(a).T, index=TFdata_pd.index, columns=TFdata_pd.columns)
        gene_err = np.sqrt(np.sum(inputdata_pd[TF].values**2))
        gene_errlist.append(gene_err)
        tmp_TF_list.append(inputdata_pd[TF])
        inputdata_pd[TF] = 0
        inputdata_pd[inputdata_pd<=0.0001]=0
        matrix_err = np.sqrt(np.sum((TFdata_pd.values - inputdata_pd.values)**2))
        matrix_err_percent = np.sqrt(np.sum((TFdata_pd.values - inputdata_pd.values)**2))/np.sqrt(np.sum((TFdata_pd.values)**2))
        matrix_errlist.append(matrix_err)
        matrix_err_percent_list.append(matrix_err_percent)
        TFdata_pd = inputdata_pd.copy()
        
        if pattern=="convergence":
            if count>3 and 
            (abs(gene_errlist[count-1]-gene_errlist[count])/gene_errlist[count-1])<0.05:
                count = count-2
  
        if pattern=="transform":
            if count ==0 and matrix_err_percent>0.15:
                count = countnumber-2
            if count ==0 and matrix_err_percent>0.1:
                count = countnumber-3
            if count ==0 and matrix_err_percent>0.05:
                count = countnumber-5   
                
        if pattern=="perturb":   
            if count>1:
                if gene_errlist[count]>=gene_errlist[count-1]:
                    count = countnumber-1
   
            
        count = count + 1
        if (count == countnumber) or (matrix_err <= min_matrix_err) :
            TFdata_pd.to_pickle(save+TF+"_TFdata"+".pkl")
        print("counts =",count,"    ","simulate_TF =",TF,"    ","matrix_err =",matrix_err,"    ",
              "matrix_err_percent =",matrix_err_percent, "    ","gene_err =",gene_err,
              "time =",time.time()-start)
    if err_save is not None :
        pd.DataFrame({"gene_err":gene_errlist,
              "matrix_err":matrix_errlist,
              "matrix_err_percent":matrix_err_percent_list}).to_pickle(save+"err/"+TF+".pkl")
    if TF_save is not None : 
        pd.DataFrame(np.array(tmp_TF_list).T, index=TFdata_pd.index).to_pickle(save+"gene/"+TF+".pkl")
        
        
        
def HVG_matrix_predict(imp_data,HVG_list,TF,grn,
                      gene_premodel,save,core=10) : 
    
    """
    Simulation prediction with a flexible number of iterations
    
    Arguments:
    
    # imp_data : imputation data which represents normalized data.  Rows represent cells, columns represent genes
    # HVG_list : List of HVG that have been screened
    # TF : the simulate knockout TF
    # grn : The dictionary form represents top TF that regulates target gene
    # gene_premodel : the folder location of pre-trained gene prediction models
    # core : The core of parallel computing
    # save : the folder location of save TF knock-out simulation matrix
    
    Returns: HVG expression matrix
    
    """

    start = time.time()
    HVGdata_pd = imp_data[HVG_list].copy()
    TFdata_pd = pd.read_pickle(save+TF+"_TFdata"+".pkl")
    pool = Pool(processes=core)
    a = list(pool.starmap(gene_predict, zip(HVG_list, repeat(gene_premodel), repeat(TFdata_pd), repeat(grn)) ))
    pool.close()
    pool.join()
    inputdata_pd = pd.DataFrame(np.array(a).T, index=HVGdata_pd.index, columns=HVGdata_pd.columns)
    inputdata_pd[inputdata_pd<=0.0001]=0
    data_pd = pd.concat([TFdata_pd,inputdata_pd],axis=1)
    data_pd.to_pickle(save+TF+".pkl")
    print("simulate_TF =",TF,"time =",time.time()-start)
    
    