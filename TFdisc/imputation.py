import numpy as np
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn import preprocessing
import os
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2.robjects.packages import importr


def imp_SAVER(data, cores=20,normalization=True):
    pandas2ri.activate()
    rscript = """
    process_saver <- function(data,cores,normalization) {
        if(require(SAVER)){
            print("load SAVER")
            print(packageVersion("SAVER"))
        }else{ 
            install.packages("SAVER")
            if(require(SAVER)){
                print("install and load SAVER")
                print(packageVersion("SAVER"))
            } else {
                stop("fail to install and load SAVER")
            }
            }
        if(require(SAVER)){
            if(normalization){
                impute_data = saver(data,ncores = cores,estimates.only = TRUE,size.factor = 1)
                return(impute_data)
            }else{
                impute_data = saver(data,ncores = cores,estimates.only = TRUE)
                return(impute_data)
            }
        }
    }
    """
    r(rscript)
    
    impute_data = r['process_saver'](pandas2ri.py2rpy(data.T),cores,normalization)
    impute_data = pd.DataFrame(impute_data.T, columns=data.columns, index=data.index)
    return impute_data

def imp_Dca(data, cores=20):
    try:
        import scanpy
    except ImportError:
        raise ImportError("Please install scanpy package")
        
    try:
        from dca.api import dca
    except ImportError:
        raise ImportError("Please install dca package")
    adata = sc.AnnData(data)    
    dca(adata, threads=cores)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    return pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

def imp_Magic(data, cores=20):
    try:
        import scanpy
    except ImportError:
        raise ImportError("Please install scanpy package")
        
    try:
        import magic
    except ImportError:
        raise ImportError("Please install magic package")
    magic_operator = magic.MAGIC(n_jobs=cores)
    X_magic = magic_operator.fit_transform(data,genes=data.columns)
    return X_magic


def imp_Knn(data, n_neighbors=10, use_rep=None, key_added='knn_average'):
    try:
        import scanpy
    except ImportError:
        raise ImportError("Please install scanpy package")
    adata = sc.AnnData(data)  
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)
    knn_graph = adata.obsp['connectivities']
    if isinstance(adata.layers["counts"], np.ndarray):
        data_matrix = adata.layers["counts"]
    else:
        data_matrix = adata.layers["counts"].toarray()

    knn_avg_matrix = knn_graph @ data_matrix 
    knn_avg_matrix /= knn_graph.sum(axis=1)  

    adata.obsm[key_added] = knn_avg_matrix
    
    imp_adata = sc.AnnData(adata.obsm[key_added])  
    sc.pp.normalize_total(imp_adata)
    sc.pp.log1p(imp_adata)
    return pd.DataFrame(imp_adata.X, index=adata.obs_names, columns=adata.var_names)

