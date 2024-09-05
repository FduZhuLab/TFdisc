import os
import glob
import pickle
import pandas as pd
import numpy as np
import time
from arboreto.algo import genie3
from arboreto.algo import grnboost2
from dask.diagnostics import ProgressBar
from distributed import LocalCluster, Client
import warnings
warnings.filterwarnings("ignore")


def TF_grn(G_data,
           tf ,
           genelist,
           method = "grnboost2",
           n_core = 16,
           n_top = 50,
           verbose = False):
    """
    Calculate gene regulatory network(GRNs) using the scRNA-seq data
    
    Arguments:
    G_data: a pandas DataFrame of gene expression data, row(cell), col(gene).
    tf: list of transcription factors. If None or 'all', the list of gene names in expression_data will be used.
    genelist: list of genes used to be target, including TFs and HVGs. If None or 'all', the list of gene names in expression_data will be used.
    method: Method of constructing gene regulatory network. "grnboost2"(fast) or "genie3"(slow)
    n_core: The number of cores to select
    n_top: The number of top regulators for target gene
    verbose: print info.
    
    Returns: a dictionary representing the inferred gene regulatory links. key: target gene. value: top regulators
    """
    
    if isinstance(G_data, pd.DataFrame):
        gene = list(G_data.columns)
    if tf is None:
        tf = gene
    elif tf == 'all':
        tf = gene
    else:
        if len(tf) == 0:
            raise ValueError('tf names is empty')
        if not set(gene).intersection(set(tf)):
            raise ValueError('Intersection of data.columns and tf names is empty.')
        if  set(gene).intersection(set(tf)):
            tf = list(set(gene) & set(tf))
    if genelist is None:
        genelist = gene
    elif genelist == 'all':
        genelist = gene
    else:
        if len(genelist) == 0:
            raise ValueError('genelist is empty')
        if not set(gene).intersection(set(genelist)):
            raise ValueError('Intersection of data.columns and genelist is empty.')
        if  set(gene).intersection(set(genelist)):
            genelist = list(set(gene) & set(genelist))
    tf = list(set(tf) & set(genelist))
    start = time.time()
    local_cluster = LocalCluster(n_workers=n_core,threads_per_worker=1)
    custom_client = Client(local_cluster)
    
    if method == "grnboost2":
        adjancencies = grnboost2(expression_data=G_data[genelist], tf_names=tf, verbose=verbose,client_or_address=custom_client)
    if method == "genie3":
        adjancencies = genie3(expression_data=G_data[genelist], tf_names=tf, verbose=verbose,client_or_address=custom_client)
    adjancencies.sort_values(by='importance', ascending=False,inplace=True)
    groupby_obj=adjancencies.groupby(by="target")
    top_select = []
    name = []
    for i in groupby_obj:
        top_select.append(list(i[1].TF[0:n_top]))
        name.append(i[0])
    fina_grn = dict(zip(name,top_select))
    print("running time : ",time.time()-start)
    return fina_grn

