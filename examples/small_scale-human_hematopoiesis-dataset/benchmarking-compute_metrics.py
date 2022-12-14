import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import scipy as sci
import sklearn as sk
import scib
from scipy.sparse.csgraph import connected_components
import re
import os


def nmi_batch(adata, cluster_key, label_key, batch_key):
  nmi = []
  for label in adata.obs[label_key].cat.categories:
    adata_sub = adata[adata.obs[label_key].isin([label])]
    nmi.append(1 - scib.metrics.nmi(adata=adata_sub, group1=cluster_key, group2=batch_key))
  
  return np.mean(nmi)


os.chdir('~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/')
md = pd.read_csv('sample-info.tsv', delimiter='\t', header=0, index_col=0)
clusters = pd.read_csv(
  'benchmarking/Clusters.csv',
  header=0, index_col=0, dtype='category'
)

emb_names = clusters.columns
eval_metrics = pd.DataFrame(
  columns=['HG', 'ASW', 'NMI', 'ARI', 'bASW', 'bNMI', 'GC', 'kBET'],
  index=emb_names
)
for emb_name in emb_names:
  emb = pd.read_csv('embeddings/' + emb_name + '.csv', header=0, index_col=0)
  obs = clusters.loc[:, emb_name].to_frame()
  obs['Cell_type'] = pd.Categorical(md.loc[obs.index, 'Cell_type'])
  obs['Donor'] = pd.Categorical(md.loc[obs.index, 'Donor'])
  obj = ad.AnnData(obs=obs, obsm={emb_name: emb.loc[obs.index].values})
  
  k = 20
  metric = 'cosine' if re.search('CellSpace', emb_name) else 'euclidean'
  dim = obj.obsm[emb_name].shape[1]
  sc.pp.neighbors(
    adata=obj, use_rep=emb_name, n_pcs=dim,
    metric=metric, n_neighbors=k, knn=True
  )
  
  cluster_key = emb_name; label_key = 'Cell_type'; batch_key = 'Donor'
  eval_metrics.loc[emb_name, 'HG'] = sk.metrics.homogeneity_score(
    labels_true=obj.obs[label_key],
    labels_pred=obj.obs[cluster_key]
  )
  eval_metrics.loc[emb_name, 'ASW'] = scib.metrics.silhouette(
    adata=obj, group_key=label_key,
    embed=emb_name, metric=metric
  )
  eval_metrics.loc[emb_name, 'NMI'] = scib.metrics.nmi(adata=obj, group1=cluster_key, group2=label_key)
  eval_metrics.loc[emb_name, 'ARI'] = scib.metrics.ari(adata=obj, group1=cluster_key, group2=label_key)
  eval_metrics.loc[emb_name, 'bASW'] = scib.metrics.silhouette_batch(
    adata=obj, batch_key=batch_key, group_key=label_key,
    embed=emb_name, metric=metric,
    return_all=False, verbose=False
  )
  eval_metrics.loc[emb_name, 'bNMI'] = nmi_batch(
    adata=obj, cluster_key=cluster_key,
    label_key=label_key, batch_key=batch_key
  )
  eval_metrics.loc[emb_name, 'GC'] = scib.metrics.graph_connectivity(adata=obj, label_key=label_key)
  eval_metrics.loc[emb_name, 'kBET'] = scib.metrics.kBET(
    adata=obj, batch_key=batch_key, label_key=label_key,
    embed=emb_name, type_='knn',
    return_df=False, verbose=False
  )

eval_metrics.to_csv('benchmarking/metrics-summary.csv')

