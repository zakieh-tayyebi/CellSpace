import palantir as pl
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn import preprocessing
from scipy.sparse import csr_matrix, find
from scipy.sparse.linalg import eigs
import os

# from palantir:
def run_diffusion_maps(
  data_df,
  n_components=10,
  metric='cosine',
  n_neighbors=30,
  alpha=0
):
    # Determine the kernel
    print('Determing nearest neighbor graph...')
    N = data_df.shape[0]
    temp = sc.AnnData(data_df.values)
    sc.pp.neighbors(temp, n_pcs=0, n_neighbors=n_neighbors, metric=metric)
    kNN = temp.obsp['distances']

    # Adaptive k
    adaptive_k = int(np.floor(n_neighbors / 3))
    adaptive_std = np.zeros(N)

    for i in np.arange(len(adaptive_std)):
        adaptive_std[i] = np.sort(kNN.data[kNN.indptr[i] : kNN.indptr[i + 1]])[
            adaptive_k - 1
        ]

    # Kernel
    x, y, dists = find(kNN)

    # X, y specific stds
    dists = dists / adaptive_std[x]
    W = csr_matrix((np.exp(-dists), (x, y)), shape=[N, N])

    # Diffusion components
    kernel = W + W.T

    # Markov
    D = np.ravel(kernel.sum(axis=1))

    if alpha > 0:
        # L_alpha
        D[D != 0] = D[D != 0] ** (-alpha)
        mat = csr_matrix((D, (range(N), range(N))), shape=[N, N])
        kernel = mat.dot(kernel).dot(mat)
        D = np.ravel(kernel.sum(axis=1))

    D[D != 0] = 1 / D[D != 0]
    T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(kernel)
    # Eigen value dcomposition
    D, V = eigs(T, n_components, tol=1e-4, maxiter=1000)
    D = np.real(D)
    V = np.real(V)
    inds = np.argsort(D)[::-1]
    D = D[inds]
    V = V[:, inds]

    # Normalize
    for i in range(V.shape[1]):
        V[:, i] = V[:, i] / np.linalg.norm(V[:, i])

    # Create are results dictionary
    res = {'T': T, 'EigenVectors': V, 'EigenValues': D}
    res['EigenVectors'] = pd.DataFrame(res['EigenVectors'])
    res['EigenVectors'].index = data_df.index
    res['EigenValues'] = pd.Series(res['EigenValues'])
    res['kernel'] = kernel

    return res

# from palantir:
def boundary_cells(ms_data, scale_components=True):
    if scale_components:
        data = pd.DataFrame(preprocessing.minmax_scale(ms_data),
                            index=ms_data.index, columns=ms_data.columns)
    else: data = copy.copy(ms_data)
    return pd.Index(set(data.idxmax()).union(data.idxmin()))


os.chdir('~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/')
cs_emb = pd.read_csv('embeddings/CellSpace-var_tiles.csv', header=0, index_col=0)

dm_res = run_diffusion_maps(cs_emb, n_components=15, metric='cosine', n_neighbors=20)
ms_data = pl.utils.determine_multiscale_space(dm_res, n_eigs=11)

boundary_cells(ms_data)
start = 'SRR5353471' # HSC boundary cell

pr_res = pl.core.run_palantir(ms_data, knn=20, early_cell=start, use_early_cell_as_start=True)
pr_res.branch_probs.to_csv('CellSpace-results/var_tiles/palantir/branch_probs.csv')
pr_res.entropy.to_csv('CellSpace-results/var_tiles/palantir/entropy.csv')
pr_res.pseudotime.to_csv('CellSpace-results/var_tiles/palantir/pseudotime.csv')

