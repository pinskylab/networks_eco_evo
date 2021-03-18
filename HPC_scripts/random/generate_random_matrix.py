from __future__ import division
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import NoNorm
import random
from scipy.linalg import circulant
import networkx as nx

def _distance_matrix(L):
    Dmax = L//2
 
    D  = range(Dmax+1)
    D += D[-2+(L%2):0:-1]
 
    return circulant(D)/Dmax
    
def _pd(d, p0, beta):
    return beta*p0 + (d <= p0)*(1-beta)
    
def watts_strogatz(L, p0, beta, directed=False, rngseed=1):
    rng = np.random.RandomState(rngseed)
 
    d = _distance_matrix(L)
    p = _pd(d, p0, beta)
 
    if directed:
        A = 1*(rng.random_sample(p.shape) < p)
        np.fill_diagonal(A, 0)
    else:
        upper = np.triu_indices(L, 1)
 
        A          = np.zeros_like(p, dtype=int)
        A[upper]   = 1*(rng.rand(len(upper[0])) < p[upper])
        A.T[upper] = A[upper]
 
    return A

def preserve_self_loops(matrix0, edges):
    # matrix0 needs to have 1's in the diagonal for this to work
    matrix1=matrix0
    for i in np.arange(0,matrix0.shape[0]):
        deg_dist = np.count_nonzero(matrix0[i,:])
        if deg_dist > edges:
            for k in np.arange(0,deg_dist-edges):
                non_zeros = (np.where(matrix0[i,:]>0)[0]).tolist()
                non_zeros.remove(i)
                random_index = non_zeros[np.random.randint(0,len(non_zeros))]
                matrix1[i,random_index] = 0
    return matrix1

def D_norm(D0):
    size = D0.shape[0]
    D1  = np.zeros((size,size)) # Preallocate matrix
    for i in np.arange(size):
        D1[:,i]  = (D0[:,i] / D0[:,i].sum()).T
    return D1

def gen_rand_matrix(seed):    
    # Network size = 20
    #! Note: directed = False produces a symmetric matrix
    random.seed(seed)

    L = 20
    K = 4
    p0 = K/(L-1)
    regular = watts_strogatz(beta=0,directed=False,L=L,p0=p0)

    G = nx.Graph()
    G.add_nodes_from(np.arange(0,20))

    for i in np.arange(0,20):
        for j in np.arange(0,20):
            if regular[i,j] == 1:
                G.add_edge(i,j)
            
    reg_G = nx.to_numpy_matrix(G)
    di = np.diag_indices(20)

    random_G = nx.double_edge_swap(nx.from_numpy_matrix(reg_G), nswap=20, max_tries=1000)
    random_G2 = nx.to_numpy_matrix(random_G)

    di = np.diag_indices(20)
    random_G2[di]=1

    random_G2_revised = np.asarray(preserve_self_loops(random_G2,edges=5))
    
    random_G2_norm = D_norm(random_G2_revised)
    
    return random_G2_norm
