import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

###
numFiles = 20

N_ave_0 = np.zeros((30,30,numFiles))
Z_ave_0 = np.zeros((30,30,numFiles))

for i in np.arange(0,numFiles):
    N = np.load("N_mat_self_rec0_"+str(i)+".npy")
    Z = np.load("Z_mat_self_rec0_"+str(i)+".npy")
    
    N_ave_0[:,:,i]=N[:,0,-1,:,:].mean(axis=0)
    Z_ave_0[:,:,i]=N[:,0,-1,:,:].mean(axis=0)

np.save("N_ave_0.npy", N_ave_0)
np.save("Z_ave_0.npy", Z_ave_0)

###
numFiles = 20

N_ave_1 = np.zeros((30,30,numFiles))
Z_ave_1 = np.zeros((30,30,numFiles))

for i in np.arange(0,numFiles):
    N = np.load("N_mat_self_rec_"+str(i)+".npy")
    Z = np.load("Z_mat_self_rec_"+str(i)+".npy")
    
    N_ave_1[:,:,i]=N[:,0,-1,:,:].mean(axis=0)
    Z_ave_1[:,:,i]=N[:,0,-1,:,:].mean(axis=0)

np.save("N_ave_1.npy", N_ave_1)
np.save("Z_ave_1.npy", Z_ave_1)
