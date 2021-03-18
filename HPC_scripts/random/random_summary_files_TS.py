import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

###
numFiles = 20

N_ts_0 = np.zeros((20,3,1500,15,15,numFiles))
Z_ts_0 = np.zeros((20,3,1500,15,15,numFiles))

for i in np.arange(0,numFiles):
    N = np.load("N_mat_self_rec0_"+str(i)+".npy")
    Z = np.load("Z_mat_self_rec0_"+str(i)+".npy")
    
    N_ts_0[:,:,:,:,:,i]=N[:,:,:,:,:]
    Z_ts_0[:,:,:,:,:,i]=Z[:,:,:,:,:]
    
N_batch_ave_ts_0 = N_ts_0.mean(axis=5)
N_batch_std_ts_0 = N_ts_0.std(axis=5)
Z_batch_ave_ts_0 = Z_ts_0.mean(axis=5)
Z_batch_std_ts_0 = Z_ts_0.std(axis=5)

np.save("./output/N_rand_ave_ts_0.npy", N_batch_ave_ts_0)
np.save("./output/N_rand_std_ts_0.npy", N_batch_std_ts_0)
np.save("./output/Z_rand_ave_ts_0.npy", Z_batch_ave_ts_0)
np.save("./output/Z_rand_std_ts_0.npy", Z_batch_std_ts_0)

###
N_ts_1 = np.zeros((20,3,500,15,15,numFiles))
Z_ts_1 = np.zeros((20,3,500,15,15,numFiles))

for i in np.arange(0,numFiles):
    N = np.load("N_mat_self_rec_"+str(i)+".npy")
    Z = np.load("Z_mat_self_rec_"+str(i)+".npy")
    
    N_ts_1[:,:,:,:,:,i]=N[:,:,:,:,:]
    Z_ts_1[:,:,:,:,:,i]=Z[:,:,:,:,:]
    
N_batch_ave_ts_1 = N_ts_1.mean(axis=5)
N_batch_std_ts_1 = N_ts_1.std(axis=5)
Z_batch_ave_ts_1 = Z_ts_1.mean(axis=5)
Z_batch_std_ts_1 = Z_ts_1.std(axis=5)

np.save("./output/N_rand_ave_ts_1.npy", N_batch_ave_ts_1)
np.save("./output/N_rand_std_ts_1.npy", N_batch_std_ts_1)
np.save("./output/Z_rand_ave_ts_1.npy", Z_batch_ave_ts_1)
np.save("./output/Z_rand_std_ts_1.npy", Z_batch_std_ts_1)