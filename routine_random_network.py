import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import NoNorm
import time
import random
from generate_random_matrix import *
from functions_revised3 import *
import random_batch_param_single1 as P
import networkx as nx

length = 30
alpha_V_mat0 = np.zeros((length,length))
cor1_mat0 = np.zeros((length,length))
cor2_mat0 = np.zeros((length,length))

alpha_V_mat = np.zeros((length,length))
cor1_mat = np.zeros((length,length))
cor2_mat = np.zeros((length,length))

alpha_vals = np.linspace(0,1.,length)
V_vals = np.linspace(0,0.2,length)

N_0_all = np.zeros((P.size,P.nsp,P.burnin,length,length))
Z_0_all = np.zeros((P.size,P.nsp,P.burnin,length,length))
SST0_all= np.zeros((P.size,P.burnin,length,length))

N_1_all = np.zeros((P.size,P.nsp,P.runtime,length,length))
Z_1_all = np.zeros((P.size,P.nsp,P.runtime,length,length))
SST_all = np.zeros((P.size,P.runtime,length,length))

anomalies_burn = np.tile(np.random.normal(0,P.temp_stoch,P.burnin),P.size).reshape((P.size,P.burnin)) 
anomalies_run = np.tile(np.random.normal(0,P.temp_stoch,P.runtime),P.size).reshape((P.size,P.runtime)) 
algaemort_full = np.random.uniform(0.15,0.15,(P.runtime+P.burnin)*P.size).reshape((P.size,P.runtime+P.burnin))


#! This function submits multiple "sub_jobs" to the cluster:
def bound(njob, job_i, array_size):
    """ Returns the range of values that the job job_i out of njob jobs
    should process.

    Args:
      njob: total number of jobs in the array.
      job_i: index of the current job.
      array_size: the size of the array to split.

    Returns:
      start: the index of array to start.
      stop: the index of the array to stop.
      Note that array[start:stop] returns array[start]...array[stop-1]; that is,
      array[stop] is excluded from the range.
      
    """
    step = np.ceil(float(array_size)/float(njob))
    print 'njob=', njob, 'job_i= ', job_i, 'array_size= ', array_size, 'step= ', step
    start = int(job_i*step)
    stop = int(min((job_i+1)*step, array_size))
    print 'start= ', start, 'stop= ', stop
    return start, stop


if __name__ == '__main__':
    # get job_i from the command line argument
    arg = sys.argv
    njob = int(arg[1])
    job_i  = int(arg[2])
    
    #SET NUMBER OF ITERATIONS

    seed_global = np.arange(0,P.iterations)
    start, stop = bound(njob, job_i, seed_global.size)
    seed_values = seed_global[start:stop] 
    
    for seed in seed_values:
        np.random.seed(seed)
    
        anomalies_burn = np.tile(np.random.normal(0,P.temp_stoch,P.burnin),P.size).reshape((P.size,P.burnin)) 
        anomalies_run = np.tile(np.random.normal(0,P.temp_stoch,P.runtime),P.size).reshape((P.size,P.runtime)) 
        algaemort_full = np.random.uniform(0.15,0.15,(P.runtime+P.burnin)*P.size).reshape((P.size,P.runtime+P.burnin))

        G = nx.fast_gnp_random_graph(n=20,p=4/20,seed=seed)
        D0 = nx.to_numpy_matrix(G)
        di = np.diag_indices(20)
        D0[di] = 1
        D = D_norm(D0)

        for i in np.arange(0,length):
            V = np.repeat(V_vals[i],P.nsp)
            for j in np.arange(0,length):
                alpha = alpha_vals[j]
                D_bias = P.bias_diagonal(D, alpha)
    
                #burnin parameters
                # The following line is for the non-cosine temperature scenario
                #SST0 = generate_temps_fun(P.size,P.mid,P.temp_range,temp_scenario='random')
                SST0 = generate_temps_cos_fun(P.size,min_SST=20,max_SST=30)
                spp_state = generate_state_fun(P.size, P.nsp, cover=0.25,random=False)
                trait_state = generate_traits_fun(P.nsp,P.size,SST0,P.mid,P.temp_range,trait_scenario='perfect_adapt')
                mpa_status = set_MPA_fun(SST0,spp_state,P.species_type,P.size,amount=0,strategy='none')
    
                time_steps = P.burnin
                num_years = P.burnin
                timemod = 0 #to offset algae mortality index
                parameters_dict = {'nsp': P.nsp, 
                                    'size': P.size, 
                                    'num_years': num_years,
                                    'step_size': P.step_size,
                                    'time_steps': P.burnin, 
                                    'species_type': P.species_type, 
                                    'V': V, 
                                    'D': D_bias, 
                                    'beta': P.beta,
                                    'r_max': P.r_max,
                                    'alphas': P.alphas,
                                    'mortality_model': P.mortality_model,
                                    'mpa_status': mpa_status,
                                    'w': P.w,
                                    'm_const': P.m_const,
                                    'maxtemp': P.maxtemp,
                                    'annual_temp_change': P.annual_temp_change,
                                    'timemod':timemod
                                    }
                            
                N_0, Z_0, SST_final0 = coral_trait_stoch_fun0(parameters_dict,spp_state,trait_state,SST0,
                                                                   anomalies_burn,algaemort_full,temp_change="constant")

                N_0_all[:,:,:,i,j] = N_0
                Z_0_all[:,:,:,i,j] = Z_0
                SST0_all[:,:,i,j] = SST_final0
        
       
                #! RUNTIME
                SST0 = SST_final0[:,-1]
                spp_state = N_0[:,:,-1]
                trait_state = Z_0[:,:,-1]
                timemod = P.burnin #to offset algae mortality index
                time_steps = P.runtime
                num_years = P.runtime

                #! RUNTIME: no MPA
                mpa_status = set_MPA_fun(SST0,spp_state,P.species_type,P.size,amount=0,strategy='none')
                parameters_dict = {'nsp': P.nsp, 
                                    'size': P.size, 
                                    'num_years': num_years,
                                    'step_size': P.step_size,
                                    'time_steps': P.runtime, 
                                    'species_type': P.species_type, 
                                    'V': V, 
                                    'D': D_bias, 
                                    'beta': P.beta,
                                    'r_max': P.r_max,
                                    'alphas': P.alphas,
                                    'mortality_model': P.mortality_model,
                                    'mpa_status': mpa_status,
                                    'w': P.w,
                                    'm_const': P.m_const,
                                    'maxtemp': P.maxtemp,
                                    'annual_temp_change': P.annual_temp_change,
                                    'timemod':timemod
                                    }
        
                N_1, Z_1, SST_final1 = coral_trait_stoch_fun0(parameters_dict,N_0[:,:,-1],Z_0[:,:,-1],
                                                                      SST_final0[:,-1],anomalies_run,algaemort_full,
                                                                      temp_change="sigmoid")
                N_1_all[:,:,:,i,j] = N_1
                Z_1_all[:,:,:,i,j] = Z_1
                SST_all[:,:,i,j] = SST_final1

                np.save("./single1/N_mat_self_rec0_" + str(seed) + ".npy", N_0_all)
                np.save("./single1/Z_mat_self_rec0_" + str(seed) + ".npy", Z_0_all)
                np.save("./single1/SST_mat_self_rec0_" + str(seed) + ".npy", SST0_all)
        
                np.save("./single1/N_mat_self_rec_" + str(seed) + ".npy", N_1_all)
                np.save("./single1/Z_mat_self_rec_" + str(seed) + ".npy", Z_1_all)
                np.save("./single1/SST_mat_self_rec_" + str(seed) + ".npy", SST_all)