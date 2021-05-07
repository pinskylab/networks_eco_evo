#! Parameters for Coral Adaptation model
import numpy as np

iterations = 20 #number of iterations 

nsp = 1 # Number of species in model
size = 20 # Number of reefs in model, size is 60 in Walsworth et al. model
burnin = 1500 # Length of burn-in period
runtime = 500 # Length of environmental change period
mid = 27. # Mean temperature across all reefs at start of simulation
temp_range = 3. # Range of temperatures across reefs at start of simulation
species_type = np.array([[1.]]) # Species type ID
species = ["C1"] # Species labels
temp_stoch = 0.6
r_max = np.array([[1.0]])
w = np.array([[1.5]])
beta = np.array([[0.05]])
m_const = 0.1 # Value used for constant mortality case
mortality_model = "temp_vary"
V = np.array([[0.]])
annual_temp_change = 0.01
maxtemp = 32 #for sigmoid temperature increase scenario
step_size = 1

#! Species interaction matrix
#species have no effect on each other or themselves
#size needs to be nsp x nsp
alphas = np.array([[1.]]) 

#! Create multivariate normal covariance matrix for temperature anomalies
mdim = size 
lindec = np.exp(np.linspace(0,-5,num=mdim)) # Decrease the second value to reduce the range of correlation
ma = np.zeros((mdim,mdim))
ma[:,0] = lindec
for i in np.arange(1,mdim):
    ma[i:,i] = lindec[0:-i]
ma_temp = ma + ma.T
np.fill_diagonal(ma_temp,np.diagonal(ma))
ma = ma_temp    
sds = np.array([np.repeat(.2*temp_stoch,size)])
b = np.multiply(sds,sds.T)
spatial_temp = b*ma

#! CONNECTIVITY MATRIX FUNCTIONS
#! Function that generates a normalized random matrix such that all columns (i.e. source sites) sum to 1
def D_random(size):
    D0 = np.random.rand(size,size)
    D_rand  = np.zeros((size,size)) # Preallocate matrix
    for i in np.arange(size):
        D_rand[:,i]  = D0[:,i] / D0[:,i].sum()
    return D_rand
    
#! Function that modifies a connectivity matrix to "bias" the diagonal
def bias_diagonal(A, alpha):
    I = np.identity(A.shape[0])
    A_bias = (1-alpha)*A + alpha*I
    return A_bias
   
#! Function that generates a matrix that approximates a linear configuration such that 
#! each patch can only send larvae to its two adjacent patches    
def D_linear(size, p_dispersal):
    D_linear = np.zeros((size,size))
    for i in np.arange(0,size):
        D_linear[i,i]= 1-p_dispersal*2
        if i != 0 and i != size-1:
            D_linear[i,i-1] = p_dispersal
            D_linear[i,i+1] = p_dispersal
        if i != 0 and i == size-1:
            D_linear[i,i-1] = p_dispersal
        if i == 0 and i != size-1:
            D_linear[i,i+1] = p_dispersal
    return D_linear  