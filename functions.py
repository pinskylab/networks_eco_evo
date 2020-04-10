#=========================================================================================
# Functions for Eco-evo network model
# Note that while the paper is ultimately about generic species (and allows for a second,
# competitor species), the script was originally based on a coral species with an optional
# algae competitor. 
# This work also allows designation of certain patches as 'marine protected areas,' where
# mortality of species 2 (algae) is elevated. This was not used in our simulations. 
#=========================================================================================

#=========================================================================================
# Load required packages
#=========================================================================================
import numpy as np

#=========================================================================================
# Growth rate function: growth_fun()
# Calculates growth rate as a function of species trait (optimum temperature) and local 
# patch temperature.
# Parameters:
#   r_max: Maximum growth rate for species
#   T: Current temperature of reef
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   species_type: Species type ID that determines type of growth
#=========================================================================================
def growth_fun(r_max,T,z,w,species_type):
    if z.shape[0] > 1: #If there is more than one patch
        T = np.repeat(T,z.shape[1]).reshape(z.shape[0],z.shape[1])
    else: #If there is a single patch
        T = np.array([np.repeat(T,z.shape[1])])
    r = np.zeros((z.shape[0],z.shape[1])) #Preallocate growth vector
    coral_col = np.where(species_type == 1)[1]
    algae_col = np.where(species_type == 2)[1]
    r[:,coral_col] =( (r_max[:,coral_col]/np.sqrt(2.*np.pi*pow(w[:,coral_col],2.)))
                    *np.exp((-pow((T[:,coral_col]-z[:,coral_col]),2.))/(2*pow(w[:,coral_col],2.))) )
    r[:,algae_col] = 0.49 * r_max[:,algae_col]
    return r

#=========================================================================================
# Mortality function: mortality_fun()
# Determines mortality rate as a function of species trait (optimum temperature) and local 
# environmental condition (reef temperature).
# Note: MPA status only affects algal mortality rate
# Parameters:
#   r_max: Maximum growth rate (potentially used for other mortality functional forms)
#   T: Current temperature of patch
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   species_type: Species type ID that determines type of mortality
#   mpa_status: Protected area status of the patch
#   alg_mort: algal mortality rate in non-mpa patch
#=========================================================================================
 def mortality_fun(r_max,T,z,w,species_type,mpa_status,alg_mort):
    #Preallocate mortality vector
    m = np.zeros((z.shape[0],z.shape[1])) 
    #Find algae columns
    algae_col = np.array([np.where(species_type == 2)[1]])
    
    if z.shape[0] > 1: # If there is more than one patch
        #Reshape temperature array
        T = np.repeat(T,z.shape[1]).reshape(z.shape[0],z.shape[1]) 
        m[z<T] = 1 - np.exp(-pow((T-z),2)/pow(w,2))[z<T]
        m[z>=T] = 0
        #Indices of mpa patches (corresponds to rows in N_all)
        is_mpa = np.array([np.where(mpa_status == 1)[1]]) 
        #Indices of non-mpa patches (corresponds to rows in N_all)
        not_mpa = np.array([np.where(mpa_status != 1)[1]])
        # Create arrays of indices that correspond to is_mpa & algae_col 
        # and not_mpa & algae_col
        is_mpa_rows = np.array([is_mpa.repeat(algae_col.shape[1])]) 
        not_mpa_rows = np.array([not_mpa.repeat(algae_col.shape[1])])
        algae_col_is_mpa = np.tile(algae_col,is_mpa.shape[1])
        algae_col_not_mpa = np.tile(algae_col,not_mpa.shape[1])
        # Macroalgae calculations for multiple reefs
        m[is_mpa_rows,algae_col_is_mpa] = 0.3
        m[not_mpa_rows,algae_col_not_mpa] = alg_mort[not_mpa_rows]
    
    else: # If there is a single patch
        T = np.array([np.repeat(T,z.shape[1])])
        #Coral calculations
        m[z<T] = 1 - np.exp(-pow((T-z),2)/pow(w,2))[z<T]
        m[z>=T] = 0
        # Macroalgae calculations for a single patch
        if mpa_status == 1:
             m[0,algae_col] = 0.3
        else:
            m[0,algae_col] = alg_mort
    
    # Apply a correction such that the minimum amount of mortality experienced is 0.03        
    m[m<0.03] = 0.03
    return m
    
#=========================================================================================
# Fitness function: fitness_fun()
# Determines fitness as a function of local growth rate, mortality rate, and species 
# interactions (competition for space).
# Calls:
#   growth_fun( )
#   mortality_fun( )
# Parameters:
#   r_max: Maximum growth rate (potentially used for other mortality functional forms)
#   T: Current temperature of reef
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   alphas: Species interaction matrix
#   species_type: Species type ID that determines type of mortality
#   mpa_status: Protected area status of the reef
#   N_all: Matrix of abundance for all species
#   m_const: Constant mortality rate
#   mortality_model: Use either constant or temperature-varying mortality be used? 
#       ("temp_vary" or "constant")
#   alg_mort: algal mortality rate in non-mpa reefs
#=========================================================================================
def fitness_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort):
    r = growth_fun(r_max,T,z,w,species_type)
     # If mortality varies with temperature
    if mortality_model == "temp_vary":
        m = mortality_fun(r_max,T,z,w,species_type,mpa_status,alg_mort)    
    # If mortality is constant
    else: 
        m = m_const
    #If there is more than one patch
    if N_all.shape[0] > 1:
        sum_interactions = np.array([np.sum(N_all[index,:] * alphas, axis=1) for index in range(N_all.shape[0])])
    else:
        sum_interactions = np.sum(alphas * N_all, axis=1)
    
    g = r * (1-sum_interactions) - m
    return g

#=========================================================================================
# Function to calculate dg/dz: dGdZ_fun()
# Calculates the partial derivative of growth rate across changes in trait 
# space.
# Calls:
#   fitness_fun()
# Parameters:
#   r_max: Maximum growth rate (potentially used for other mortality functional forms)
#   T: Current temperature of reef
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   alphas: Species interaction matrix
#   species_type: Species type ID that determines type of mortality
#   mpa_status: Protected area status of the reef
#   N_all: Vector of abundance for all species
#   m_const: Constant mortality rate
#   mortality_model: Use either constant or temperature-varying mortality be used? 
#       ("temp_vary" or "constant")
#   alg_mort: algal mortality rate in non-mpa reefs
#=========================================================================================        
def dGdZ_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort):
    h = 1e-5
    dGdZ = np.zeros(z.shape)
    
    #If there is more than one patch
    if N_all.shape[0] > 1:
        # For each patch
        for i in np.arange(z.shape[0]): 
            # For each species
            for j in np.arange(z.shape[1]):
                h_matrix = np.zeros(z.shape)
                h_matrix[i,j] = h
                # Take the symmetric difference quotient at point z[i,j]
                term1 = fitness_fun(r_max,T,z+h_matrix,w,alphas,species_type,mpa_status,
                                     N_all,m_const,mortality_model,alg_mort)
                term2 = fitness_fun(r_max,T,z-h_matrix,w,alphas,species_type,mpa_status,
                                         N_all,m_const,mortality_model,alg_mort)
                delta = (term1-term2)/(2*h)
                dGdZ[i,j] = delta[i,j]
    else:
        for j in np.arange(z.shape[1]):
            h_array = np.zeros(z.shape)
            h_array[0,j] = h
            # Take the symmetric difference quotient at point z[i,j]
            term1 = fitness_fun(r_max,T,z+h_array,w,alphas,species_type,mpa_status,
                                 N_all,m_const,mortality_model,alg_mort)
            term2 = fitness_fun(r_max,T,z-h_array,w,alphas,species_type,mpa_status,
                                     N_all,m_const,mortality_model,alg_mort)
            delta = (term1-term2)/(2*h)
            dGdZ[0,j] = delta[0,j]
            
    return dGdZ

#=========================================================================================
# Function to calculate d2g/dz2: dGdZ2_fun()
# Calcualtes the second partial derivative of growth rate across changes in 
# trait space.
# Calls: 
#   fitness_fun()
# Parameters:
#   r_max: Maximum growth rate (potentially used for other mortality functional forms)
#   T: Current temperature of reef
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   alphas: Species interaction matrix
#   species_type: Species type ID that determines type of mortality
#   mpa_status: Protected area status of the reef
#   N_all: Matrix of abundance for all species
#   m_const: Constant mortality rate
#   mortality_model: Use either constant or temperature-varying mortality be used? 
#       ("temp_vary" or "constant")
#   alg_mort: algal mortality rate in non-mpa reefs
#=========================================================================================   
def dGdZ2_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort):
    h = 1e-5
    dGdZ2 = np.zeros(z.shape)
    
    #If there is more than one patch
    if N_all.shape[0] > 1:
        # For each patch
        for i in np.arange(z.shape[0]): 
            # For each species
            for j in np.arange(z.shape[1]):
                h_matrix = np.zeros(z.shape)
                h_matrix[i,j] = h
                # Take the symmetric difference quotient at point z[i,j]
                term1 = fitness_fun(r_max,T,z+h_matrix,w,alphas,species_type,mpa_status,
                                     N_all,m_const,mortality_model,alg_mort)
                term2 = fitness_fun(r_max,T,z-h_matrix,w,alphas,species_type,mpa_status,
                                         N_all,m_const,mortality_model,alg_mort)
                term3 = fitness_fun(r_max,T,z,w,alphas,species_type,mpa_status,
                                         N_all,m_const,mortality_model,alg_mort)
                delta = (term1+term2-2*term3)/pow(h,2)
                dGdZ2[i,j] = delta[i,j]
    else:
        for j in np.arange(z.shape[1]):
            h_array = np.zeros(z.shape)
            h_array[0,j] = h
            # Take the symmetric difference quotient at point z[i,j]
            term1 = fitness_fun(r_max,T,z+h_array,w,alphas,species_type,mpa_status,
                                 N_all,m_const,mortality_model,alg_mort)
            term2 = fitness_fun(r_max,T,z-h_array,w,alphas,species_type,mpa_status,
                                     N_all,m_const,mortality_model,alg_mort)
            term3 = fitness_fun(r_max,T,z,w,alphas,species_type,mpa_status,
                                     N_all,m_const,mortality_model,alg_mort)
            delta = (term1+term2-2*term3)/pow(h,2)
            dGdZ2[0,j] = delta[0,j]
            
    return dGdZ2
    
#=========================================================================================
# Function to calculate dN/dt: dNdt_fun()
# Calculate the derivative of relative abundance across time. 
# Calls:
#   fitness_fun()
# Parameters:
#   r_max: Maximum growth rate (potentially used for other mortality functional forms)
#   T: Current temperature of reef
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   alphas: Species interaction matrix
#   species_type: Species type ID that determines type of mortality
#   mpa_status: Protected area status of the reef
#   N_all: Matrix of abundance for all species
#   m_const: Constant mortality rate
#   mortality_model: Use either constant or temperature-varying mortality be used? ("temp_vary" or "constant")
#   alg_mort: Algal mortality rate in non-mpa reefs
#   V: Genetic variation
#   D: Dispersal matrix (same for all species)
#   beta: Vector of combined dipsersal and effective fecundity rate
#=========================================================================================     
def dNdt_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort,V,D,beta):    
    if N_all.shape[0] > 1:
        V = np.tile(V, N_all.shape[0]).reshape(N_all.shape[0],N_all.shape[1])
    g = fitness_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort)
    dGdZ2 = dGdZ2_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort)      
    popdy = np.multiply(N_all,g)
    genload = 0.5 * V * dGdZ2 * N_all
    dispersal = beta * np.dot(D,N_all)
    free_space = 1 - N_all.sum(axis=1)
    larval_input = np.array([dispersal[index,:] * x for index, x in enumerate(free_space)])
    algae_ID = np.where(species_type==2)[1] #find algae columns
    larval_input[:,algae_ID] = 0
    dNdt = popdy + genload + larval_input
    #! Prevent NaN or population values below 1e-6 in output
    if np.isnan(dNdt).any():
        ID = np.where(np.isnan(dNdt))
        dNdt[ID] = 1e-6
    if (dNdt+N_all < 1e-6).any():
        ID = np.where(dNdt+N_all < 1e-6)
        dNdt[ID] = 1e-6

    return dNdt
    
#=========================================================================================
# Function to calculate q: q_fun()
# Prevents directional selection of virtually extinct populations and enhance numerical 
# stability. 
# Parameters:
#   N_all: Matrix of abundance for all species
#   N_min: Minimum value for percent cover at any given location
#=========================================================================================
def q_fun(N_all, N_min=1e-6):
    q = np.maximum(0, 1- N_min/(np.maximum(N_min,2*N_all)))
    return q
    
#=========================================================================================
# Function to calculate dZ/dt: dZdt_fun()
# Calculates the derivative of average trait value across time. 
# Calls:
#   fitness_fun()
#   dGdZ_fun()
#Parameters:
#   r_max: Maximum growth rate (potentially used for other mortality functional forms)
#   T: Current temperature of reef
#   z: Current optimum temperature for species
#   w: Temperature tolerance
#   alphas: Species interaction matrix
#   species_type: Species type ID that determines type of mortality
#   mpa_status: Protected area status of the reef
#   N_all: Matrix of abundance for all species
#   m_const: Constant mortality rate
#   mortality_model: Use either constant or temperature-varying mortality be used? 
#       ("temp_vary" or "constant")
#   alg_mort: Algal mortality rate in non-mpa reefs
#   V: Genetic variation
#   D: Dispersal matrix (same for all species)
#   beta: Vector of combined dipsersal and effective fecundity rate
#=========================================================================================    
def dZdt_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort,V,D,beta):
    if N_all.shape[0] > 1:
        V = np.tile(V, N_all.shape[0]).reshape(N_all.shape[0],N_all.shape[1])
    
    q = q_fun(N_all)
    dGdZ = dGdZ_fun(r_max,T,z,w,alphas,species_type,mpa_status,N_all,m_const,mortality_model,alg_mort)
    directional_selection = q * V * dGdZ
    gene_flow_term1 = (np.dot(D, N_all * z) / np.dot(D, N_all)) - z
    gene_flow_term2 = (beta * np.dot(D, N_all)) / (beta * np.dot(D, N_all) + N_all)
    free_space = 1 - N_all.sum(axis=1)
    gene_flow = np.array([(gene_flow_term1*gene_flow_term2)[index,:] * x for index, x in enumerate(free_space)])
    algae_ID = np.where(species_type==2)[1] #find algae columns
    gene_flow[:,algae_ID] = 0
    dZdt = directional_selection + gene_flow
    return dZdt
        
#=========================================================================================
# Function to generate initial temperatures: generate_temps_cos_fun()
# Generates temperatures across the ring network. 
# Parameters:
#   size: Number of patches in the network
#   min_SST: minimum temperature in the network
#   max_SST: maximum temperature in the network
#=========================================================================================  
def generate_temps_cos_fun(size,min_SST=20,max_SST=30):
    range_SST = max_SST - min_SST
    mid = range_SST/2
    mid_SST = (max_SST + min_SST)/2
    x_vals = np.linspace(0,2,size)
    cos_vals = mid * np.cos(np.pi*x_vals) + mid_SST
    return cos_vals

#! GENERATE INITIAL TRAITS FUNCTION         
def generate_traits_fun(nsp,size,temps,mid=25,temp_range=2.5,trait_scenario='perfect_adapt'):
    if trait_scenario == 'u_constant':
        a = np.linspace(mid-(temp_range/4), mid+(temp_range/4), nsp+2)[1:-nsp+2]
        b = np.repeat(a,size)
        traits = np.reshape(b,(nsp,size))
    if trait_scenario == 'same_constant':
        traits = np.full((nsp,size),mid)
    elif trait_scenario == 'perfect_adapt':
        a = np.tile(temps,nsp)
        traits = np.reshape(a,(nsp,size))
    return traits.T
    
#! GENERATE INITIAL STATE FUNCTION     
def generate_state_fun(size, nsp, cover=0.01,random=False):
    state = np.full((size,nsp),cover)
    if random:
        state = np.full(nsp*size,np.random.uniform(1e-6,.33,nsp*size)).reshape(size,nsp)
    return state
    
#! SET MPA FUNCTION        
def set_MPA_fun(temps,N_all,species_type,size,amount=0.2,strategy='random'):
    mpa = np.zeros(size)
    ID = np.where(species_type==1)[1]
    corals = N_all[:,ID].sum(axis=1) # Get the subset of N_all that corresponds to coral cover per patch
    #The following line applies to the 'portfolio' strategy that I haven't coded here
    ncoral = np.asarray(np.where(species_type==1)).sum() # How many coral species?
    
    if strategy == 'none':
        index = np.asarray([])
    elif strategy=='hot':
        index = (-temps).argsort()[0:np.int(amount*size)]
    elif strategy=='cold':
        index = temps.argsort()[0:np.int(amount*size)]
    elif strategy=='hotcold':
        index = np.r_[0:np.int(amount*size/2),np.int(size-(amount*size/2)):np.int(size)]
    elif strategy=='space':
        index = np.round(np.linspace(0,size-1,np.int(amount*size)))
    elif strategy=='highcoral':
        index = (-corals).argsort()[0:np.int(amount*size)]
    elif strategy=='lowcoral':
        index = corals.argsort()[0:np.int(amount*size)]
    elif strategy=='random':
        index = np.random.choice(np.arange(0,size),np.int(amount*size), replace=False)
    
    mpa[index.astype(int)]=1
    return np.array([mpa])    
    
#! MAIN ROUTINE       
#! BURNIN + RUNTIME
def coral_trait_stoch_fun0(param,spp_state,trait_state,temps,anomalies,algaemort_full,temp_change="constant"):
   
    nsp = param['nsp']
    size = param['size']
    num_years = param['num_years']
    step_size = param['step_size']
    time_steps = param['time_steps']
    species_type = param['species_type']
    r_max = param['r_max']
    V = param['V']
    D = param['D']
    beta = param['beta']
    m_const = param['m_const']
    w = param['w']
    alphas = param['alphas']
    mpa_status = param['mpa_status']
    mortality_model = param['mortality_model']
    maxtemp = param['maxtemp']
    annual_temp_change = param['annual_temp_change']
    timemod = param['timemod']

    SST_matrix = np.zeros([size,num_years])
    SST_matrix[:,0] = temps    
        
    for i in np.arange(1, num_years):
        if temp_change == "sigmoid":
            temps =  SST_matrix[:,i-1] 
            dtemp1 = annual_temp_change*temps.mean()
            dtemp2 = 1-(temps.mean()/maxtemp)
            dtemp = np.repeat(dtemp1*dtemp2, size)
            SST_matrix[:,i] = SST_matrix[:,i-1] + dtemp
            
        if temp_change == "constant":
            dtemp = np.repeat(0, size)
            SST_matrix[:,i] = SST_matrix[:,i-1] + dtemp
            
        if temp_change == "linear":
            dtemp = np.repeat(annual_temp_change, size)
            SST_matrix[:,i] = SST_matrix[:,i-1] + dtemp
        
    #if temp_change == "sigmoid":
    SST_matrix =  SST_matrix + anomalies

    N_ALL = np.zeros((size,nsp,time_steps))
    Z_ALL = np.zeros((size,nsp,time_steps))
    N_ALL[:,:,0] = spp_state
    Z_ALL[:,:,0] = trait_state

    algaemort_sub = algaemort_full[:,timemod:timemod+num_years]
    tick=0
    # Second-order Runge Kutta solver
    for i in np.arange(0,num_years-1):
        alg_mort = algaemort_sub[:,i]
        for j in np.arange(0,len(np.arange(0,1,step_size))):    
            dN1 = dNdt_fun(r_max,SST_matrix[:,i],Z_ALL[:,:,tick],w,alphas,species_type,mpa_status,
                             N_ALL[:,:,tick],m_const,mortality_model,alg_mort,V,D,beta)
            dZ1 = dZdt_fun(r_max,SST_matrix[:,i],Z_ALL[:,:,tick],w,alphas,species_type,mpa_status,
                            N_ALL[:,:,tick],m_const,mortality_model,alg_mort,V,D,beta)

            N_ALL_1 = N_ALL[:,:,tick] + dN1
            Z_ALL_1 = Z_ALL[:,:,tick]  + dZ1

            dN2 = dNdt_fun(r_max,SST_matrix[:,i],Z_ALL_1,w,alphas,species_type,mpa_status,
                             N_ALL_1,m_const,mortality_model,alg_mort,V,D,beta)
            dZ2 = dZdt_fun(r_max,SST_matrix[:,i],Z_ALL_1,w,alphas,species_type,mpa_status,
                            N_ALL_1,m_const,mortality_model,alg_mort,V,D,beta)

            N_ALL[:,:,tick+1] = N_ALL[:,:,tick] + step_size*(dN1 + dN2)/2
            Z_ALL[:,:,tick+1] = Z_ALL[:,:,tick] + step_size*(dZ1 + dZ2)/2

            tick += 1
    
    return N_ALL, Z_ALL, SST_matrix