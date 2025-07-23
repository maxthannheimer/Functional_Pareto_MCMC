#inputs for simulating observation data



gridsize=9#lenght of fine grid
N_fine=gridsize^2 #number of fine grid points
N_coarse=8 #number of weather stations/ conditioning points, obersavation points
num_sim=200 #number of simulated realizations
num_runs=30000 #use mcmc approach via simu_specfcts_MCMC, let chain run for num_runs steps for each observation

#true params for simulation
beta_true=0.5
c_true=3.0
alpha_true = 2.0

#threshold and threshold method
threshold_method="fixed"
threshold=1.0

#inputs for MCMC
N_MCMC=10000
N_est_c=20000
N_est_c_approx=20000
N_cond_sim=100
#N_burn_in=1000
n_trial_print=2000 #print every n-th step



#Test:
#=
gridsize=7#lenght of fine grid
N_fine=gridsize^2 #number of fine grid points
N_coarse=8 #number of weather stations/ conditioning points, obersavation points
num_sim=200 #number of simulated realizations
num_runs=100 #use mcmc approach via simu_specfcts_MCMC, let chain run for num_runs steps for each observation

#true params for simulation
beta_true=1.5
c_true=2.0
alpha_true = 2.0

#threshold and threshold method
threshold_method="fixed"
threshold=1.0

#inputs for MCMC
N_MCMC=100
N_est_c=200
N_est_c_approx=100
N_cond_sim=10
#N_burn_in=0
n_trial_print=500 #print every n-th step
=#