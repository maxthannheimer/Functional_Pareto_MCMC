include("functions.jl")

#inputs
gridsize=10#lenght of fine grid
N_fine=gridsize^2 #number of fine grid points
N_coarse=10 #number of weather stations/ conditioning points, obersavation points
num_sim=250 #number of simulated realizations

#true params for simulation
alpha_true = 4.0
beta_true=1.0
c_true=1.0
param=[ c_true , beta_true]
alpha=4.0

#Threshold definition as quantile
p=0.99
threshold= (1-p)^(-1/alpha)
#threshold=1.0

#MCMC params
N_MCMC=1000
param_start=[1.0,1.0]
alpha_start=1.0
N_est_c=1000
N_est_cond=5
N_burn_in=1000

#create grids
#create grid and observation points, also plots if last argument is true
Create_Grid_and_Observation(gridsize,N_coarse, true)

#create grid and observation points, also plots if last argument is true
#here the coarse observation points are already on the fine grid
#also safe coords 
(coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation_on_fine_grid(gridsize,N_coarse, true)

#gives the nearest fine gride coordinates for each coarse grid obsrvation 
#last argument is the coarse coordinates rounded to the nearest gridpoint coordinates
coord_cond_rows=get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
#redefine it to make sure its in the right order
coord_coarse=coord_fine[coord_cond_rows,:]
coord_x0=coord_fine[row_x0,:]

##############
#variograms and matrices
##############

#vario
#variogram vario of the process 
vario([1.0 , 1.0 ], [1/2,1/2]) 
vario([1.0,1.0],param)
#vec_vario
##variogram vario of the process for 2 d coordinates, written as  matrix with normalization in point coord_x0
vec_vario(param,coord_fine,coord_fine[row_x0,:])
vec_vario(param,coord_fine,[0.0, 0.0])
vec_vario([1/2,1/2],[0.0 1.0;1.0 1.0],[0.0,0.0] )

#generates covariance function depending on variogram with reference point 
cov_fun_vario([1/2,1/2],[0.1,0.1],[0.2,0.2],[0.0, 0.0])
cov_fun_vario(param,coord_fine[1,:],coord_fine[2,:],coord_fine[row_x0,:])
cov_fun_vario(param,coord_fine[1,:],coord_fine[2,:],[0.0, 0.0])

#Simulate data on grid
#simulate data on all points and reduce it to observation data (  coarse observations)
#use mcmc approach via simu_specfcts_MCMC, let chain run for 2000 steps for each observation
num_runs=2000
@time(sim_data= [simu_specfcts_MCMC(num_runs, alpha, coord_fine,param,row_x0 ) for i in 1:num_sim] )
sim_data=reduce(hcat,sim_data)' #just make vector of vectors a matrix (same below for observations)
#observation_data=reduce(hcat,[sim_data[i,N_fine+1:N_fine+N_coarse] for i in 1:num_sim])' #first argument is number of sim, second is coordinate

#get threshold empirically:
@time( quantile_threshold(num_runs,alpha,coord_fine,param,row_x0,0.8,1000))

#reduce to the coarse sites, treat as observation data
observation_data=reduce(hcat,[sim_data[i,coord_cond_rows] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
observation_x0=reduce(hcat,[sim_data[i,row_x0] for i in 1:num_sim])'



#NOW WE HAVE OBSERVATIONS: NEXT STEP: FUNCTIONS TO DO INFERENCE

#chol_mat (not used in the moment)
#calculates cholesky matrix without normalization for given coordinates and the distance based variogram with parameters
cholmat=chol_mat(coord_fine, x->vario(x,param))


#Cov Mat 
#calculate cov matrix for to matrices of coordinates and the variogramm given via param and normalized at coord_x0 (or [0, 0]) 
cov_mat_for_vectors(coord_coarse, coord_coarse, param, coord_x0)
cov_mat_for_vectors(coord_coarse, coord_coarse, param, [0.0,0.0])




# Circulant embedding testing will follow in seperate file
n_test=10000
@time begin
tmp=[FBM_simu_fast_vec_dependent(param,gridsize,num_sim) for i in 1:n_test]
println("dependend simulation")
end


@time begin
tmp=[FBM_simu_fast_vec(param,gridsize,num_sim) for i in 1:n_test]
println("independet simulation")
end



#gaussian process simulation

#r_gaussian_vec
#simulate numrep many 1/α [G(s)-G(x0)-γ(s-x0)] 
num_rep=11
r_gaussian_vec(coord_fine,param,row_x0,num_rep,alpha)

#simulate numrep many exp(1/α[G(s)-G(x0)-γ(s-x0)])
r_log_gaussian_vec(coord_fine,param,row_x0,num_rep,alpha) 
#simulate one rep of exp(1/α[G(s)-G(x0)-γ(s-x0)])
r_log_gaussian(coord_fine,param,row_x0,alpha) 


#conditional simulation of exp(1/α[G(s)-G(x0)-γ(s-x0)]) | observation_data(i)./observation_x0(i)
r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0,alpha) 
#conditional simulation of exp(1/α[G(s)-G(x0)-γ(s-x0)]) | observation_data./observation_x0 as vector
#first dim is number of simulations, second is number of wanted simulations, thrisd number of fine grid sites
r_cond_log_gaussian_vec(observation_data,observation_x0, coord_fine,coord_coarse,param,row_x0,num_rep,alpha)[1][1]



[r_gaussian_vec(coord_fine,param,row_x0,num_rep,alpha) for j in 1:size(observation_data,1)][1][1]


#prior: log of log gaussian and log gaussian:
x,mu,sigma=7.0,0.0,3.0
log_likehood_log_gauss_1d(x,mu,sigma)
log(log_gauss_1d(x,mu,sigma))


# Estimation of 1/c via samples
number_of_exceed = 200
plots1 = Vector{}(undef, 4)
i=1
for N_est_c in [100,1000,10000,40000]
println("time for $N_est_c:")
@time(tmp=[l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) for i in 1:100])
#pl=scatter(1:100,sort(tmp),label="l_2_fun",title="l_2_fun for $N_est_c samples")
#hline!(pl, [mean(tmp)],label="Mean of estimates for l_2_fun: $(mean(tmp))",linewidth=3)
plots1[i]=histogram(tmp,title="$N_est_c samples, mean: $(mean(tmp)) ")
i=i+1
end
plot(plots1..., layout = (2,2),size=(1000,1000))

number_of_exceed = 200
plots2 = Vector{}(undef, 4)
i=1
for N_est_c in [100,1000,10000,40000]
println("time for $N_est_c:")
@time(tmp=[l_2_fun_dependent(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) for i in 1:100])
#pl=scatter(1:100,sort(tmp),label="l_2_fun",title="l_2_fun for $N_est_c samples")
#hline!(pl, [mean(tmp)],label="Mean of estimates for l_2_fun: $(mean(tmp))",linewidth=3)
plots2[i]=histogram(tmp,title="$N_est_c samples, mean: $(mean(tmp)) ")
i=i+1
end
plot(plots2..., layout = (2,2),size=(1000,1000))


#EVERYTHING IS PUT TOGETHER IN MCMC function which estimates parameters via MCMC algorithm for observaed data


println(param_start)
n_trial_print=100 #print every n-th step
@time (result= MCMC(N_MCMC,observation_data,observation_x0,threshold, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_est_cond))



#mean of params with burn in period canceled out
[mean(result["beta"][N_burn_in+1:N_MCMC+1]),mean(result["c"][N_burn_in+1:N_MCMC+1]),mean(result["alpha"][N_burn_in+1:N_MCMC+1])]

#plot and hist for beta,zeta,lambda_1,lambda_2
histogram(result["beta"][N_burn_in+1:N_MCMC+1],title="Histogram for β, Burn in: $N_burn_in")
p_beta=scatter(1:N_MCMC+1,result["beta"],label="$N_MCMC samples of beta",title="Markov chain for β")
hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β: $(round((mean(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)

histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c, Burn in: $N_burn_in")
p_c=scatter(1:N_MCMC+1,result["c"],label="$N_MCMC samples of c",title="Markov chain for c")
hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c: $(round((mean(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)

histogram(result["alpha"][N_burn_in+1:N_MCMC+1],title="Histogram for α, Burn in: $N_burn_in")
p_alpha=scatter(1:N_MCMC+1,result["alpha"],label="$N_MCMC samples of α",title="Markov chain for α")
hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α: $(round((mean(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)


# with exp(c) 
histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
p_c=scatter(1:N_MCMC+1,exp.(result["c"]),label="$N_MCMC samples of c",title="Markov chain for c, Burn in: $N_burn_in")
hline!(p_c, [mean(exp.(result["c"][N_burn_in+1:N_MCMC+1]))],label="Mean estimate for c: $(round((mean(exp.(result["c"][N_burn_in+1:N_MCMC+1]))),digits=3))",linewidth=3)





#code to safe everything

import Pkg;Pkg.add("JLD2")
using JLD2
repetition=1
save_object("$(N_MCMC)_simulation_$(gridsize)_grid_0$repetition.jld2", result)
