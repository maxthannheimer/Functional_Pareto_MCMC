using Plots
using Random, Distributions
include("functions.jl")


a=rand(10000)
var(a)
1/size(a,1)*sum([(x-mean(a))^2 for x in a])








###################################
#copy of simulated_evaluation.jl
##################################
include("functions.jl")
using Measures
N_burn_in=100

number_sim_int=5
number_parallel_int=18
sim_res_folder_name=readdir("results_dir")[number_parallel_int]

path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
(result_cond_sim,result_approx,input_MCMC_cond_sim,input_MCMC_approx)=load(path)["single_stored_object"]
(N_MCMC,observation_data,observation_x0,threshold,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim)=input_MCMC_cond_sim
#plot grid 
t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid Simulations",color=:red,alpha=0.4)
scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations",color=:purple)
scatter!([coord_fine[row_x0,1]],[coord_fine[row_x0,2]],label="",color=:purple)
title!("Observation and Simulation Points")
display(t)

(N_MCMC,observation_data,observation_x0,threshold, alpha_start, coord_x0,coord_coarse,param_start,n_trial_print,N_est_c_approx)=input_MCMC_approx
#plot grid 
t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid Simulations",color=:red,alpha=0.4)
scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations",color=:purple)
scatter!([coord_fine[row_x0,1]],[coord_fine[row_x0,2]],label="",color=:purple)
title!("Observation and Simulation Points")
display(t)


#einzelcheck
#for number_parallel_int in acceptable_simulations
    number_parallel_int=18
    sim_res_folder_name=readdir("results_dir")[number_parallel_int]
    number_sim_int=5
        path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
        (result_cond_sim,result_approx,input_MCMC_cond_sim,input_MCMC_approx)=load(path)["single_stored_object"]
    path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
 (N_MCMC,observation_data,observation_x0,threshold,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim)=input_MCMC_cond_sim


plot_interval=9400:9500
println("Minimum exceed number: $(minimum(result_cond_sim["Number of exceedance"])) and maximum exceed number: $(maximum(result_cond_sim["Number of exceedance"]))")
println("Minimum likelihood: $(minimum(result_cond_sim["log_likelihood"])) and maximum likelihood: $(maximum(result_cond_sim["log_likelihood"]))")
#end

scatter(plot_interval,result_cond_sim["Number of exceedance"][plot_interval])
scatter(plot_interval,result_cond_sim["log_likelihood"][plot_interval])
result_cond_sim["log_likelihood"][plot_interval]

scatter(plot_interval,result_cond_sim["alpha"][plot_interval],label="α",title="Markov chain for α")
scatter(plot_interval,result_cond_sim["c"][plot_interval],label="c",title="Markov chain for c")
scatter(plot_interval,result_cond_sim["beta"][plot_interval],label="β",title="Markov chain for β")
println("Minimum beta: $(minimum(result_cond_sim["beta"])) and maximum beta: $(maximum(result_cond_sim["beta"]))")



scatter(N_burn_in:N_MCMC,result_cond_sim["alpha"][N_burn_in+1:N_MCMC+1],label="α",title="Markov chain for α, Burn in: $N_burn_in")
scatter(N_burn_in:N_MCMC,result_cond_sim["c"][N_burn_in+1:N_MCMC+1],label="c",title="Markov chain for c, Burn in: $N_burn_in")
scatter(N_burn_in:N_MCMC,result_cond_sim["beta"][N_burn_in+1:N_MCMC+1],label="β",title="Markov chain for β, Burn in: $N_burn_in")

n=9408
N_cond_sim=100
N_est_c=1000
result_cond_sim["log_likelihood"][n]
param=[result_cond_sim["c"][n],result_cond_sim["beta"][n]]
alpha_par=result_cond_sim["alpha"][n]
  num_obs=size(observation_data,1)
(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha_par, coord_fine,coord_coarse,param,row_x0 )

(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, 2.0, coord_fine,coord_coarse,[3.0,0.5],row_x0 )
param=[100,0.1]
alpha_par=50

 N=N_est_c
    tmp = vcat(tmp,r_log_gaussian_vec_dependent(coord_fine,param,row_x0, N,alpha_par))
    r_W_alpha_sample =        [mean( tmp[i] )^(alpha_par) for i in 1:N]
    mean_r_W_alpha_sample = mean(r_W_alpha_sample)

    1/mean_r_W_alpha_sample*sqrt(1/N*var(r_W_alpha_sample))
    
    N=Int(round(1/(mean_r_W_alpha_sample)^2*var(r_W_alpha_sample)/ 0.05^2)+1)
    




#calculate new log likelihood
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha_par)
        #l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        @time( l2=l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha_par,N_est_c) * size(modified_observation,1)) #number of exceedances
        l2
        @time(l2_og=l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha_par,N_est_c))
        l2_og
        @time(l2_new=l_2_fun_dependent_no_number_of_exceed_new(coord_fine, param,row_x0,alpha_par,N_est_c) * size(modified_observation,1))  
        l2_new
        l3=l_3_fun(modified_observation_x0, alpha_par, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha_par,0.0,1.0,50.0)
        log_likelihood_new =sum([l1,l2,l3,prior]) 
        log_likelihood_new_og =sum([l1,l2_og,l3,prior]) 
        log_likelihood_new_new =sum([l1,l2_new,l3,prior])
param=[3.0,0.5]
alpha_par=2.0
  num_obs=size(observation_data,1)
(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha_par, coord_fine,coord_coarse,param,row_x0 )
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha_par) 
        l2=l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha_par,10000) * size(modified_observation,1) 
        l2_new=l_2_fun_dependent_no_number_of_exceed_new(coord_fine, param,row_x0,alpha_par,N_est_c) * size(modified_observation,1)
        l2_og=l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha_par,10000)
        l3=l_3_fun(modified_observation_x0, alpha_par, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha_par,0.0,1.0,50.0)
        log_likelihood_new =sum([l1,l2,l3,prior]) 
        log_likelihood_new_og =sum([l1,l2_og,l3,prior]) 
        log_likelihood_new_new =sum([l1,l2_new,l3,prior])
