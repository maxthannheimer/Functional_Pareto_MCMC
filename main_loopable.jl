include("functions.jl")
include("input_for_simulation.jl")
include("simulate_observation_data.jl")
  


    #For running multiple scripts in parallel we make different folders called sim_res_folder_xx for each parallel run
    # Make sure the results_dir exists
    if !isdir("results_dir")
    mkdir("results_dir")
    end
    number_parallel=highest_folder_number( readdir("results_dir")  )
    path=mkdir("results_dir/sim_res_folder_$(number_parallel)")
#run simu forever
while(true)
    #Simulate observation data
    number_str,(coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, param_start, alpha_start)=simulate_observations_fixed_seven_by_seven(gridsize,N_coarse,num_sim,num_runs,beta_true,c_true, alpha_true,path)
    #our method with conditional simulation
    #run MCMC algorithm
    input_MCMC_cond_sim=(N_MCMC,observation_data,observation_x0,threshold,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim)
    @time (result_cond_sim= MCMC(N_MCMC,observation_data,observation_x0,threshold,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim))
    #reload "simulated_observations_$number_str.jld2" to be sure the staring values are the same
    (coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, param_start, alpha_start )=load(path*"/simulated_observations_$number_str.jld2")["single_stored_object"]
    #standard method with approx risk functional
    input_MCMC_approx=(N_MCMC,observation_data,observation_x0,threshold, alpha_start, coord_fine[row_x0,:],coord_coarse,param_start,n_trial_print,N_est_c_approx)
    @time (result_approx= MCMC_approx(N_MCMC,observation_data,observation_x0,threshold, alpha_start, coord_fine[row_x0,:],coord_coarse,param_start,n_trial_print,N_est_c_approx) )
    #code to safe everything
    save_object(path*"/sim_$(number_str).jld2", (result_cond_sim,result_approx,input_MCMC_cond_sim,input_MCMC_approx))
end
