    include("functions.jl")

    #inputs
    gridsize=5#lenght of fine grid
    N_fine=gridsize^2 #number of fine grid points
    N_coarse=5 #number of weather stations/ conditioning points, obersavation points
    num_sim=1000 #number of simulated realizations
    num_runs=10000

    #true params for simulation
    alpha_true = 2.0
    beta_true=1.5
    c_true=2.0
    param=[ c_true , beta_true]
    alpha=2.0

    #threshold_method can be "quantile", "fixed" or "all_exceed
    threshold_method="quantile"

    #Threshold definition as quantile
    threshold_quantile=0.80
    #threshold= (1-threshold_quantile)^(-1/(alpha_true+1))
    #threshold=2.2
    #threshold=1.0
    

    #MCMC params
    N_MCMC=1000
    param_start=[2.0,1.5]
    alpha_start=2.0
    N_est_c=5000
    N_burn_in=0
    N_cond_sim=100

    #create grids


    #safe grid and observation points, also plots if last argument is true
    #(coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation(gridsize,N_coarse, true)
    (coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation_on_fine_grid(gridsize,N_coarse, true)
    #gives the nearest fine gride coordinates for each coarse grid obsrvation 
    #last argument is the coarse coordinates rounded to the nearest gridpoint coordinates
    coord_cond_rows=get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
    coord_coarse=coord_fine[coord_cond_rows,:]



    #Simulate data on grid



    #simulate data on all points and reduce it to observation data (  coarse observations)
    #cholmat=chol_mat(vcat(coord_fine, coord_coarse), x->vario(x,param))
    cholmat=chol_mat(coord_fine, x->vario(x,param))
    #@time(sim_data= [simu_specfcts(vcat(coord_fine, coord_coarse), x->vario(x,param), cholmat, alpha_true)  for i in 1:num_sim])
    @time(sim_data= [simu_specfcts_MCMC(num_runs, alpha, coord_fine,param,row_x0 ) for i in 1:num_sim] )
    #@time(sim_data=[simu_specfcts_verynew(coord_fine, x->vario(x,param), cholmat, alpha_true)  for i in 1:num_sim])
    sim_data=reduce(hcat,sim_data)' #just make vector of vectors a matrix (same below for observations)
    #observation_data=reduce(hcat,[sim_data[i,N_fine+1:N_fine+N_coarse] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
    observation_data=reduce(hcat,[sim_data[i,coord_cond_rows] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
    observation_x0=reduce(hcat,[sim_data[i,row_x0] for i in 1:num_sim])'

    #run MCMC algorithm

   n_trial_print=500 #print every n-th step
   #@time (result1= MCMC(N_MCMC,observation_data,observation_x0,threshold_quantile,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim))
   
   param_start=[2.0,1.5]
   alpha_start=2.0
   threshold_method="fixed"
   threshold=quantile_threshold(10000,alpha_start,coord_fine,param_start,row_x0,0.8,10000)
   @time (result2= MCMC(N_MCMC,observation_data,observation_x0,threshold,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim))
   #threshold_method="all_exceed"
   #@time (result3= MCMC(N_MCMC,observation_data,observation_x0,0.0000001,threshold_method, alpha_start, coord_fine,coord_coarse,param_start,row_x0,n_trial_print,N_est_c,N_cond_sim))

   #result=result1
   #result=result2
   #result=result3
  #=   #mean of params with burn in period canceled out
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
 =#
 




    #code to safe everything
 
    import Pkg;Pkg.add("JLD2")
    using JLD2
    
    sim_name=Int(floor(rand()*10000000))
    save_object("sim_$(sim_name).jld2", (result1,result2,coord_coarse, coord_fine, row_x0,sim_data))

    
   # repetition=0
   # save_object("$(N_MCMC)_simulation_$(gridsize)_grid_0$repetition.jld2", (result1,result2))
   # save_object("$(num_sim)_simulated_DATA_$(gridsize)_grid_0$repetition.jld2", sim_data)


