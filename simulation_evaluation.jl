include("functions.jl")
using Measures
N_burn_in=2000
N_MCMC=10000
true_param=Dict("beta" => 0.5, "c" => 3.0, "alpha" => 2.0)
number_parallel_int_max=parse(Int,highest_folder_number( readdir("results_dir")  ))-1 #number of parallel runs, here we assume that the folders are named sim_res_folder_01, sim_res_folder_02, etc.
number_sim_int_max=parse(Int,highest_sim_result_number(readdir("results_dir/sim_res_folder_01")))-1

#two simulations stopped since there were no exceedances
#acceptable_simulations=[i for i in 1:number_parallel_int_max if i!=9 && i!=22] #the simulations 10 and 36 were stopped since there where no exceedances
acceptable_simulations=[i for i in 1:number_parallel_int_max]

est_cond_sim_mean=Dict( "beta" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "c" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max])
est_cond_sim_median=Dict( "beta" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "c" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max])
est_approx_mean=Dict( "beta" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "c" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max])
est_approx_median=Dict( "beta" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max],
                    "c" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in acceptable_simulations, j in 1:number_sim_int_max])


for number_parallel_int in 1:size(acceptable_simulations,1)

    println(acceptable_simulations[number_parallel_int])
    sim_res_folder_name=readdir("results_dir")[acceptable_simulations[number_parallel_int]]
    for number_sim_int in 1:number_sim_int_max
        path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
        (result_cond_sim,result_approx)=load(path)["single_stored_object"]
        if number_parallel_int==22 ||number_parallel_int==18 || number_parallel_int==9 #ignore simulations with warnings
            for key in keys(est_approx_mean)
                est_cond_sim_mean[key][number_parallel_int,number_sim_int]=true_param[key]
                est_cond_sim_median[key][number_parallel_int,number_sim_int]=true_param[key]
                est_approx_mean[key][number_parallel_int,number_sim_int]=true_param[key]
                est_approx_median[key][number_parallel_int,number_sim_int]=true_param[key]
            end
        else
            for key in keys(est_approx_mean)
                est_cond_sim_mean[key][number_parallel_int,number_sim_int]=mean(result_cond_sim[key][N_burn_in:end])
                est_cond_sim_median[key][number_parallel_int,number_sim_int]=median(result_cond_sim[key][N_burn_in:end])
                est_approx_mean[key][number_parallel_int,number_sim_int]=mean(result_approx[key][N_burn_in:end])
                est_approx_median[key][number_parallel_int,number_sim_int]=median(result_approx[key][N_burn_in:end])
            end
        end
    end 
end



key="beta"
key="c"
key="alpha"
            est_cond_sim_mean[key].-true_param[key]
            est_cond_sim_median[key].-true_param[key]
            est_approx_mean[key].-true_param[key]
            est_approx_median[key].-true_param[key]


RMSE_cond_sim_mean=Dict( "beta" => 0.0, 
                    "c" => 0.0, 
                    "alpha" => 0.0)
RMSE_cond_sim_median=Dict( "beta" => 0.0,
                    "c" => 0.0, 
                    "alpha" => 0.0)
RSME_est_approx_mean=Dict( "beta" => 0.0, 
                    "c" => 0.0, 
                    "alpha" => 0.0)
RMSE_est_approx_median=Dict( "beta" => 0.0,
                    "c" => 0.0, 
                    "alpha" => 0.0)

for key in keys(est_approx_mean)
            RMSE_cond_sim_mean[key]=round(mean((est_cond_sim_mean[key].-true_param[key]).^2)^0.5 , digits=3)
            RMSE_cond_sim_median[key]=round(mean((est_cond_sim_median[key].-true_param[key]).^2)^0.5 , digits=3)
            RSME_est_approx_mean[key]=round(mean((est_approx_mean[key].-true_param[key]).^2)^0.5 , digits=3)
            RMSE_est_approx_median[key]=round(mean((est_approx_median[key].-true_param[key]).^2)^0.5 , digits=3)
end

RMSE_est_approx_median
RSME_est_approx_mean

RMSE_cond_sim_median
RMSE_cond_sim_mean


scalefontsizes(1*1.1) 
for number_parallel_int in 17:1:19
plots = Vector{}(undef, 6)
for number_sim_int in 5:1:5

sim_res_folder_name=readdir("results_dir")[number_parallel_int]

path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
(result_cond_sim,result_approx)=load(path)["single_stored_object"]  
max_beta=max(maximum(result_cond_sim["beta"]),maximum(result_approx["beta"]))+0.2
max_c=max(maximum(result_cond_sim["c"]),maximum(result_approx["c"]))+0.2
max_alpha=max(maximum(result_cond_sim["alpha"]),maximum(result_approx["alpha"]))+0.2
for res in [(result_cond_sim,"two-step-exact",1),(result_approx,"one-step-approximate",4)]
    (result,method,i)=res


#histogram(result["beta"][N_burn_in+1:N_MCMC+1],title="Histogram for β, Burn in: $N_burn_in")
p_beta = scatter(
    1:2:N_MCMC+1, result["beta"];
    label="",#"$N_MCMC samples of beta",
    title="Markov chain for β",#, Burn in: $N_burn_in",
    ylims=(0, max_beta),
    marker=:circle,
    markercolor=:white,           # non-filled
    markerstrokecolor=:black,     # circle outline color
    markerstrokewidth=0.2,        # thin outline
    markersize=4,
    legend=:topright, 
    margin=5mm
)
vline!(p_beta, [N_burn_in], label="",linewidth=3, linestyle=:dot)
#hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β",linewidth=3)
#hline!(p_beta, [median(result["beta"][N_burn_in+1:N_MCMC+1])],label="Median estimate for β: $(round((median(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_beta, [true_param["beta"] for i in N_burn_in+1:N_MCMC+1], label="True value for β",linewidth=3, linestyle=:dot)
#vline!(p_beta, [N_burn_in], label="Burn in: $N_burn_in",linewidth=3, linestyle=:dot)
hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β: $(round((mean(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_beta, [median(result["beta"][N_burn_in+1:N_MCMC+1])],label="Median estimate for β: $(round((median(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
#hline!(p_beta, [true_param["beta"] for i in N_burn_in+1:N_MCMC+1], label="True value for β: $(true_param["beta"]) ",linewidth=3, linestyle=:dash)

#histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
p_c = scatter(
    1:2:N_MCMC+1, result["c"];
    label="", #"$N_MCMC samples of c",
    title="Markov chain for c",#, Burn in: $N_burn_in",
    ylims=(1.0, max_c),
    marker=:circle,
    markercolor=:white,           # non-filled
    markerstrokecolor=:black,     # circle outline color
    markerstrokewidth=0.2,        # thin outline
    markersize=4,
    legend=:bottomright, 
    margin=5mm
)

vline!(p_c, [N_burn_in], label="",linewidth=3, linestyle=:dot)
#hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c",linewidth=3)
hline!(p_c, [(true_param["c"]) for i in N_burn_in+1:N_MCMC+1], label="True value for c",linewidth=3, linestyle=:dot)
#vline!(p_c, [N_burn_in], label="Burn in: $N_burn_in",linewidth=3, linestyle=:dot)
hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c: $(round((mean(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_c, [median(result["c"][N_burn_in+1:N_MCMC+1])],label="Median estimate for c: $(round((median(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
#hline!(p_c, [(true_param["c"]) for i in N_burn_in+1:N_MCMC+1], label="True value for c: $(true_param["c"]) ",linewidth=3, linestyle=:dash)



#histogram(result["alpha"][N_burn_in+1:N_MCMC+1],title="Histogram for α, Burn in: $N_burn_in")
p_alpha = scatter(
    1:2:N_MCMC+1, result["alpha"];
    label="",#"$N_MCMC samples of α",
    title="Markov chain for α", #, Burn in: $N_burn_in",
    ylims=(0.8, max_alpha),
    marker=:circle,
    markercolor=:white,           # non-filled
    markerstrokecolor=:black,     # circle outline color
    markerstrokewidth=0.2,        # thin outline
    markersize=4,
    legend=:bottomright,
    margin=5mm
)

vline!(p_alpha, [N_burn_in], label="",linewidth=3, linestyle=:dot)
#hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α",linewidth=3)
hline!(p_alpha, [(true_param["alpha"]) for i in N_burn_in+1:N_MCMC+1], label="True value for α",linewidth=3, linestyle=:dot)
#vline!(p_alpha, [N_burn_in], label="Burn in: $N_burn_in",linewidth=3, linestyle=:dot)
hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α: $(round((mean(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_alpha, [median(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Median estimate for α: $(round((median(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
#hline!(p_alpha, [(true_param["alpha"]) for i in N_burn_in+1:N_MCMC+1], label="True value for α: $(true_param["alpha"]) ",linewidth=3, linestyle=:dash)


plots[i+0],plots[i+1],plots[i+2]=p_beta,p_c,p_alpha

end
combined = plot(plots[1],plots[4],plots[2],plots[5],plots[3],plots[6],plot_title="      Two-step-exact                One-step-approximate", layout = (3,2),size=(1500,1000))
#combined = plot(plots[1],plots[2],plots[3],layout = (3,1),size=(1500,500))
savefig(combined, "3_plots_$(number_parallel_int)_$(number_sim_int).pdf")
end
end

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
N_est_c=10000
result_cond_sim["log_likelihood"][n]
param=[result_cond_sim["c"][n],result_cond_sim["beta"][n]]
alpha_par=result_cond_sim["alpha"][n]
  num_obs=size(observation_data,1)
(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha_par, coord_fine,coord_coarse,param,row_x0 )

(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, 2.0, coord_fine,coord_coarse,[3.0,0.5],row_x0 )
(modified_observation, modified_observation_x0)=(modified_observation[2:end,:], modified_observation_x0[2:end])
param=[100,0.1]
alpha_par=50
#calculate new log likelihood
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha_par)
        #l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        @time( l2=l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha_par,N_est_c) * size(modified_observation,1)) #number of exceedances
        l2
        @time(l2_og=l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha_par,N_est_c))
        l2_og
        l3=l_3_fun(modified_observation_x0, alpha_par, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha_par,0.0,1.0,50.0)
        log_likelihood_new =sum([l1,l2,l3,prior]) 
        log_likelihood_new_og =sum([l1,l2_og,l3,prior]) 

param=[3.0,0.5]
alpha_par=2.0
  num_obs=size(observation_data,1)
(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha_par, coord_fine,coord_coarse,param,row_x0 )
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha_par) 
        l2=l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha_par,N_est_c) * size(modified_observation,1) 
        l2_new=l_2_fun_dependent_no_number_of_exceed_new(coord_fine, param,row_x0,alpha_par,N_est_c) * size(modified_observation,1)
        l2_og=l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha_par,N_est_c)
        l3=l_3_fun(modified_observation_x0, alpha_par, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha_par,0.0,1.0,50.0)
        log_likelihood_new =sum([l1,l2,l3,prior]) 
        log_likelihood_new_og =sum([l1,l2_og,l3,prior]) 
        
        #l_2_fun_dependent
        N_est_c=40000
    tmp = r_log_gaussian_vec_dependent(coord_fine,param,row_x0, N_est_c,alpha_par) 
    
    -size(modified_observation,1) * log(mean([mean(tmp[i] 
        )^(alpha_par) for i in 1:N_est_c]))  #* size(modified_observation,1) #number of exceedances
     # minus for 1/c_l (in log)
    surface(coord_fine[:,1],coord_fine[:,2],tmp[1],title="FBM simulation",xlabel="x",ylabel="y",zlabel="z")
sum(tmp[1])
    tmp[1][row_x0]
    mean(tmp[1] )
    [mean(tmp[i] ) for i in 1:N_est_c]
    [mean(tmp[i] )^(alpha_par) for i in 1:N_est_c]
    mean([mean(tmp[i] )^(alpha_par) for i in 1:N_est_c])
    log(mean([mean(tmp[i] )^(alpha_par) for i in 1:N_est_c]))
     -size(modified_observation,1) *  log(mean([mean(tmp[i] )^(alpha_par) for i in 1:N_est_c]))

        gridsize = Int(sqrt(size(coord_fine,1)))
    res = FBM_simu_fast_vec(param, gridsize,N_est_c)
    trend=vec_vario(param,coord_fine,coord_fine[row_x0,:])
    for i in 1:N_est_c
            res[i] = exp.(1/alpha_par*(res[i] - trend .-res[i][row_x0])) #variogram
    end
    res
    mean(res[22])


    l_1_fun(coord_fine,coord_coarse,coarse_observation,param, observation_x0, row_x0,alpha_par)
    l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha_par)
       
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.* gridsize)./gridsize)
    cov_mat_coarse_inv= inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    inv_determinant = det(cov_mat_coarse_inv)
    trend = -vec_vario(param,coord_fine[coord_cond_rows,:],coord_fine[row_x0,:])
    sum([(log_d_gaussian(trend ,cov_mat_coarse_inv , alpha_par*log.(modified_observation[i,:]./modified_observation_x0[i]), inv_determinant))+log(alpha_par)*length(trend) for i in 1:size(modified_observation,1)])

  