include("functions.jl")

N_burn_in=1000
N_MCMC=8000
true_param=Dict("beta" => 1.5, "c" => 2.0, "alpha" => 2.0)
number_parallel_int_max=parse(Int,highest_folder_number( readdir("results_dir")  ))-1
number_parallel=number_parallel_int_max #number of parallel runs, here we assume that the folders are named sim_res_folder_01, sim_res_folder_02, etc.
number_sim_int_max=parse(Int,highest_sim_result_number(readdir("results_dir/sim_res_folder_01")))-1

est_cond_sim_mean=Dict( "beta" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "c" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max])
est_cond_sim_median=Dict( "beta" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "c" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max])
est_approx_mean=Dict( "beta" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "c" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max])
est_approx_median=Dict( "beta" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max],
                    "c" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max], 
                    "alpha" => [0.0 for i in 1:number_parallel, j in 1:number_sim_int_max])


for number_parallel_int in 1:number_parallel_int_max
    println(number_parallel_int)
    sim_res_folder_name=readdir("results_dir")[number_parallel_int]
    for number_sim_int in 1:number_sim_int_max
        path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
        (result_cond_sim,result_approx)=load(path)["single_stored_object"]
        for key in keys(est_approx_mean)
            est_cond_sim_mean[key][number_parallel_int,number_sim_int]=mean(result_cond_sim[key][N_burn_in:end])
            est_cond_sim_median[key][number_parallel_int,number_sim_int]=median(result_cond_sim[key][N_burn_in:end])
            est_approx_mean[key][number_parallel_int,number_sim_int]=mean(result_approx[key][N_burn_in:end])
            est_approx_median[key][number_parallel_int,number_sim_int]=median(result_approx[key][N_burn_in:end])
        end
    end 
end

key="beta"
key="c"
key="alpha"
            est_cond_sim_mean[key]
            est_cond_sim_median[key]
            est_approx_mean[key]
            est_approx_median[key]


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
            RMSE_cond_sim_mean[key]=mean((est_cond_sim_mean[key].-true_param[key]).^2)^0.5
            RMSE_cond_sim_median[key]=mean((est_cond_sim_median[key].-true_param[key]).^2)^0.5
            RSME_est_approx_mean[key]=mean((est_approx_mean[key].-true_param[key]).^2)^0.5
            RMSE_est_approx_median[key]=mean((est_approx_median[key].-true_param[key]).^2)^0.5
end

RMSE_est_approx_median
RSME_est_approx_mean

RMSE_cond_sim_median
RMSE_cond_sim_mean




plots = Vector{}(undef, 6)
for number_sim_int in 1:3
number_parallel_int=3
sim_res_folder_name=readdir("results_dir")[number_parallel_int]

path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
(result_cond_sim,result_approx)=load(path)["single_stored_object"]  
for res in [(result_cond_sim,"cond_sim",1),(result_approx,"approx_risk",4)]
    (result,method,i)=res


#histogram(result["beta"][N_burn_in+1:N_MCMC+1],title="Histogram for β, Burn in: $N_burn_in")
p_beta=scatter(1:N_MCMC+1,result["beta"],label="$N_MCMC samples of beta",title="Markov chain for β($method), Burn in: $N_burn_in")
hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β: $(round((mean(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_beta, [median(result["beta"][N_burn_in+1:N_MCMC+1])],label="Median estimate for β: $(round((median(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_beta, [true_param["beta"] for i in N_burn_in+1:N_MCMC+1], label="True value for β: $(true_param["beta"]) ",linewidth=3)

#histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
p_c=scatter(1:N_MCMC+1,result["c"],label="$N_MCMC samples of c",title="Markov chain for c ($method), Burn in: $N_burn_in")
hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c: $(round((mean(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_c, [median(result["c"][N_burn_in+1:N_MCMC+1])],label="Median estimate for c: $(round((median(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_c, [(true_param["c"]) for i in N_burn_in+1:N_MCMC+1], label="True value for c: $(true_param["c"]) ",linewidth=3)

#histogram(result["alpha"][N_burn_in+1:N_MCMC+1],title="Histogram for α, Burn in: $N_burn_in")
p_alpha=scatter(1:N_MCMC+1,result["alpha"],label="$N_MCMC samples of α",title="Markov chain for α ($method), Burn in: $N_burn_in")
hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α: $(round((mean(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_alpha, [median(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Median estimate for α: $(round((median(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_alpha, [(true_param["alpha"]) for i in N_burn_in+1:N_MCMC+1], label="True value for α: $(true_param["alpha"]) ",linewidth=3)

plots[i+0],plots[i+1],plots[i+2]=p_beta,p_alpha,p_c

end
combined = plot(plots[1],plots[4],plots[2],plots[5],plots[3],plots[6], layout = (3,2),size=(1500,1000))
savefig(combined, "6_plots_$number_sim_int.pdf")
end


number_sim_int=6
number_parallel_int=17
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
    number_parallel_int=17
    sim_res_folder_name=readdir("results_dir")[number_parallel_int]
    number_sim_int=6
        path="results_dir/"*sim_res_folder_name*"/sim_"* @sprintf("%04d", number_sim_int) *".jld2"
        (result_cond_sim,result_approx)=load(path)["single_stored_object"]
scatter(1:8000,result_cond_sim["Number of exceedance"])


plot(1:1001,result_cond_sim["alpha"][6000:1:7000])