include("functions.jl")
include("input_for_simulation.jl")



function simulate_observations(gridsize,N_coarse,num_sim,num_runs,beta_true,c_true, alpha_true,path)
#create grid and observation points, also plots if last argument is true
#here the coarse observation points are already on the fine grid
#also safe coords 
(coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation_on_fine_grid(gridsize,N_coarse, false)

#gives the nearest fine gride coordinates for each coarse grid obsrvation 
#last argument is the coarse coordinates rounded to the nearest gridpoint coordinates
coord_cond_rows=get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
#redefine it to make sure its in the right order
coord_coarse=coord_fine[coord_cond_rows,:]
#coord_x0=coord_fine[row_x0,:]

#Simulate data on grid
#simulate data on all points and reduce it to observation data (  coarse observations)
#use mcmc approach via simu_specfcts_MCMC, let chain run for 2000 steps for each observation

@time(sim_data= [simu_specfcts_MCMC(num_runs, alpha_true, coord_fine,[ c_true , beta_true],row_x0 ) for i in 1:num_sim] )
sim_data=reduce(hcat,sim_data)' #just make vector of vectors a matrix (same below for observations)
observation_data=reduce(hcat,[sim_data[i,coord_cond_rows] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
observation_x0=reduce(hcat,[sim_data[i,row_x0] for i in 1:num_sim])'

#If threshold shall be estimated, use this 
#=
threshold_quantile=0.80
@time( estimated_threshold=quantile_threshold(10000,alpha,coord_fine,param,row_x0,threshold_quantile,num_sim))
=#

param_start=[gaussian_proposal(1.0,2.0),rand()*2.0]
alpha_start=gaussian_proposal(1.0,2.0)
#save simulated observation data
number_str=highest_sim_number(readdir(path))
save_object(path*"/simulated_observations_$(number_str).jld2", (coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, param_start, alpha_start))#,estimated_threshold)) 

return(number_str,(coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, param_start, alpha_start))
end




#function to simulate observations on a 7x7 grid, with fixed 9 coarse observation points
function simulate_observations_fixed_seven_by_seven(gridsize,N_coarse,num_sim,num_runs,beta_true,c_true, alpha_true,path)
#create grid and observation points, also plots if last argument is true
#here the coarse observation points are already on the fine grid
#also safe coords 
(coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation_seven_times_seven(false)

#gives the nearest fine gride coordinates for each coarse grid obsrvation 
#last argument is the coarse coordinates rounded to the nearest gridpoint coordinates
coord_cond_rows=get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
#redefine it to make sure its in the right order
coord_coarse=coord_fine[coord_cond_rows,:]
#coord_x0=coord_fine[row_x0,:]

#Simulate data on grid
#simulate data on all points and reduce it to observation data (  coarse observations)
#use mcmc approach via simu_specfcts_MCMC, let chain run for 2000 steps for each observation

@time(sim_data= [simu_specfcts_MCMC(num_runs, alpha_true, coord_fine,[ c_true , beta_true],row_x0 ) for i in 1:num_sim] )
sim_data=reduce(hcat,sim_data)' #just make vector of vectors a matrix (same below for observations)
observation_data=reduce(hcat,[sim_data[i,coord_cond_rows] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
observation_x0=reduce(hcat,[sim_data[i,row_x0] for i in 1:num_sim])'

#If threshold shall be estimated, use this 
#=
threshold_quantile=0.80
@time( estimated_threshold=quantile_threshold(10000,alpha,coord_fine,param,row_x0,threshold_quantile,num_sim))
=#

param_start=[gaussian_proposal(1.0,2.0),rand()*2.0]
alpha_start=gaussian_proposal(1.0,2.0)
#save simulated observation data
number_str=highest_sim_number(readdir(path))
save_object(path*"/simulated_observations_$(number_str).jld2", (coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, param_start, alpha_start))#,estimated_threshold)) 

return(number_str,(coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, param_start, alpha_start))
end