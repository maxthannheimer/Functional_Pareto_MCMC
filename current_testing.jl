include("functions.jl")


#inputs
gridsize=5#lenght of fine grid
N_fine=gridsize^2 #number of fine grid points
N_coarse=5 #number of weather stations/ conditioning points, obersavation points
num_sim=1000 #number of simulated realizations

#true params for simulation
alpha_true = 2.0
beta_true=1.5
c_true=2.0
param=[ c_true , beta_true]
alpha=2.0

#Threshold definition as quantile
#p=0.80
#threshold= (1-p)^(-1/(alpha_true+1))
#threshold=0.01
threshold_quantile=0.7

#MCMC params
N_MCMC=1000
param_start=[2.0,1.5]
alpha_start=2.0
N_est_c=1000
N_est_cond=10
N_burn_in=0

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

num_runs=1000
num_sim=1000
alpha=1



r_log_gaussian_vec(coord_fine,param,row_x0,5,alpha) 




function simu_specfcts_MCMC(num_runs, alpha, coord_fine,param,row_x0 )
    tmp = r_log_gaussian_vec(coord_fine,param,row_x0,num_runs+1,alpha) 
    old_value = tmp[1]
    #println(old_value)
   #old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
    for trial in 1:num_runs
       #direkt num_obs viele simulations
            # if (trial==1)
            #         old_value = tmp[i][trial+1]
            # end
            proposal = tmp[trial+1]
            acceptance_rate = min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
            #res_ell_X[i]=observation_x0[i]*mean(old_value)
    end
    proposal=proposal/mean(proposal)
    proposal*=(1/(1-rand()))^(1/alpha)
    #proposal=proposal*(1/(1-rand()))^(1/alpha)
end

@time ( res=[simu_specfcts_MCMC_single_sim(num_runs, alpha, coord_fine,param,row_x0 ) for i in 1:num_sim] )


@time(sim_data= [simu_specfcts_verynew(coord_fine, x->vario(x,param), cholmat, alpha_true)  for i in 1:num_sim])
sim_data


r_log_gaussian_vec(coord_fine,param,row_x0,1,alpha)[1]
r_cond_log_gaussian_vec(coarse_observation,observation_x0, coord_fine,coord_coarse,param,row_x0,num_rep,alpha) #coord_x0 (hier c egal)
   
@time (simu_specfcts_MCMC_multisim(1000,1000,  alpha, coord_fine,param,row_x0 ))

a=[2, 2]
a*=7