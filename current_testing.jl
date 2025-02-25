include("functions.jl")



#inputs
gridsize=20 #lenght of fine grid
N_fine=gridsize^2 #number of fine grid points
N_coarse=10 #number of weather stations/ conditioning points, obersavation points
num_sim=100#number of simulated realizations

#true params for simulation
alpha_true = 1.0
beta_true=1.0
c_true=1.0
param=[ c_true , beta_true]
alpha=1.0

#Threshold definition as quantile
#p=0.0
#threshold= (1-p)^(-1/alpha)
threshold=1

#MCMC params
N_MCMC=20
param_start=[2.0,0.1]
N_est_c=1000
N_est_cond=5
N_burn_in=0

#create grids


#safe grid and observation points, also plots if last argument is true
#(coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation(gridsize,N_coarse, true)
(coord_coarse, coord_fine, row_x0)=Create_Grid_and_Observation_on_fine_grid(gridsize,N_coarse, true)

#gives the nearest fine gride coordinates for each coarse grid obsrvation 
#last argument is the coarse coordinates rounded to the nearest gridpoint coordinates
coord_cond_rows=get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)



#Simulate data on grid



#simulate data on all points and reduce it to observation data (  coarse observations)
#cholmat=chol_mat(vcat(coord_fine, coord_coarse), x->vario(x,param))
cholmat=chol_mat(coord_fine, x->vario(x,param))
#@time(sim_data= [simu_specfcts(vcat(coord_fine, coord_coarse), x->vario(x,param), cholmat, alpha_true)  for i in 1:num_sim])
@time(sim_data= [simu_specfcts(coord_fine, x->vario(x,param), cholmat, alpha_true)  for i in 1:num_sim])
sim_data=reduce(hcat,sim_data)' #just make vector of vectors a matrix (same below for observations)
#observation_data=reduce(hcat,[sim_data[i,N_fine+1:N_fine+N_coarse] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
observation_data=reduce(hcat,[sim_data[i,coord_cond_rows] for i in 1:num_sim])' #first argument is number of sim, second is coordinate
observation_x0=reduce(hcat,[sim_data[i,row_x0] for i in 1:num_sim])'
#no normalization
#normalized_observation_data=reduce(hcat,[observation_data[i,1:N_coarse]./observation_x0[i,1] for i in 1:num_sim])'




#play with conditioning
coord_coarse=hcat(collect(0:(gridsize-1)),collect(0:(gridsize-1)))./gridsize

coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)

coord_coarse=coord_fine[coord_cond_rows,:]
cond_vec=vcat(repeat([1],50),repeat([-1],47))
cond_vec_complete=vcat(repeat([1],50),repeat([-1],50))
field=log.(r_log_gaussian(coord_fine, [1.0,1.2],row_x0))
field=log.(r_cond_log_gaussian(exp.(cond_vec), 1, coord_fine, coord_fine[coord_cond_rows,:],[1.0,1.6], row_x0 ))
field2=reshape(field, gridsize,gridsize)
tx=ty=(1:gridsize)/gridsize
i=60


plotd=surface(tx,ty,field2,title="$gridsize × $gridsize FBM simulation, α=$(param[2]), c= $(param[1])")
#plotd=surface(tx,ty,field2,title="$gridsize × $gridsize FBM simulation, α=$(param[2]), c= $(param[1])")
plot!(vec(tx),vec(ty), vec((zeros(gridsize).+1)))

#import Pkg; Pkg.add("PyPlot")
using PyPlot

gr()
anim = Animation()
for i in range(0, stop = 90, step = 1)
    p = surface(tx,ty,field2,title="$gridsize × $gridsize FBM simulation, α=$(param[2]), c= $(param[1])",camera=(i, i/4))
    #plot!(vec(tx),vec(ty), vec((zeros(gridsize).+1)))
    plot!(vec(tx),vec(ty), cond_vec_complete, line = (:black, 5, 0.2))
    #Plots.plot(sol, vars=(1, 2, 3), camera=(i, i))
    frame(anim, p)
end
gif(anim, "gr1.gif", fps=24)

pyplot()
anim = Animation()




scatter(rand(10),rand(10),rand(10))

for i in coord_fine[coord_cond_rows,:]
    println(i)
end

r_cond_log_gaussian(coarse_observation)

Diagonal()


# @time(vec(FBM_simu_fast(param,gridsize,10000)[1]'))
# @time(FBM_simu_fast_vec(param,gridsize,10000)[1])


# @time(
#     for i in 1:1000
#         r_log_gaussian(coord_fine,param,row_x0)
#     end)
# r_log_gaussian_vec(coord_fine,param,row_x0,1000)

# @time( for i in 1:10000
#         r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
#         end) 
#@time (r_cond_log_gaussian_vec(observation_data,observation_x0, coord_fine,coord_coarse,param,row_x0,100) )#coord_x0 (hier c egal)

# size(observation_data,1)



@time(r_log_gaussian_vec(coord_fine,param,row_x0,1))
# #cov_mat_for_vectors(coord_fine, coord_coarse, param, coord_fine[row_x0,:])
 @time (
    

    
   a= [    [ r_cond_log_gaussian(observation_data[j,:],observation_x0[j], coord_fine,coord_coarse,param,row_x0) for rep in 1:1000 ]  for j in 1:size(observation_x0,1) ];#coord_x0 (hier c egal)
  


 )

 @time (b=r_cond_log_gaussian_vec(observation_data,observation_x0, coord_fine,coord_coarse,param,row_x0, 1000); ) #coord_x0 (hier c egal)

b[1][1]


mean(mean(mean(a.-b))) #@time (modified_observation, modified_observation_x0) =exceed_cond_sim(10,num_sim,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )

@time (a=exceed_cond_sim(100,num_sim,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 ) )
  
(modified_observation, modified_observation_x0) = exceed_cond_sim(100,num_sim,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )

threshold=1
(-alpha-1)*sum(log.(modified_observation_x0./threshold))

@time (b=exceed_cond_sim_vec(100,num_sim,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 ) )

N_est_c=100000
number_of_excced=size(b[1],1)
@time( a=l_2_fun(coord_fine, param,row_x0, number_of_excced,alpha,N_est_c) )
@time( b=l_2_fun_vec(coord_fine, param,row_x0, number_of_excced,alpha,N_est_c) )
 #@time exceed_cond_sim_with_all_info(10,num_sim,normalized_log_observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )

# (modified_observation, modified_observation_x0) = exceed_cond_sim(10,num_sim,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )
N_est_cond=1000
@time ( a = l_3_fun(coord_fine, coord_coarse, param, row_x0, observation_data, observation_x0, alpha, N_est_cond) )
@time ( b = l_3_fun_vec(coord_fine, coord_coarse, param, row_x0, observation_data, observation_x0, alpha, N_est_cond) )

#exceed_cond_sim_fast(10,num_sim,normalized_log_observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )
                        

#N=size(coord_fine,1)
#uncond_rows=setdiff(1:gridsize*gridsize,coord_cond_rows)
#coord_cond=coord_fine[coord_cond_rows,:]
#coord_uncond=coord_fine[uncond_rows,:]
#cov_mat_coarse_inv_alt=inv(cov_mat_for_vectors(coord_cond,coord_cond, param,  coord_fine[row_x0,:])) 

#  @time l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, observation_x0, row_x0)

#  N_est_c=1000000
# @time l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)


#  N_est_cond=50

# @time l3=l_3_fun(coord_fine, coord_coarse, param, row_x0, modified_observation, observation_x0, alpha, N_est_cond)   

#  log_likelihood_old=sum([l1,l2,l3])















#some FBM testing
#Random.seed!(1234)
@time(FBM_simu_fast_vec(param, gridsize,1))
#Random.seed!(1234)
FBM_res=FBM_simu_fast_vec(param,gridsize,num_sim)
FBM_res[1][1]
r_log_gaussian(coord_fine,param,row_x0)




#empirical and true cov matrices to test FBM 
emp_cov_mat=zeros(N_fine,N_fine)
true_cov_mat=zeros(N_fine,N_fine)
mean_vec=zeros(N_fine)
for i in 1:N_fine
    mean_vec[i]=mean([FBM_res[rep][i] for rep in 1:num_sim])
end


for row in 1:N_fine
    for col in 1:N_fine
        emp_cov_mat[row,col]=1/(num_sim-1)*sum([(FBM_res[rep][row]-mean_vec[row])*(FBM_res[rep][col]-mean_vec[col]) for rep in 1:num_sim])
    end
end


for row in 1:N_fine
    for col in 1:N_fine
        true_cov_mat[row,col]=param[1]*(norm(coord_fine[row,:])^param[2]+norm(coord_fine[col,:])^param[2]-norm(coord_fine[row,:]-coord_fine[col,:])^param[2])
    end
end


maximum(
    abs.((emp_cov_mat-true_cov_mat)./(true_cov_mat.+10^(-8)))
    )