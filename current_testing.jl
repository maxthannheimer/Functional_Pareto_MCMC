include("functions.jl")
using JLD2


(coord_coarse, coord_fine, row_x0,sim_data, observation_data, observation_x0, alpha, param)=load("simulated_observations.jld2")["single_stored_object"]

#inputs
gridsize=size(coord_fine,1)#lenght of fine grid
N_fine=gridsize^2 #number of fine grid points
N_coarse=size(coord_coarse,1) #number of weather stations/ conditioning points, obersavation points
num_sim=size(observation_data,1) #number of simulated realizations

#Numbers for simulation based estimation

N_est_c=1000
N_est_cond=5
N_burn_in=1000



coord_x0=coord_fine[row_x0,:]
coord_coarse_plus=vcat(coord_coarse,coord_x0')
num_rep=200
r_gaussian_vec_coarse(coord_coarse_plus,coord_x0,param,row_x0,num_rep,alpha)
cov_mat=cov_mat_for_vectors(coord_coarse_plus, coord_coarse_plus,param,coord_x0)
cov_mat-cov_mat'

function r_gaussian_vec_coarse(coord_coarse_plus,coord_x0,param,row_x0,num_rep,alpha) 
    N = size(coord_coarse_plus,1)
    cov_mat = cov_mat_for_vectors(coord_coarse_plus, coord_coarse_plus,param,coord_x0).+1e-6 
    res = rand(MvNormal([0.0 for i in 1:N],cov_mat),num_rep)
    trend=vec_vario(param,coord_coarse_plus,coord_coarse_plus[row_x0,:])
    for i in 1:num_rep
            res[:,i] = 1/alpha*(res[:,i] - trend .-res[:,end]) #variogram
    end
    res
end







#DIFFERENCE BETWEEN DEPENDEND AND INDEPENDENT SIMULATION IN FBM SIMULATION
# Circulant embedding testing will follow in seperate file
n_test=100
num_sim=100
@time begin
tmp=[mean(FBM_simu_fast_vec_dependent(param,gridsize,num_sim)) for i in 1:n_test]
println("dependend simulation")
end


@time begin
tmp=[mean(FBM_simu_fast_vec(param,gridsize,num_sim)) for i in 1:n_test]
println("independet simulation")
end


# Estimation of 1/c via samples
number_of_exceed = 200
plots1 = Vector{}(undef, 4)
i=1
for N_est_c in [10,100,1000,5000]
println("time for $N_est_c:")
@time(tmp=[c_estimation(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) for i in 1:100])
#pl=scatter(1:100,sort(tmp),label="l_2_fun",title="l_2_fun for $N_est_c samples")
#hline!(pl, [mean(tmp)],label="Mean of estimates for l_2_fun: $(mean(tmp))",linewidth=3)
plots1[i]=histogram(tmp,title="$N_est_c samples, mean: $(mean(tmp)) ")
i=i+1
end
combined=plot(plots1..., layout = (2,2),size=(1000,1000))
savefig(combined, "c_estimation_independent.pdf")

number_of_exceed = 200
plots2 = Vector{}(undef, 4)
i=1
for N_est_c in [10,100,1000,5000]
println("time for $N_est_c:")
@time(tmp=[c_estimation_dependent(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) for i in 1:100])
#pl=scatter(1:100,sort(tmp),label="l_2_fun",title="l_2_fun for $N_est_c samples")
#hline!(pl, [mean(tmp)],label="Mean of estimates for l_2_fun: $(mean(tmp))",linewidth=3)
plots2[i]=histogram(tmp,title="$N_est_c samples, mean: $(mean(tmp)) ")
i=i+1
end
combined2=plot(plots2..., layout = (2,2),size=(1000,1000))
savefig(combined2, "c_estimation_dependent.pdf")

# Estimation of 1/c via samples
number_of_exceed = 200
plots1 = Vector{}(undef, 4)
i=1
for N_est_c in [10,100,1000,5000]
println("time for $N_est_c:")
@time(tmp=[l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) for i in 1:100])
#pl=scatter(1:100,sort(tmp),label="l_2_fun",title="l_2_fun for $N_est_c samples")
#hline!(pl, [mean(tmp)],label="Mean of estimates for l_2_fun: $(mean(tmp))",linewidth=3)
plots1[i]=histogram(tmp,title="$N_est_c samples, mean: $(mean(tmp)) ")
i=i+1
end
combined=plot(plots1..., layout = (2,2),size=(1000,1000))
savefig(combined, "l2_fun_independent.pdf")

number_of_exceed = 200
plots2 = Vector{}(undef, 4)
i=1
for N_est_c in [10,100,1000,10000]
println("time for $N_est_c:")
@time(tmp=[l_2_fun_dependent(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) for i in 1:100])
#pl=scatter(1:100,sort(tmp),label="l_2_fun",title="l_2_fun for $N_est_c samples")
#hline!(pl, [mean(tmp)],label="Mean of estimates for l_2_fun: $(mean(tmp))",linewidth=3)
plots2[i]=histogram(tmp,title="$N_est_c samples, mean: $(mean(tmp)) ")
i=i+1
end
plot(plots2..., layout = (2,2),size=(1000,1000))
savefig(combined, "l_2_fun_dependent.pdf")