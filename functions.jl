using LinearAlgebra
#import Pkg; Pkg.add("Distributions")
using Random, Distributions
using Plots
using FFTW


# creating the grid of fine points and the coarse observation points
# boolean plot makes plot of grid if true
function Create_Grid_and_Observation(gridsize,N_coarse, plot::Bool=false )
    #N_all=gridsize^2+N_coarse+1
    coord_fine=ones(gridsize*gridsize,2)
    for x in 0:(gridsize-1) 
        for y in 0:(gridsize-1)
            coord_fine[y+x*gridsize+1,2]=x/(gridsize)
            coord_fine[y+x*gridsize+1,1]=y/(gridsize)
        end
    end
    coord_coarse=rand(N_coarse,2).*(gridsize-1)/gridsize
    #coord_all=vcat(coord_fine, coord_coarse)
    #coord_x0=rand(1,2).*(gridsize-1)/gridsize
    x0 = rand(1:gridsize^2)

    if plot
        t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid to simulate on")
        scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations")
        scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="Normalizing Observation")
        title!("Observation and simulation points")
        display(t)
    end
    return coord_coarse, coord_fine, x0
end


# creating the grid of fine points and the coarse observation points
# boolean plot makes plot of grid if true
function Create_Grid_and_Observation_on_fine_grid(gridsize,N_coarse, plot::Bool=false )
    coord_fine=ones(gridsize*gridsize,2)
    for x in 0:(gridsize-1) 
        for y in 0:(gridsize-1)
            coord_fine[y+x*gridsize+1,2]=x/(gridsize)
            coord_fine[y+x*gridsize+1,1]=y/(gridsize)
        end
    end
    coord_sample = sample(1:(gridsize*gridsize), N_coarse+1; replace=false)
    coord_cond_rows = coord_sample[2:end]
    coord_coarse = coord_fine[coord_cond_rows,:]
    x0 = coord_sample[1]

    if plot
        t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid to simulate on")
        scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations")
        scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="Normalizing Observation")
        title!("Observation and simulation points")
        display(t)
    end
    return coord_coarse, coord_fine, x0
end




#get all the nearest gridpoints to the stations, treat the observations as observations there
#functions gives the row indices from the common rows from the reference matrix
function get_common_rows_indices(ReferenceMatrix,ComparisonMatrix)

    # Initialize an empty array to store the indices of common rows
    common_rows_indices = Integer[]
    
    # Iterate through the rows of ReferenceMatrix
    for i in 1:size(ReferenceMatrix, 1)
        # Get the current row of ReferenceMatrix
        row1 = ReferenceMatrix[i, :]
        
        # Iterate through the rows of ComparisonMatrix
        for j in 1:size(ComparisonMatrix, 1)
            # Get the current row of ComparisonMatrix
            row2 = ComparisonMatrix[j, :]
            
            # Check if the rows are equal
            if row1 == row2
                # If they are equal, store the index of the common row
                push!(common_rows_indices, i)
                break  # Break the inner loop since a common row is found
            end
        end
    end
    common_rows_indices
end


#variogram of the process (2 d coordinates, written as vector or matrix)
vario(x,param)=param[1]*sqrt(
                                            (x[1])^2
                                            + (x[2])^2
                                            )^param[2]

vec_vario(param,coord_fine,coord_x0)=param[1].*norm.(eachrow(coord_fine.-coord_x0')).^param[2]


function cov_fun_vario(param,x,y,coord_x0)
    vario(x-coord_x0,param) + vario(y-coord_x0,param)-vario(x-y,param)
end




#l-pareto simulation for creating data from oesting et al
function simu_specfcts(coord, vario_with_params, chol_mat, alpha) 
    N=size(coord,1)
    shift=rand(1:N)
    trend=[vario_with_params(coord[i,:]-coord[shift,:]) for i in 1:N]
    res=chol_mat * rand(Normal(),N)
    res=exp.((res.-res[shift]-trend)) # W_ell (maybe without normalizing constant)
    res=res/mean(res) # W_ell/ l(w_ell)
    res*=(1/(1-rand()))^(1/alpha)
end


#cov and chol matrices for gaussian processes

#cholesky matrix for given vario_with_params and coords (with small numeric constant addition for stability)
function chol_mat(coord, vario_with_params)
    N=size(coord,1)
    cov_mat=ones(N,N)
    for i in 1:N
        for j in 1:N
            cov_mat[i,j]=vario_with_params(coord[i,:])+vario_with_params(coord[j,:])-vario_with_params(coord[i,:]-coord[j,:])
        end
    end
    cov_mat.+=1e-6 
    chol_mat=cholesky(cov_mat).L
end


#Cov Mat 
function cov_mat_for_vectors(coord_x, coord_y, param, coord_x0)
    Nx=size(coord_x,1)
    Ny=size(coord_y,1)
    cov_mat=ones(Nx,Ny)
    for ix in 1:Nx
        for jy in 1:Ny
            cov_mat[ix,jy]=cov_fun_vario(param,coord_x[ix,:], coord_y[jy,:], coord_x0 )  #vario(coord[i,:])+vario(coord[j,:])-vario(coord[i,:]-coord[j,:])
        end
    end
    cov_mat
end






#################################################
#################################################

#create functions for embedding simulation, one for 2dfft, one for modified cov function and one for the simulation, further one to apply simulation and one to apply conditional simulation

#mathlab fft2 analog function
function fft2(A)
    # Apply 1D FFT along each dimension
    fft_rows = fft(A, 1)
    fft_cols = fft(fft_rows, 2)

    # Return the result
    return fft_cols
end

function getFourierMatrix(n)
    F=rand(n,n)*im
    for i in 0:(n-1),j in 0:(n-1)
        F[i+1,j+1]=exp(2*pi*im*i*j/n)
    end
    F   
end

function fft2_hardcode(A)
    F=getFourierMatrix(size(A,1))
    fft_rows = F*A
    fft_cols = (F*fft_rows')'
    # Return the result
    return fft_cols
end


#modified cov function
function rho(x,y,R,alpha)
    #embedding of covariance function on a larger [0,R] × [0,R] Grid
    if alpha<=1.1 #alpha=2 Hurst param
        beta=0
        c_2=alpha/2
        c_0=1-alpha/2
    else
        beta=alpha*(2-alpha)/(3*R*(R^2-1))
        c_2=(alpha-beta*(R-1)^2*(R+2))/2
        c_0=beta*(R-1)^3+1-c_2
    end
    #create cont isotropic cov function
    r=sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
    if r<=1
        out=c_0-r^alpha+c_2*r^2
    elseif R==1 #caution?
        out=0
    elseif r<=R
        out=beta*(R-r)^3/r
    else
        out=0
    end
    return (out,c_0,c_2)
end

function FBM_simu_fast(param,gridsize,num_sim) 
    alpha=param[2]
    c_sqrt=sqrt(param[1])
    H=alpha/2 #Hurst Param
    if alpha<=1.1
        R=1
    else
        R=2 #expanded gridsize, region of interest is [0,1]×[0.1]
    end
    n=m=R*gridsize #size of grid is m ∗ n, cov mat size is n^2*m^2
    tx=(1:n)/n*R; ty=(1:m)/m*R #grid points in x and y 
    Rows=zeros(m,n);
    for i in 1:n 
        for j in 1:m
            Rows[j,i]=rho([tx[i],ty[j]],[tx[1],ty[1]],R,2*H)[1]
        end
    end
    BlkCirc_row=[Rows  Rows[:,end-1:-1:2] ;  Rows[end-1:-1:2,:]  Rows[end-1:-1:2,end-1:-1:2 ]  ]
    #calculate eigenvalues via fft
    eig_vals=real(fft2(BlkCirc_row)/(4*(m-1)*(n-1)))
    #optional:
    #set small values to zero:
    #eps = 10^-8
    #eig_vals=eig_vals.*(abs.(eig_vals) .> eps)
    #eig_vals=eig_vals.*(eig_vals .>= 0)
    eig_vals=sqrt.(eig_vals)
   
    res1=[Matrix{Float64}(undef,gridsize,gridsize) for i in 1:num_sim]
    #one can get two times as many simulations for free by using the imaginary and real part 
    #of the complex gaussian, but they are dependend
    
    #res2=[Matrix{Float64}(undef,gridsize,gridsize) for i in 1:num_sim]
    
    for trial in 1:num_sim
        #generate field with covariance given by block circulant matrix
        Z= randn(2*(m-1),2*(n-1)) + im* randn(2*(m-1),2*(n-1)) 
        #fft to calculate diag*Z
        F=fft2(eig_vals.*Z)
        #extract subblock with desired cov variance
        F=F[1:gridsize,1:gridsize]
        (out,c_0,c_2)=rho([0,0],[0,0],R,2*H)
        #optional two (dependend) real fields
        field1=real(F)
        #field2=imag(F)
        #set field zero at origin
        field1=field1.-field1[1,1]
        #field2=field2.-field2[1,1]
        #make correction for embedding with a term c_2*r^2
        X_grid_comp = [i for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
        Y_grid_comp = [j for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
    
        res1[trial]=field1 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
        #res2[trial]=field2 +  (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
    end

    #(c_sqrt.*res1,c_sqrt.*res2)
    c_sqrt.*res1
end

function FBM_simu_fast_vec(param,gridsize,num_sim) 
    alpha=param[2]
    c_sqrt=sqrt(param[1])
    H=alpha/2 #Hurst Param
    if alpha<=1.1
        R=1
    else
        R=2 #expanded gridsize, region of interest is [0,1]×[0.1]
    end
    n=m=R*gridsize #size of grid is m ∗ n, cov mat size is n^2*m^2
    tx=(1:n)/n*R; ty=(1:m)/m*R #grid points in x and y 
    Rows=zeros(m,n);
    for i in 1:n 
        for j in 1:m
            Rows[j,i]=rho([tx[i],ty[j]],[tx[1],ty[1]],R,2*H)[1]
        end
    end
    BlkCirc_row=[Rows  Rows[:,end-1:-1:2] ;  Rows[end-1:-1:2,:]  Rows[end-1:-1:2,end-1:-1:2 ]  ]
    #calculate eigenvalues via fft
    eig_vals=real(fft2(BlkCirc_row)/(4*(m-1)*(n-1)))
    #optional:
    #set small values to zero:
    #eps = 10^-8
    #eig_vals=eig_vals.*(abs.(eig_vals) .> eps)
    #eig_vals=eig_vals.*(eig_vals .>= 0)
    eig_vals=sqrt.(eig_vals)
    res1=[Vector{Float64}(undef,gridsize*gridsize) for i in 1:num_sim]
    #res2=[Matrix{Float64}(undef,gridsize,gridsize) for i in 1:num_sim]
    
    for trial in 1:num_sim

        
        #generate field with covariance given by block circulant matrix
        Z= randn(2*(m-1),2*(n-1)) + im* randn(2*(m-1),2*(n-1)) 

        #fft to calculate diag*Z
        F=fft2(eig_vals.*Z)
        #extract subblock with desired cov variance
        F=F[1:gridsize,1:gridsize]
        (out,c_0,c_2)=rho([0,0],[0,0],R,2*H)
        #optional two (dependend) real fields
        field1=real(F)
        #field2=imag(F)
        #set field zero at origin
        field1=field1.-field1[1,1]
        #field2=field2.-field2[1,1]
        #make correction for embedding with a term c_2*r^2
        X_grid_comp = [i for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
        Y_grid_comp = [j for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
    
        res1[trial]=vec((field1 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2))')
        #res2[trial]=field2 +  (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
    end

    #(c_sqrt.*res1,c_sqrt.*res2)
    c_sqrt.*res1
end

#gaussian process simulation
function r_log_gaussian(coord_fine,param,row_x0) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    field1 = FBM_simu_fast(param, gridsize,1)[1]
    res = vec(field1')
    res = res - vec_vario(param,coord_fine,coord_fine[row_x0,:]) .-res[row_x0] #variogram
    exp.(res)
    ###           #this is just the semi variogram which is the trend we correct with
end

function r_log_gaussian_vec(coord_fine,param,row_x0,num_rep) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    res = FBM_simu_fast_vec(param, gridsize,num_rep)
    trend=vec_vario(param,coord_fine,coord_fine[row_x0,:])
    for i in 1:num_rep
            res[i] = exp.(res[i] - trend .-res[i][row_x0]) #variogram
    end
    res
end

function r_gaussian_vec(coord_fine,param,row_x0,num_rep) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    res = FBM_simu_fast_vec(param, gridsize,num_rep)
    trend=vec_vario(param,coord_fine,coord_fine[row_x0,:])
    for i in 1:num_rep
            res[i] = (res[i] - trend .-res[i][row_x0]) #variogram
    end
    res
end

function r_cond_log_gaussian(coarse_observation,observation_x0, coord_fine,coord_coarse,param,row_x0) #coord_x0 (hier c egal)
    log_normalized_coarse_observation=log.(coarse_observation./observation_x0)
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
    coord_cond=coord_fine[coord_cond_rows,:]
    
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    sigma_zy= cov_mat_for_vectors(coord_cond, coord_fine,param,coord_fine[row_x0,:])'  #hier
    res=log.(r_log_gaussian(coord_fine,param,row_x0)) 
    res= res .+ sigma_zy*(sigma_yy_inv*(log_normalized_coarse_observation.-res[coord_cond_rows]))
    exp.(res)
end

function r_cond_log_gaussian_vec(coarse_observation,observation_x0, coord_fine,coord_coarse,param,row_x0,num_rep) #coord_x0 (hier c egal)
    #first dim sparse observation, second dim repetition, third dim site
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
    coord_cond=coord_fine[coord_cond_rows,:]
    
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    sigma_zy= cov_mat_for_vectors(coord_cond, coord_fine,param,coord_fine[row_x0,:])'   
    res=[r_gaussian_vec(coord_fine,param,row_x0, num_rep) for j in 1:size(coarse_observation,1)]
    for j in 1:size(coarse_observation,1)
        for i in 1:num_rep
            log_normalized_coarse_observation = log.(coarse_observation[j,:]./observation_x0[j])
            res[j][i] = exp.(res[j][i] .+ sigma_zy*(sigma_yy_inv*(log_normalized_coarse_observation.-res[j][i][coord_cond_rows]))) #variogram
        end
    end
    res
end






# #make conditional simulations for fixed params and then test which realization are risk functional exceedances, use them as data points
# function exceed_cond_sim(num_runs,num_obs,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )
#     res_ell_X=[0.0 for i in 1:num_obs] 
#     old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
#     for trial in 1:num_runs
#         for i in 1:num_obs #direkt num_obs viele simulations
#             if (trial==1)
#                 old_value=r_cond_log_gaussian(observation_data[i,:],observation_x0[i], coord_fine,coord_coarse,param,row_x0)
#             end
#             proposal=r_cond_log_gaussian(observation_data[i,:],observation_x0[i], coord_fine,coord_coarse,param,row_x0)
#             acceptance_rate=min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
#             if (rand()< acceptance_rate)
#                 old_value=proposal
#             end
#             res_ell_X[i]=observation_x0[i]*mean(old_value)
#         end
#     end
#     if sum(res_ell_X.>threshold)==0
#         println("not a single threshold exceedance")
#     else
#         #likelihood calculation and param updates
#         #find all threshold exccedances and calculate the log of them
#         ind=findall(res_ell_X.>threshold)
#         (observation_data[ind,:],observation_x0[ind])
#     end
# end

function exceed_cond_sim(num_runs,num_obs,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )
    tmp = r_cond_log_gaussian_vec(observation_data, observation_x0, coord_fine, coord_coarse,param,row_x0,num_runs+1)
    res_ell_X = [0.0 for i in 1:num_obs] 
    old_value = tmp[1][1]
   #old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
    for trial in 1:num_runs
        for i in 1:num_obs #direkt num_obs viele simulations
            # if (trial==1)
            #         old_value = tmp[i][trial+1]
            # end
            proposal = tmp[i][trial+1]
            acceptance_rate = min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
            res_ell_X[i]=observation_x0[i]*mean(old_value)
        end
    end
    if sum(res_ell_X.>threshold)==0
        println("not a single threshold exceedance")
    else
        #likelihood calculation and param updates
        #find all threshold exccedances and calculate the log of them
        ind=findall(res_ell_X.>threshold)
        (observation_data[ind,:],observation_x0[ind])
    end
end




#likelihood functions for parameter updates

#first gaussian density
function log_d_gaussian(trend,cov_mat_inv,x,inv_determinant)
    (-0.5*transpose(x-trend) * cov_mat_inv * (x-trend) )+0.5*log(inv_determinant)#+log(sqrt((2*pi)^(-N)))
end

#log likelihood of log gaussian density of observations 
function l_1_fun(coord_fine,coord_coarse,coarse_observation,param, observation_x0, row_x0)
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.* gridsize)./gridsize)
    cov_mat_coarse_inv= inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    inv_determinant = det(cov_mat_coarse_inv)
    trend = -vec_vario(param,coord_fine[coord_cond_rows,:],coord_fine[row_x0,:])
    sum([(log_d_gaussian(trend ,cov_mat_coarse_inv , log.(coarse_observation[i,:]./observation_x0[i]), inv_determinant)) for i in 1:size(coarse_observation,1)])
end



#reciprocal mean estimation
# one estimate of c is: mean(exp.(r_gaussian(res_fine["mu"],res_fine["chol_mat"])))^(-alpha)
#therefore we N_est_c many of them and take there mean, in addition we multiply with the number of threshold exccedances
# function l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
#     -number_of_exceed * log(mean([mean(r_log_gaussian(coord_fine,param,row_x0 )  
#         )^(alpha) for i in 1: N_est_c]))  
#      # minus for 1/c_l (in log)
#  end

 function l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
    tmp = r_log_gaussian_vec(coord_fine,param,row_x0, N_est_c) 
    -number_of_exceed * log(mean([mean(tmp[i] 
        )^(alpha) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end



#  #conditional mean estimator for intergal part #add param zu l_3
# function l_3_fun(coord_fine, coord_coarse, param, row_x0, coarse_observation, observation_x0, alpha, N_est_cond)   
#     res_cond_int = zeros(N_est_cond,size(coarse_observation,1))
#     for row in 1:N_est_cond
#         res_cond_int[row,:]=[(mean(
#                                 r_cond_log_gaussian(coarse_observation[i,:],observation_x0[i], coord_fine, coord_coarse, param, row_x0)    
#                                 )^alpha) for i in 1:size(coarse_observation,1) ]
#     end
#     sum(log.(mean( res_cond_int, dims=1)))
# end


 #conditional mean estimator for intergal part #add param zu l_3
function l_3_fun(coord_fine, coord_coarse, param, row_x0, coarse_observation, observation_x0, alpha, N_est_cond)   
    tmp = r_cond_log_gaussian_vec(coarse_observation,observation_x0,coord_fine, coord_coarse, param, row_x0, N_est_cond)
    res_cond_int = zeros(N_est_cond,size(coarse_observation,1))
    for row in 1:N_est_cond
        res_cond_int[row,:]=[(mean(
                                 tmp[i][row]
                                )^alpha) for i in 1:size(coarse_observation,1) ]
    end
    sum(log.(mean( res_cond_int, dims=1)))
end


#r_cond_log_gaussian(coarse_observation,observation_x0, coord_fine,coord_coarse,param,row_x0) 

#this function proposes a new value uniformly distributed in a 2 eps window in between min_val and max_val arround the old_param
function uniform_proposal(old_param,eps,min_val,max_val)
    if (old_param>min_val+eps && old_param<max_val-eps)
        return rand(Uniform(old_param-eps,old_param+eps))
    elseif (old_param<min_val+eps)
        return rand(Uniform(min_val,min_val+2*eps))
    else
        return rand(Uniform(max_val-2*eps,max_val))
    end
end





#this functions calculates the acceotance rate according to the likelihoods and updates the parameter accordingly
function parameter_update(param_old,param_new,log_likelihood_old,log_likelihood_new)
    #calculate MCMC acceptance rate a
    if log_likelihood_new<log_likelihood_old
        a=exp(log_likelihood_new-log_likelihood_old)
    else
        a=1
    end
    rand_sample=rand()
    #accept new value with probability set as acceptance rate a
    if (rand_sample<a)
        #we update param aka return the new ones
        return(param_new,log_likelihood_new)
    else
        #else we just keep the old values and return them
        return(param_old,log_likelihood_old)
    end
end




#MCMC Algorithm
function MCMC(N_MCMC,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0,n_trial_print,N_est_c,N_est_cond)

    num_obs=size(observation_data,1)

    par_beta_vec  = repeat([param[2]],N_MCMC+1)
    par_c_vec = [param[1] for i=1:N_MCMC+1]
    par_alpha_vec = [alpha for i=1:N_MCMC+1]
    for trial in 1:N_MCMC
        #println("current trial: $trial and current param: $param")
        param[2]=par_beta_vec[trial]
        param[1]=par_c_vec[trial]
        alpha=par_alpha_vec[trial]
        par_beta_old=par_beta_vec[trial]
        par_c_old=par_c_vec[trial]
        #fixed param
        (modified_observation, modified_observation_x0) = exceed_cond_sim(5,num_sim,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )

        #fixed param
        #gridsize = Int(sqrt(size(coord_fine,1)))
        #coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0)
        l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l3=l_3_fun(coord_fine, coord_coarse, param, row_x0, modified_observation, modified_observation_x0, alpha, N_est_cond)   
        log_likelihood_old=sum([l1,l2,l3])
        #new param c

        #new param beta
        par_beta_old=param[2]
        #propose new param beta
        eps_beta=0.05
        param[2]=uniform_proposal(par_beta_old,eps_beta,0.0,2.0)
        par_beta=param[2]

        #calculate new log likelihood
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0)
        l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l3=l_3_fun(coord_fine, coord_coarse, param, row_x0, modified_observation, modified_observation_x0, alpha, N_est_cond)   
        log_likelihood_new=sum([l1,l2,l3])  

        #MCMC update of param according to acceptance rate calculated with old and new likelihoods
        param[2],log_likelihood_old=parameter_update(par_beta_old,par_beta,log_likelihood_old,log_likelihood_new)
        #safe param after MCMC step
        par_beta_vec[trial+1]=param[2]
        
        #new param c
        par_c_old=param[1]
        #propose new param c
        eps_c=0.1
        param[1]=uniform_proposal(par_c_old,eps_c,0.0,10.0)
        par_c=param[1]
        #calculate new log likelihood
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0)
        l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l3=l_3_fun(coord_fine, coord_coarse, param, row_x0, modified_observation, modified_observation_x0, alpha, N_est_cond)   
        log_likelihood_new=sum([l1,l2,l3])  
        #MCMC update of param according to acceptance rate calculated with old and new likelihoods
        param[1],log_likelihood_old=parameter_update(par_c_old,par_c,log_likelihood_old,log_likelihood_new)
        #safe param c after MCMC step
        par_c_vec[trial+1]=param[1]     
        
        
        #new param alpha
        #calculate old log likelihood
        number_of_exceed=size(modified_observation,1)
        log_likelihood_old_alpha=number_of_exceed*log(alpha) + l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) + (-alpha-1)*sum(log.(modified_observation_x0./threshold))
        alpha_old=alpha
        #propose new alpha
        eps_alpha=0.05
        alpha=uniform_proposal(alpha_old,eps_alpha,0.0,10.0)
        #calculate new log likelihood
        log_likelihood_new_alpha=number_of_exceed*log(alpha) + l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c) + (-alpha-1)*sum(log.(modified_observation_x0./threshold))
        alpha,log_likelihood_old_alpha=parameter_update(alpha_old,alpha,log_likelihood_old_alpha,log_likelihood_new_alpha)
        #safe param c after MCMC step
        par_alpha_vec[trial+1]=alpha 
        #just print every n_trial_print's trial number to see how fast programme runs
        if trial%n_trial_print==0
            println(trial)
        end
    end
    Dict( "beta" => par_beta_vec, "c" => par_c_vec, "alpha" => par_alpha_vec)
end  