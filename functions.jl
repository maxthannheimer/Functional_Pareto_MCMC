using LinearAlgebra
#import Pkg; Pkg.add("Distributions")
using Random, Distributions
using Plots
using FFTW
using JLD2
using Printf

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
        t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid Simulations",color=:red,alpha=0.4)
        scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations",color=:purple)
        scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="",color=:purple)
        title!("Observation and Simulation Points")
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
        t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid Simulations",color=:red,alpha=0.4)
        scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations",color=:purple)
        scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="", color=:purple)
        title!("Observation and Simulation Points")
        display(t)
    end
    return coord_coarse, coord_fine, x0
end



#create a 7 x 7 grid with 9 points as observations nicely distributed
function Create_Grid_and_Observation_seven_times_seven( plot::Bool=false )
    gridsize=7
    coord_fine=ones(gridsize*gridsize,2)
    for x in 0:(gridsize-1) 
        for y in 0:(gridsize-1)
            coord_fine[y+x*gridsize+1,2]=x/(gridsize)
            coord_fine[y+x*gridsize+1,1]=y/(gridsize)
        end
    end
    #coord_sample = sample(1:(gridsize*gridsize), N_coarse+1; replace=false)
    coord_cond_rows = [9,11,13,23,27,37,39,41]
    coord_coarse = coord_fine[coord_cond_rows,:]
    x0 = 25 # index of the normalization point

    if plot
        t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid Simulations",color=:red,alpha=0.4)
        scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations",color=:purple)
        scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="", color=:purple)
        title!("Observation and Simulation Points")
        display(t)
    end
    return coord_coarse, coord_fine, x0
end




#create a 9 x 9 grid with 9 points as observations nicely distributed
function Create_Grid_and_Observation_9x9( plot::Bool=false )
    gridsize=9
    coord_fine=ones(gridsize*gridsize,2)
    for x in 0:(gridsize-1) 
        for y in 0:(gridsize-1)
            coord_fine[y+x*gridsize+1,2]=x/(gridsize)
            coord_fine[y+x*gridsize+1,1]=y/(gridsize)
        end
    end
    #coord_sample = sample(1:(gridsize*gridsize), N_coarse+1; replace=false)
    coord_cond_rows = [11,14,17,38,44,65,68,71]
    coord_coarse = coord_fine[coord_cond_rows,:]
    x0 = 41 # index of the normalization point

    if plot
        t=scatter(coord_fine[:,1],coord_fine[:,2], label="Fine Grid Simulations",color=:red,alpha=0.4)
        scatter!(coord_coarse[:,1],coord_coarse[:,2],label="Coarse Observations",color=:purple)
        scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="", color=:purple)
        title!("Observation and Simulation Points")
        display(t)
    end
    return coord_coarse, coord_fine, x0
end


#9x9 grid, 9 points as observations nicely distributed, just plot the grid, points are translated to be centered in the unit square, just for illustration
function Create_Grid_and_Observation_9x9_translated_plot()
    gridsize=9
    coord_fine=ones(gridsize*gridsize,2)
    for x in 0:(gridsize-1) 
        for y in 0:(gridsize-1)
            coord_fine[y+x*gridsize+1,2]=x/(gridsize)
            coord_fine[y+x*gridsize+1,1]=y/(gridsize)
        end
    end
    #coord_sample = sample(1:(gridsize*gridsize), N_coarse+1; replace=false)
    coord_fine=coord_fine.+1/18 #translate the grid a bit to see the points better
    coord_cond_rows = [11,14,17,38,44,65,68,71]
    coord_coarse = coord_fine[coord_cond_rows,:]
    x0 = 41 # index of the normalization point
   # Define tick positions and labels for 0, 1/3, 2/3, 1
    ticks = [0.0, 1/3, 2/3, 1.0]
    ticklabels = ["0", "1/3", "2/3", "1"]

    t = scatter(
        coord_fine[:,1], coord_fine[:,2],
        color=:red, alpha=0.4,
        ylims=(0,1), xlims=(0,1),
        xticks=(ticks, ticklabels),
        yticks=(ticks, ticklabels),
        label=""
    )
    scatter!(coord_coarse[:,1],coord_coarse[:,2],color=:purple,label="")
    scatter!([coord_fine[x0,1]],[coord_fine[x0,2]],label="", color=:purple)
    display(t)
    
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


#variogram of the process 
# c⋅||x||^β
vario(x,param)=param[1]*sqrt(
                                            (x[1])^2
                                            + (x[2])^2
                                            )^param[2]
#variogram of the process (2 d coordinates, written as vector or matrix)
# c⋅||x-x0||^β for each location in coord_fine
vec_vario(param,coord_fine,coord_x0)=param[1].*norm.(eachrow(coord_fine.-coord_x0')).^param[2]

#cov_function for two locations x and y 
# c⋅||x-x0||^β + c⋅||y-x0||^β -c⋅||x-y||^β 
function cov_fun_vario(param,x,y,coord_x0)
    vario(x-coord_x0,param) + vario(y-coord_x0,param)-vario(x-y,param)
end




#thre different l-pareto simulation approaches for creating data from oesting et al, not sure if valied for α≠1
function simu_specfcts_new(coord, vario_with_params, chol_mat, alpha) 
    N=size(coord,1)
    shift=rand(1:N)
    trend=[vario_with_params(coord[i,:]-coord[shift,:]) for i in 1:N]
    res=chol_mat * rand(Normal(),N)
    res=exp.(res.-res[shift]-trend) # W_ell (maybe without normalizing constant)
    res=res/mean(res) # W_ell/ l(w_ell)
    res*=(1/(1-rand()))
    res.^(1/alpha)
end

function simu_specfcts_old(coord, vario_with_params, chol_mat, alpha) 
    N=size(coord,1)
    shift=rand(1:N)
    trend=[vario_with_params(coord[i,:]-coord[shift,:]) for i in 1:N]
    res=chol_mat * rand(Normal(),N)
    res=exp.(1/alpha*(res.-res[shift]-trend)) # W_ell (maybe without normalizing constant)
    res=res/mean(res) # W_ell/ l(w_ell)
    res*=(1/(1-rand()))^(1/alpha)
end

function simu_specfcts_verynew(coord, vario_with_params, chol_mat, alpha)
    N=size(coord,1)
    shift=rand(1:N)
    trend=[vario_with_params(coord[i,:]-coord[shift,:]) for i in 1:N]
    res=chol_mat * rand(Normal(),N)
    res=exp.(1/alpha.*(res.-res[shift]-trend)) # W_ell (maybe without normalizing constant)
    res=res/(sum(res.^alpha)^(1/alpha)) # W_ell/ l(w_ell) for l=alpha-norm
    res*=(1/(1-rand()))^(1/alpha) / N^(1-alpha)
    while(sum(res) < 1)    # acceptance-rejection to obtain l=1-norm
        shift=rand(1:N)
        trend=[vario_with_params(coord[i,:]-coord[shift,:]) for i in 1:N]
        res=chol_mat * rand(Normal(),N)
        res=exp.(1/alpha.*(res.-res[shift]-trend))
        res=res/(sum(res.^alpha)^(1/alpha))
        res*=(1/(1-rand()))^(1/alpha) / N^(1-alpha)
    end
    return(res/mean(res))
end 


#using MCMC approach from dombry et all 2015 to simulate chain with pareto process as stationary distrib
function simu_specfcts_MCMC(num_runs, alpha, coord_fine,param,row_x0 )
    #sample num_runs+1 many realizations of exp(1/α[G(s)-G(x0)-γ(s-x0)])
    tmp = r_log_gaussian_vec(coord_fine,param,row_x0,num_runs+1,alpha) 
    old_value = tmp[1]
    #use samples in MCMC approach to converge to stationary distrib W^(r)
    for trial in 1:num_runs
            proposal = tmp[trial+1]
            acceptance_rate = min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
    end
    #normalize sample and multiply pareto intesity to get pareto process sample
    old_value = old_value/mean(old_value)
    old_value*=(1/(1-rand()))^(1/alpha)
    #proposal=proposal*(1/(1-rand()))^(1/alpha)
end

#get threshold empirically:
function quantile_threshold(num_runs,alpha,coord_fine,param,row_x0,threshold_quantile,num_sim)
    #calculate threshold for each simulation
    sim_data= [simu_specfcts_MCMC(num_runs, alpha, coord_fine,param,row_x0 ) for i in 1:num_sim]
    sim_data=reduce(hcat,sim_data)'
    risk_functional_evaluation=[mean(sim_data[I,:]) for I in 1:num_sim]
    #calculate threshold for each simulation
    quantile(sort(risk_functional_evaluation),threshold_quantile)
end

#cov and chol matrices for gaussian processes

#cholesky matrix 
#calculates cholesky matrix without normalization for given coordinates and the distance based variogram with parameters
# (with small numeric constant addition for stability)
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
#calculate cov matrix for to matrices of coordinates and the variogramm given via param and normalized at x0
function cov_mat_for_vectors(coord_x, coord_y, param, coord_x0)
    Nx=size(coord_x,1)
    Ny=size(coord_y,1)
    cov_mat=ones(Nx,Ny)
    for ix in 1:Nx
        for jy in 1:Ny
            cov_mat[ix,jy]=cov_fun_vario(param,coord_x[ix,:], coord_y[jy,:], coord_x0 )  
            #second implementation to check if everything works
            #cov_mat[ix,jy]=vario(coord_x[ix,:]-coord_x0,param)+vario(coord_y[jy,:]-coord_x0,param)-vario(coord_x[ix,:]-coord_y[jy,:],param)
        end
    end
    cov_mat
end






#################################################
#################################################

#create functions for embedding simulation, one for 2dfft, one for modified cov function and one for the simulation, further one to apply simulation and one to apply conditional simulation

#matlab fft2 analog function
function fft2(A)
    # Apply 1D FFT along each dimension
    fft_rows = fft(A, 1)
    fft_cols = fft(fft_rows, 2)

    # Return the result
    return fft_cols
end

#alternative implementation of fft2 for comparision
#frist create FourierMatrix
function getFourierMatrix(n)
    F=rand(n,n)*im
    for i in 0:(n-1),j in 0:(n-1)
        F[i+1,j+1]=exp(2*pi*im*i*j/n)
    end
    F   
end

#Then peform Fourier Transform
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


function FBM_simu_fast_vec_dependent(param,gridsize,num_sim) 
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
    if isodd(num_sim)
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
    
        res1[end]=vec((field1 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2))') 
    end
    if num_sim==1
        return c_sqrt.*res1
    end  

    for trial in 1:(num_sim÷2)

        
        #generate field with covariance given by block circulant matrix
        Z= randn(2*(m-1),2*(n-1)) + im* randn(2*(m-1),2*(n-1)) 

        #fft to calculate diag*Z
        F=fft2(eig_vals.*Z)
        #extract subblock with desired cov variance
        F=F[1:gridsize,1:gridsize]
        (out,c_0,c_2)=rho([0,0],[0,0],R,2*H)
        #optional two (dependend) real fields
        field1=real(F)
        field2=imag(F)
        #set field zero at origin
        field1=field1.-field1[1,1]
        field2=field2.-field2[1,1]
        #make correction for embedding with a term c_2*r^2
        X_grid_comp = [i for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
        Y_grid_comp = [j for i in tx[1:gridsize], j in tx[1:gridsize]].-tx[1]
    
        res1[trial]=vec((field1 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2))')
        res1[trial+num_sim÷2]=vec((field2 + (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2))')
        #res2[trial]=field2 +  (X_grid_comp*randn()+Y_grid_comp*randn() ) .*sqrt(2*c_2)
    end
    #(c_sqrt.*res1,c_sqrt.*res2)

    
    c_sqrt.*res1
end

#gaussian process simulation
#simulate exp(1/α[G(s)-G(x0)-γ(s-x0)])
function r_log_gaussian(coord_fine,param,row_x0,alpha) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    field1 = FBM_simu_fast(param, gridsize,1)[1]
    res = vec(field1')
    res = res - vec_vario(param,coord_fine,coord_fine[row_x0,:]) .-res[row_x0] #variogram
    exp.(1/alpha*res)
    ###           #this is just the semi variogram which is the trend we correct with
end
#simulate numrep many exp(1/α[G(s)-G(x0)-γ(s-x0)])
function r_log_gaussian_vec(coord_fine,param,row_x0,num_rep,alpha) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    res = FBM_simu_fast_vec(param, gridsize,num_rep)
    trend=vec_vario(param,coord_fine,coord_fine[row_x0,:])
    for i in 1:num_rep
            res[i] = exp.(1/alpha*(res[i] - trend .-res[i][row_x0])) #variogram
    end
    res
end
#simulate numrep many exp(1/α[G(s)-G(x0)-γ(s-x0)])
function r_log_gaussian_vec_dependent(coord_fine,param,row_x0,num_rep,alpha) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    res = FBM_simu_fast_vec_dependent(param, gridsize,num_rep)
    trend=vec_vario(param,coord_fine,coord_fine[row_x0,:])
    for i in 1:num_rep
            res[i] = exp.(1/alpha*(res[i] - trend .-res[i][row_x0])) #variogram
    end
    res
end

#simulate numrep many 1/α [G(s)-G(x0)-γ(s-x0)] 
function r_gaussian_vec(coord_fine,param,row_x0,num_rep,alpha) 
    gridsize = Int(sqrt(size(coord_fine,1)))
    res = FBM_simu_fast_vec(param, gridsize,num_rep)
    trend=vec_vario(param,coord_fine,coord_fine[row_x0,:])
    for i in 1:num_rep
            res[i] = 1/alpha*(res[i] - trend .-res[i][row_x0]) #variogram
    end
    res
end

#simulate numrep many 1/α [G(s)-G(x0)-γ(s-x0)] on coarse grid
function r_gaussian_vec_coarse(coord_coarse_plus,coord_x0,param,num_rep,alpha) 
    N = size(coord_coarse_plus,1)
    cov_mat = cov_mat_for_vectors(coord_coarse_plus, coord_coarse_plus,param,coord_x0).+1e-6 
    res = rand(MvNormal([0.0 for i in 1:N],cov_mat),num_rep)'
    trend=vec_vario(param,coord_coarse_plus,coord_coarse_plus[end,:])
    for i in 1:num_rep
            res[i,:] = 1/alpha*(res[i,:] - trend .-res[i,end]) #variogram
    end
    res
end

#conditional simulation of exp(1/α[G(s)-G(x0)-γ(s-x0)]) | observation_data(i)./observation_x0(i)
function r_cond_log_gaussian(coarse_observation,observation_x0, coord_fine,coord_coarse,param,row_x0,alpha) #coord_x0 (hier c egal)
    log_normalized_coarse_observation=log.(coarse_observation./observation_x0)
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
    coord_cond=coord_fine[coord_cond_rows,:]
    
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    sigma_zy= cov_mat_for_vectors(coord_cond, coord_fine,param,coord_fine[row_x0,:])'  #hier
    res=log.(r_log_gaussian(coord_fine,param,row_x0,alpha)) 
    res= res + sigma_zy*(sigma_yy_inv*(log_normalized_coarse_observation-res[coord_cond_rows]))
    exp.(res)
end

function r_cond_log_gaussian_vec(coarse_observation,observation_x0, coord_fine,coord_coarse,param,row_x0,num_rep,alpha) #coord_x0 (hier c egal)
    #first dim number of simulated or observed data repetitions, second dim num_rep (how many simulations are wanted), third dim site in fine grid
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
    coord_cond=coord_fine[coord_cond_rows,:]
    
    sigma_yy_inv = inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    sigma_zy= cov_mat_for_vectors(coord_cond, coord_fine,param,coord_fine[row_x0,:])'   
    res=[r_gaussian_vec(coord_fine,param,row_x0,num_rep,alpha) for j in 1:size(coarse_observation,1)]
    for j in 1:size(coarse_observation,1)
        for i in 1:num_rep
            log_normalized_coarse_observation = log.(coarse_observation[j,:]./observation_x0[j])
            res[j][i] = exp.(res[j][i] + sigma_zy*(sigma_yy_inv*(log_normalized_coarse_observation-res[j][i][coord_cond_rows]))) #variogram
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
    tmp = r_cond_log_gaussian_vec(observation_data, observation_x0, coord_fine, coord_coarse,param,row_x0,num_runs+1,alpha)
    res_ell_X = [0.0 for i in 1:num_obs] 

   #old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
   for i in 1:num_obs #direkt num_obs viele simulations 
        old_value = tmp[i][1]
        for trial in 1:num_runs
            proposal = tmp[i][trial+1]
            acceptance_rate = min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
        end
        res_ell_X[i]=observation_x0[i]*mean(old_value)
    end
    if sum(res_ell_X.>threshold)==0
        println("not a single threshold exceedance:")
        println("parameter: $param, alpha: $alpha")
        Base._throw_argerror("not a single threshold exceedance")
    else
        #likelihood calculation and param updates
        #find all threshold exccedances and calculate the log of them
        ind=findall(res_ell_X.>threshold)
        (observation_data[ind,:],observation_x0[ind])
    end
end

function exceed_cond_sim_approx(num_obs,observation_data,observation_x0,threshold )
#    tmp = r_cond_log_gaussian_vec(observation_data, observation_x0, coord_fine, coord_coarse,param,row_x0,num_runs+1,alpha)
    res_ell_X = [0.0 for i in 1:num_obs] 

   #old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
   for i in 1:num_obs #direkt num_obs viele simulations 
        #old_value = tmp[i][1]
        #for trial in 1:num_runs
         #   proposal = tmp[i][trial+1]
          #  acceptance_rate = min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
           # if (rand()< acceptance_rate)
             #   old_value=proposal
            #end
        #end
        res_ell_X[i]=mean(observation_data[i,:]) #add x0
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

function exceed_cond_sim_quantile(num_runs,num_obs,observation_data,observation_x0,threshold_quantile, alpha, coord_fine,coord_coarse,param,row_x0 )
    tmp = r_cond_log_gaussian_vec(observation_data, observation_x0, coord_fine, coord_coarse,param,row_x0,num_runs+1,alpha)
    res_ell_X = [0.0 for i in 1:num_obs] 
   #old_value=r_cond_log_gaussian(observation_data[1,:],observation_x0[1], coord_fine,coord_coarse,param,row_x0)
   for i in 1:num_obs #direkt num_obs viele simulations
        old_value = tmp[i][1]    
        for trial in 1:num_runs
            proposal = tmp[i][trial+1]
            acceptance_rate = min(1,mean(proposal)^alpha/mean(old_value)^alpha)   
            if (rand()< acceptance_rate)
                old_value=proposal
            end
        end
        res_ell_X[i]=observation_x0[i]*mean(old_value)
    end

        #likelihood calculation and param updates
        #find all threshold exccedances and calculate the log of them
#define threshold as quantile
        threshold = quantile(sort(res_ell_X),threshold_quantile)
        ind = findall(res_ell_X.> threshold)
        #(observation_data,observation_x0,res_ell_X)
        (observation_data[ind,:],observation_x0[ind],threshold)
end



#likelihood functions for parameter updates

#first gaussian density
function log_gauss_1d(x,mu,sigma)
    exp(-0.5*(log(x)-mu)^2/sigma^2)/(x*sqrt(2*pi)*sigma)
end

function log_likehood_log_gauss_1d(x,mu,sigma)
    (-0.5*(log(x)-mu)^2/sigma^2)-log(x*sqrt(2*pi)*sigma)
end

function log_likehood_log_gauss_1d_cut(x,mu,sigma,cut_value)
    if x>cut_value
         @warn "should not happen, parameter proposal is over $cut_value, please resample"
        return -Inf
    else
    (-0.5*(log(x)-mu)^2/sigma^2)-log(x*sqrt(2*pi)*sigma)
    end
end


function log_d_gaussian(trend,cov_mat_inv,x,inv_determinant)
    (-0.5*transpose(x-trend) * cov_mat_inv * (x-trend) )+0.5*log(inv_determinant)#+log(sqrt((2*pi)^(-N)))
end

function d_gaussian(trend,cov_mat_inv,x,inv_determinant)
    sqrt(det(inv_determinant))*exp(-1/2*(transpose(x-trend)*cov_mat_inv*(x-trend)))#*sqrt((2*pi)^(-length(trend)))
    
end

#log likelihood of log gaussian density of observations 
function l_1_fun(coord_fine,coord_coarse,coarse_observation,param, observation_x0, row_x0,alpha)
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.* gridsize)./gridsize)
    cov_mat_coarse_inv= inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    inv_determinant = det(cov_mat_coarse_inv)
    trend = -vec_vario(param,coord_fine[coord_cond_rows,:],coord_fine[row_x0,:])
    sum([(log_d_gaussian(trend ,cov_mat_coarse_inv , alpha*log.(coarse_observation[i,:]./observation_x0[i]), inv_determinant))+log(alpha)*length(trend) for i in 1:size(coarse_observation,1)])
end

function l_1_fun_approx(coord_coarse_plus,coarse_observation,param, observation_x0,alpha)
    #N= size(coord_coarse_plus,1)
    #trend=vec_vario(param,coord_coarse_plus,coord_coarse_plus[end,:])
    trend = -vec_vario(param,coord_coarse_plus[1:(end-1),:],coord_coarse_plus[end,:])
    cov_mat_coarse_inv= inv(cov_mat_for_vectors(coord_coarse_plus[1:(end-1),:],coord_coarse_plus[1:(end-1),:],  param, coord_coarse_plus[end,:])) #hier 
    inv_determinant = det(cov_mat_coarse_inv)
    #trend = -vec_vario(param,coord_fine[coord_cond_rows,:],coord_fine[row_x0,:])
    sum([(log_d_gaussian(trend ,cov_mat_coarse_inv , alpha*log.(coarse_observation[i,:]./observation_x0[i]), inv_determinant))+log(alpha)*length(trend) for i in 1:size(coarse_observation,1)])
end


function l_1_fun_alternative(coord_fine,coord_coarse,coarse_observation,param, observation_x0, row_x0, alpha)
    gridsize = Int(sqrt(size(coord_fine,1)))
    coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.* gridsize)./gridsize)
    cov_mat_coarse_inv= inv(cov_mat_for_vectors(coord_fine[coord_cond_rows,:],coord_fine[coord_cond_rows,:],  param, coord_fine[row_x0,:])) #hier 
    inv_determinant = det(cov_mat_coarse_inv)
    trend = -vec_vario(param,coord_fine[coord_cond_rows,:],coord_fine[row_x0,:])
    #sum([log(d_gaussian(trend ,cov_mat_coarse_inv , alpha*log.(coarse_observation[i,:]./observation_x0[i]), inv_determinant)*alpha^length(trend)) for i in 1:size(coarse_observation,1)])
    sum([log(d_gaussian(trend ,cov_mat_coarse_inv , alpha*log.(coarse_observation[i,:]./observation_x0[i]), inv_determinant))+log(alpha)*length(trend) for i in 1:size(coarse_observation,1)])

end


#reciprocal mean estimation
# one estimate of c is: mean(exp.(r_gaussian(res_fine["mu"],res_fine["chol_mat"])))^(-alpha)
#therefore we N_est_c many of them and take there mean, in addition we multiply with the number of threshold exccedances
# function l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
#     -number_of_exceed * log(mean([mean(r_log_gaussian(coord_fine,param,row_x0 )  
#         )^(alpha) for i in 1: N_est_c]))  
#      # minus for 1/c_l (in log)
#  end

    
function l_2_fun_approx(coord_coarse_plus,coord_x0,param,number_of_exceed,alpha,N_est_c)#(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
    tmp= exp.(r_gaussian_vec_coarse(coord_coarse_plus,coord_x0,param,N_est_c,alpha)) 
    #tmp = r_log_gaussian_vec(coord_fine,param,row_x0, N_est_c,alpha) 
    -number_of_exceed * log(mean([mean(tmp[i,:] 
        )^(alpha) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end

 function l_2_fun(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
    tmp = r_log_gaussian_vec(coord_fine,param,row_x0, N_est_c,alpha) 
    -number_of_exceed * log(mean([mean(tmp[i] 
        )^(alpha) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end

 

 function c_estimation(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
    tmp = r_log_gaussian_vec(coord_fine,param,row_x0, N_est_c,alpha) 
    1/(mean([mean(tmp[i] 
        )^(alpha) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end

 function l_2_fun_dependent(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
    tmp = r_log_gaussian_vec_dependent(coord_fine,param,row_x0, N_est_c,alpha) 
    -number_of_exceed * log(mean([mean(tmp[i] 
        )^(alpha) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
 end

function l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha,N_est_c)
    tmp = r_log_gaussian_vec_dependent(coord_fine,param,row_x0, N_est_c,alpha) 
    -log(mean([mean( tmp[i] 
        )^(alpha) for i in 1:N_est_c]))  
     # minus for 1/c_l (in log)
end

 function c_estimation_dependent(coord_fine, param,row_x0, number_of_exceed,alpha,N_est_c)
    tmp = r_log_gaussian_vec_dependent(coord_fine,param,row_x0, N_est_c,alpha) 
    1/(mean([mean(tmp[i] 
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
#= function l_3_fun_old(coord_fine, coord_coarse, param, row_x0, coarse_observation, observation_x0, alpha, N_est_cond)   
    tmp = r_cond_log_gaussian_vec(coarse_observation,observation_x0,coord_fine, coord_coarse, param, row_x0, N_est_cond)
    res_cond_int = zeros(N_est_cond,size(coarse_observation,1))
    for row in 1:N_est_cond
        res_cond_int[row,:]=[(mean(
                                 tmp[i][row]
                                )^alpha) for i in 1:size(coarse_observation,1) ]
    end
    sum(log.(mean( res_cond_int, dims=1)))
end =#

function l_3_fun(observation_x0, alpha, threshold)   
    #sum([log(alpha*(observation_x0[i]/threshold)^(-alpha-1)) for i in 1:size(observation_x0,1)])
    sum([log(alpha)+log((observation_x0[i]/threshold))*(-alpha-1) for i in 1:size(observation_x0,1)])
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

function gaussian_proposal(old_param,eps)
    return exp(rand(Normal(log(old_param),eps)))
end

function gaussian_proposal_cut(old_param,eps,cut_val)
    while true
        val=exp(rand(Normal(log(old_param),eps)))
        if (val<cut_val)
            return val
        end
        @warn "should not happen, gaussian proposal is over $cut_val, please resample"  
    end
end


#this functions calculates the acceotance rate according to the likelihoods and updates the parameter accordingly
function parameter_update_old(param_old,param_new,log_likelihood_old,log_likelihood_new)
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
function parameter_update(beta_old,c_old,alpha_old,beta_new,c_new,alpha_new,log_likelihood_old,log_likelihood_new)
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
        return(beta_new,c_new,alpha_new,log_likelihood_new)
    else
        #else we just keep the old values and return them
        return(beta_old,c_old,alpha_old,log_likelihood_old)
    end
end
function parameter_update(beta_old,c_old,alpha_old,beta_new,c_new,alpha_new,log_likelihood_old,log_likelihood_new,l2_no_exceed_old,l2_no_exceed_new)
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
        return(beta_new,c_new,alpha_new,log_likelihood_new,l2_no_exceed_new)
    else
        #else we just keep the old values and return them
        return(beta_old,c_old,alpha_old,log_likelihood_old,l2_no_exceed_old)
    end
end


# MCMC Algorithm with Dependent 1/c estimation
function MCMC(N_MCMC,observation_data,observation_x0,threshold,threshold_method, alpha, coord_fine,coord_coarse,param,row_x0,n_trial_print,N_est_c,N_cond_sim)
    if threshold_method=="quantile"
        threshold_quantile = threshold
    elseif !(threshold_method=="fixed" || threshold_method=="all_exceed")
        println("threshold_method not defined, please use fixed or quantile or all_exceed")
        Base._throw_argerror("threshold_method not defined, please use fixed or quantile or all_exceed")
    end
    num_obs=size(observation_data,1)

    par_beta_vec  = repeat([param[2]],N_MCMC+1)
    par_c_vec = [param[1] for i=1:N_MCMC+1]
    par_alpha_vec = [alpha for i=1:N_MCMC+1]
    threshold_vec = [threshold for i=1:N_MCMC]
    number_exceed_vec = [-1 for i=1:N_MCMC]
    log_likelihood_vec = [0.0 for i=1:N_MCMC]
    for trial in 1:N_MCMC
        #println("current trial: $trial and current param: $param")
        param[2]=par_beta_vec[trial]
        param[1]=par_c_vec[trial]
        alpha=par_alpha_vec[trial]
        
        #quantile threshold, changes threshold every trial with fixed quantile threshold_quantile
        if threshold_method=="quantile"
            (modified_observation, modified_observation_x0,threshold) = exceed_cond_sim_quantile(N_cond_sim,num_obs,observation_data,observation_x0,threshold_quantile, alpha, coord_fine,coord_coarse,param,row_x0 )
        #fixed threshold, threshold exceedances are evaluated every step with changing parameters
        elseif threshold_method=="fixed"
            (modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )
        #fixed threshold, everything is an eceedance
        elseif threshold_method=="all_exceed"
            (modified_observation, modified_observation_x0) = (observation_data,observation_x0)
        end
        #safe threshold for comparison
        threshold_vec[trial]=threshold
        #eceed sim
        #(modified_observation, modified_observation_x0,threshold) = exceed_cond_sim_quantile(20,num_sim,observation_data,observation_x0,threshold_quantile, alpha, coord_fine,coord_coarse,param,row_x0 )
        #(modified_observation, modified_observation_x0) = exceed_cond_sim(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )

        
        #fixed threshold, all exceed 
        #threshold=1.0
        #(modified_observation, modified_observation_x0) = (observation_data,observation_x0)
        

        #gridsize = Int(sqrt(size(coord_fine,1)))
        #coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha)
        l2=l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha,N_est_c) * size(modified_observation,1) #number of exceedances
        #l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l3=l_3_fun(modified_observation_x0, alpha, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha,0.0,1.0,50.0)
        log_likelihood_old=sum([l1,l2,l3,prior])
        #new param 

        #new params
        par_beta_old = param[2]
        par_c_old=param[1]
        alpha_old=alpha
        #propose new params
        eps_beta=0.05
        eps_c=0.05
        eps_alpha=0.05


        param[2]=uniform_proposal(par_beta_old,eps_beta,0.0,2.0)
        par_beta=param[2]
        param[1]=gaussian_proposal(par_c_old,eps_c)
        #param[1]=uniform_proposal(par_c_old,eps_c,0.0,10.0)
        par_c =param[1]
        alpha=gaussian_proposal(alpha_old,eps_alpha)
        #alpha=uniform_proposal(alpha_old,eps_alpha,0.0,10.0)


        #calculate new log likelihood
        l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha)
        #l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l2=l_2_fun_dependent_no_number_of_exceed(coord_fine, param,row_x0,alpha,N_est_c) * size(modified_observation,1) #number of exceedances
        l3=l_3_fun(modified_observation_x0, alpha, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha,0.0,1.0,50.0)
        log_likelihood_new =sum([l1,l2,l3,prior])  

        #MCMC update of param according to acceptance rate calculated with old and new likelihoods
        #param[2],log_likelihood_old=parameter_update(par_beta_old,par_beta,log_likelihood_old,log_likelihood_new)
        param[2],param[1],alpha,log_likelihood_old =parameter_update(par_beta_old,par_c_old,alpha_old,par_beta,par_c,alpha,log_likelihood_old,log_likelihood_new)
        #safe param after MCMC step
        par_beta_vec[trial+1]=param[2]
        par_c_vec[trial+1]=param[1]  
        par_alpha_vec[trial+1]=alpha 
        number_exceed_vec[trial]=size(modified_observation,1)
        log_likelihood_vec[trial]=log_likelihood_old
        #just print every n_trial_print's trial number to see how fast programme runs
        if trial%n_trial_print==0
            println("Trial: $trial")
            println("Number of exceedance: $(size(modified_observation,1)) " )
            println("Threshold: $threshold")
            println("alpha: $alpha")
            println("c: $(param[1])")
            println("beta: $(param[2])")
        end
    end
    Dict( "beta" => par_beta_vec, "c" => par_c_vec, "alpha" => par_alpha_vec, "threshold" => threshold_vec, "Number of exceedance"=> number_exceed_vec, "log_likelihood" => log_likelihood_vec)
end  

 


#MCMC Algorithm using approx risk functional
function MCMC_approx(N_MCMC,observation_data,observation_x0,threshold, alpha, coord_x0,coord_coarse,param,n_trial_print,N_est_c)
    num_obs=size(observation_data,1)
    coord_coarse_plus=vcat(coord_coarse,coord_x0')
    par_beta_vec  = repeat([param[2]],N_MCMC+1)
    par_c_vec = [param[1] for i=1:N_MCMC+1]
    par_alpha_vec = [alpha for i=1:N_MCMC+1]
    threshold_vec = [threshold for i=1:N_MCMC]
    number_exceed_vec = [-1 for i=1:N_MCMC]
    for trial in 1:N_MCMC
        #println("current trial: $trial and current param: $param")
        param[2]=par_beta_vec[trial]
        param[1]=par_c_vec[trial]
        alpha=par_alpha_vec[trial]
        
#fixed threshold, threshold exceedances are evaluated every step with changing parameters WITH THE APPROXIMATE RISK FUNCTIONAL

        (modified_observation, modified_observation_x0)=exceed_cond_sim_approx(num_obs,observation_data,observation_x0,threshold)  
        # exceed_cond_sim_approx(N_cond_sim,num_obs,observation_data,observation_x0,threshold, alpha, coord_fine,coord_coarse,param,row_x0 )
                    
        #safe threshold for comparison
        threshold_vec[trial]=threshold

        #gridsize = Int(sqrt(size(coord_fine,1)))
        #coord_cond_rows = get_common_rows_indices(coord_fine,floor.(coord_coarse.*gridsize)./gridsize)

        #l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha)
        l1=l_1_fun_approx(coord_coarse_plus,modified_observation,param, modified_observation_x0,alpha)
        #l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l2 = l_2_fun_approx(coord_coarse_plus,coord_x0,param,size(modified_observation,1),alpha,N_est_c)
        l3=l_3_fun(modified_observation_x0, alpha, threshold)  
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha,0.0,1.0,50.0)
        log_likelihood_old=sum([l1,l2,l3,prior])
        #new param 

        #new params
        par_beta_old=param[2]
        par_c_old=param[1]
        alpha_old=alpha
        #propose new params
        eps_beta=0.05
        eps_c=0.05
        eps_alpha=0.05


        param[2]=uniform_proposal(par_beta_old,eps_beta,0.0,2.0)
        par_beta=param[2]
        param[1]=gaussian_proposal(par_c_old,eps_c)
        #param[1]=uniform_proposal(par_c_old,eps_c,0.0,10.0)
        par_c =param[1]
        alpha=gaussian_proposal(alpha_old,eps_alpha)
        #alpha=uniform_proposal(alpha_old,eps_alpha,0.0,10.0)


        #calculate new log likelihood
        #l1=l_1_fun(coord_fine,coord_coarse,modified_observation,param, modified_observation_x0, row_x0,alpha)
        l1=l_1_fun_approx(coord_coarse_plus,modified_observation,param, modified_observation_x0,alpha)
        #l2= l_2_fun(coord_fine, param,row_x0, size(modified_observation,1),alpha,N_est_c)
        l2 = l_2_fun_approx(coord_coarse_plus,coord_x0,param,size(modified_observation,1),alpha,N_est_c)
        l3=l_3_fun(modified_observation_x0, alpha, threshold) 
        prior=log_likehood_log_gauss_1d_cut(param[1],0.0,1.5,500.0)+log_likehood_log_gauss_1d_cut(alpha,0.0,1.0,50.0)
        log_likelihood_new =sum([l1,l2,l3,prior])  

        #MCMC update of param according to acceptance rate calculated with old and new likelihoods
        #param[2],log_likelihood_old=parameter_update(par_beta_old,par_beta,log_likelihood_old,log_likelihood_new)
        param[2],param[1],alpha,log_likelihood_old =parameter_update(par_beta_old,par_c_old,alpha_old,par_beta,par_c,alpha,log_likelihood_old,log_likelihood_new)
        #safe param after MCMC step
        par_beta_vec[trial+1]=param[2]
        par_c_vec[trial+1]=param[1]  
        par_alpha_vec[trial+1]=alpha 
        number_exceed_vec[trial]=size(modified_observation,1)
        #just print every n_trial_print's trial number to see how fast programme runs
        if trial%n_trial_print==0
            println("Trial: $trial")
            println("Number of exceedance: $(size(modified_observation,1)) " )
            println("Threshold: $threshold")
            println("alpha: $alpha")
        end
    end
    Dict( "beta" => par_beta_vec, "c" => par_c_vec, "alpha" => par_alpha_vec, "threshold" => threshold_vec, "Number of exceedance"=> number_exceed_vec)
end  




#Look for the highest simulation number in the current directory and give the next number as string
function highest_sim_number(files::Vector{String})
    pattern = r"^simulated_observations_(\d+)\.jld2$"
    numbers = Int[]

    for f in files
        m = match(pattern, f)
        if m !== nothing
            push!(numbers, parse(Int, m.captures[1]))
        end
    end

    n= isempty(numbers) ? 1 : (maximum(numbers)+1)
    str = @sprintf("%04d", n)
end
function highest_sim_result_number(files::Vector{String})
    pattern = r"^sim_(\d+)\.jld2$"
    numbers = Int[]

    for f in files
        m = match(pattern, f)
        if m !== nothing
            push!(numbers, parse(Int, m.captures[1]))
        end
    end

    n= isempty(numbers) ? 1 : (maximum(numbers)+1)
    str = @sprintf("%04d", n)
end
#Look for the highest simulation number in the current directory and give the next number as string
function highest_folder_number(files::Vector{String})
    pattern = r"^sim_res_folder_(\d+)$"
    numbers = Int[]

    for f in files
        m = match(pattern, f)
        if m !== nothing
            push!(numbers, parse(Int, m.captures[1]))
        end
    end

    n= isempty(numbers) ? 1 : (maximum(numbers)+1)
    str = @sprintf("%02d", n)
end
