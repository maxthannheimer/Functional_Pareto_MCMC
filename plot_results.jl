include("functions.jl")

#import Pkg;Pkg.add("JLD2")
using JLD2

N_MCMC=2000
gridsize=20
repetition=1
result=load("$(N_MCMC)_simulation_$(gridsize)_grid_0$repetition.jld2")["single_stored_object"]

N_burn_in=0


#mean of params with burn in period canceled out
[mean(result["beta"][N_burn_in+1:N_MCMC+1]),mean(result["c"][N_burn_in+1:N_MCMC+1]),mean(result["alpha"][N_burn_in+1:N_MCMC+1])]

#mean of params with burn in period canceled out WITH EXP
#[mean(result["beta"][N_burn_in+1:N_MCMC+1]),mean(exp.(result["c"][N_burn_in+1:N_MCMC+1])),mean(result["alpha"][N_burn_in+1:N_MCMC+1])]



#plot and hist for beta,zeta,lambda_1,lambda_2
histogram(result["beta"][N_burn_in+1:N_MCMC+1],title="Histogram for β, Burn in: $N_burn_in")
p_beta=scatter(1:N_MCMC+1,result["beta"],label="$N_MCMC samples of beta",title="Markov chain for β")
hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β: $(round((mean(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)

# with exp(c) 
#histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
#p_c=scatter(1:N_MCMC+1,exp.(result["c"]),label="$N_MCMC samples of c",title="Markov chain for c, Burn in: $N_burn_in")
#hline!(p_c, [mean(exp.(result["c"][N_burn_in+1:N_MCMC+1]))],label="Mean estimate for c: $(round((mean(exp.(result["c"][N_burn_in+1:N_MCMC+1]))),digits=3))",linewidth=3)

histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
p_c=scatter(1:N_MCMC+1,result["c"],label="$N_MCMC samples of c",title="Markov chain for c, Burn in: $N_burn_in")
hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c: $(round((mean(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)

histogram(result["alpha"][N_burn_in+1:N_MCMC+1],title="Histogram for α, Burn in: $N_burn_in")
p_alpha=scatter(1:N_MCMC+1,result["alpha"],label="$N_MCMC samples of α",title="Markov chain for α")
hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α: $(round((mean(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)


