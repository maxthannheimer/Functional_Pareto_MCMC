include("functions.jl")

#import Pkg;Pkg.add("JLD2")
using JLD2

alpha_true = 2.0
beta_true=1.5
c_true=2.0

N_MCMC=8000
gridsize=5
#repetition=1 #very bad alpha and c in non quantile threshold case
#repetition=3
#(result1,result2)=load("sim_0$(repetition).jld2")["single_stored_object"]
N_burn_in=0




#result=result1
#result=result2

#mean of params with burn in period canceled out
#[mean(result["beta"][N_burn_in+1:N_MCMC+1]),mean(result["c"][N_burn_in+1:N_MCMC+1]),mean(result["alpha"][N_burn_in+1:N_MCMC+1])]

plots = Vector{}(undef, 6)


for repetition in 1:5
(result1,result2)=load("sim_0$(repetition).jld2")["single_stored_object"]
for res in [(result1,"quantile",1),(result2,"fixed",4)]
    (result,method,i)=res


#histogram(result["beta"][N_burn_in+1:N_MCMC+1],title="Histogram for β, Burn in: $N_burn_in")
p_beta=scatter(1:N_MCMC+1,result["beta"],label="$N_MCMC samples of beta",title="Markov chain for β($method), Burn in: $N_burn_in")
hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β: $(round((mean(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_beta, [median(result["beta"][N_burn_in+1:N_MCMC+1])],label="Median estimate for β: $(round((median(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_beta, [beta_true for i in N_burn_in+1:N_MCMC+1], label="True value for β: $beta_true ",linewidth=3)

#histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
p_c=scatter(1:N_MCMC+1,result["c"],label="$N_MCMC samples of c",title="Markov chain for c ($method), Burn in: $N_burn_in")
hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c: $(round((mean(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_c, [median(result["c"][N_burn_in+1:N_MCMC+1])],label="Median estimate for c: $(round((median(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_c, [c_true for i in N_burn_in+1:N_MCMC+1], label="True value for c: $c_true ",linewidth=3)

#histogram(result["alpha"][N_burn_in+1:N_MCMC+1],title="Histogram for α, Burn in: $N_burn_in")
p_alpha=scatter(1:N_MCMC+1,result["alpha"],label="$N_MCMC samples of α",title="Markov chain for α ($method), Burn in: $N_burn_in")
hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α: $(round((mean(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_alpha, [median(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Median estimate for α: $(round((median(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)
hline!(p_alpha, [alpha_true for i in N_burn_in+1:N_MCMC+1], label="True value for α: $alpha_true ",linewidth=3)

plots[i+0],plots[i+1],plots[i+2]=p_beta,p_alpha,p_c

end
combined = plot(plots[1],plots[4],plots[2],plots[5],plots[3],plots[6], layout = (3,2),size=(1000,1500))
savefig(combined, "six_plots_$repetition.pdf")
end







include("functions.jl")

#import Pkg;Pkg.add("JLD2")
using JLD2

N_MCMC=3000
gridsize=10
repetition=0
result=load("$(N_MCMC)_simulation_$(gridsize)_grid_0$repetition.jld2")["single_stored_object"]
N_burn_in=0


#mean of params with burn in period canceled out
[mean(result["beta"][N_burn_in+1:N_MCMC+1]),mean(result["c"][N_burn_in+1:N_MCMC+1]),mean(result["alpha"][N_burn_in+1:N_MCMC+1])]



#plot and hist for beta,zeta,lambda_1,lambda_2
histogram(result["beta"][N_burn_in+1:N_MCMC+1],title="Histogram for β, Burn in: $N_burn_in")
p_beta=scatter(1:N_MCMC+1,result["beta"],label="$N_MCMC samples of beta",title="Markov chain for β")
hline!(p_beta, [mean(result["beta"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for β: $(round((mean(result["beta"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)

histogram(result["c"][N_burn_in+1:N_MCMC+1],title="Histogram for c")
p_c=scatter(1:N_MCMC+1,result["c"],label="$N_MCMC samples of c",title="Markov chain for c, Burn in: $N_burn_in")
hline!(p_c, [mean(result["c"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for c: $(round((mean(result["c"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)

histogram(result["alpha"][N_burn_in+1:N_MCMC+1],title="Histogram for α, Burn in: $N_burn_in")
p_alpha=scatter(1:N_MCMC+1,result["alpha"],label="$N_MCMC samples of α",title="Markov chain for α")
hline!(p_alpha, [mean(result["alpha"][N_burn_in+1:N_MCMC+1])],label="Mean estimate for α: $(round((mean(result["alpha"][N_burn_in+1:N_MCMC+1])),digits=3))",linewidth=3)


