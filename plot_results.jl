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

(res_200,res_100,result_ap)=load("sim_21.jld2")["single_stored_object"]
res_200_beta=res_200["beta"]
res_200_alpha=res_200["alpha"]
res_200_c=res_200["c"]
res_100_beta=res_100["beta"]
res_100_alpha=res_100["alpha"]
res_100_c=res_100["c"]
result_ap_beta=result_ap["beta"]
result_ap_alpha=result_ap["alpha"]
result_ap_c=result_ap["c"]

for repetition in 2:5
    (res200,res100,resap)=load("sim_2$(repetition).jld2")["single_stored_object"]
res_200_beta=vcat(res_200_beta,res200["beta"])
res_200_alpha=vcat(res_200_alpha,res200["alpha"])
res_200_c=vcat(res_200_c,res200["c"])
res_100_beta=vcat(res_100_beta,res100["beta"])
res_100_alpha=vcat(res_100_alpha,res100["alpha"])
res_100_c=vcat(res_100_c,res100["c"])
result_ap_beta=vcat(result_ap_beta,resap["beta"])
result_ap_alpha=vcat(result_ap_alpha,resap["alpha"])
result_ap_c=vcat(result_ap_c,resap["c"])
end

plots = Vector{}(undef, 9)

plots[1]=histogram(res_200_beta,title="Histogram for β, 200 N_cond_sim")
vline!(plots[1], [beta_true],label=" β: $(beta_true)",linewidth=3)
vline!([mean(res_200_beta)],label="Mean: $(round(mean(res_200_beta),digits=3))",linecolor=:brown,width=3)

plots[2]=histogram(res_100_beta,title="Histogram for β, 100 N_cond_sim")
vline!(plots[2], [beta_true],label=" β: $(beta_true)",linewidth=3)
vline!([mean(res_100_beta)],linecolor=:brown,label="Mean: $(round(mean(res_100_beta),digits=3))",linewidth=3)

plots[3]=histogram(result_ap_beta,title="Histogram for β, approx_risk")
vline!(plots[3], [beta_true],label=" β: $(beta_true)",linewidth=3)
vline!([mean(result_ap_beta)],linecolor=:brown,label="Mean: $(round(mean(result_ap_beta[N_burn_in+1:N_MCMC+1]),digits=3))",linewidth=3)

plots[4]=histogram(res_200_alpha,title="Histogram for α, 200 N_cond_sim")
vline!(plots[4], [alpha_true],label=" α: $(alpha_true)",linewidth=3)
vline!([mean(res_200_alpha)],linecolor=:brown,label="Mean: $(round(mean(res_200_alpha),digits=3))",linewidth=3)

plots[5]=histogram(res_100_alpha,title="Histogram for α, 100 N_cond_sim")
vline!(plots[5], [alpha_true],label=" α: $(alpha_true)",linewidth=3)
vline!([mean(res_100_alpha)],linecolor=:brown,label="Mean: $(round(mean(res_100_alpha),digits=3))",linewidth=3)

plots[6]=histogram(result_ap_alpha,title="Histogram for α, approx_risk")
vline!(plots[6], [alpha_true],label=" α: $(alpha_true)",linewidth=3)
vline!([mean(result_ap_alpha)],linecolor=:brown,label="Mean: $(round(mean(result_ap_alpha[N_burn_in+1:N_MCMC+1]),digits=3))",linewidth=3)

plots[7]=histogram(res_200_c,title="Histogram for c, 200 N_cond_sim")
vline!(plots[7], [c_true],label=" c: $(c_true)",linewidth=3)
vline!([mean(res_200_c)],linecolor=:brown,label="Mean: $(round(mean(res_200_c),digits=3))",linewidth=3)

plots[8]=histogram(res_100_c,title="Histogram for c, 100 N_cond_sim")
vline!(plots[8], [c_true],label=" c: $(c_true)",linewidth=3)
vline!([mean(res_100_c)],linecolor=:brown,label="Mean: $(round(mean(res_100_c),digits=3))",linewidth=3)

plots[9]=histogram(result_ap_c,title="Histogram for c, approx_risk")
vline!(plots[9], [c_true],label=" c: $(c_true)",linewidth=3)
vline!([mean(result_ap_c)],linecolor=:brown,label="Mean: $(round(mean(result_ap_c[N_burn_in+1:N_MCMC+1]),digits=3))",linewidth=3)

combined = plot(plots[1],plots[4],plots[7],plots[2],plots[5],plots[8],plots[3],plots[6],plots[9], layout = (3,3),size=(1500,1500))
savefig(combined, "nine_plots_Histogram.pdf")





#result=result1
#result=result2

#mean of params with burn in period canceled out
#[mean(result["beta"][N_burn_in+1:N_MCMC+1]),mean(result["c"][N_burn_in+1:N_MCMC+1]),mean(result["alpha"][N_burn_in+1:N_MCMC+1])]


############################################################################
#plot with three simulations



plots = Vector{}(undef, 9)



for repetition in 1:5
(result1,result2,result_approx)=load("sim_2$(repetition).jld2")["single_stored_object"]
for res in [(result1,"200 N_cond_sim",1),(result2,"100 N_cond_sim",4),(result_approx,"approx_risk",7)]
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
combined = plot(plots[1],plots[4],plots[7],plots[2],plots[5],plots[8],plots[3],plots[6],plots[9], layout = (3,3),size=(1500,1500))
savefig(combined, "nine_plots_$repetition.pdf")
end




############################################################################
#plot with three simulations

#=

plots = Vector{}(undef, 9)


for repetition in 1:5
(result1,result2,result_approx)=load("sim_2$(repetition).jld2")["single_stored_object"]
for res in [(result1,"200 N_cond_sim",1),(result2,"100 N_cond_sim",4),(result_approx,"approx_risk",7)]
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
combined = plot(plots[1],plots[4],plots[7],plots[2],plots[5],plots[8],plots[3],plots[6],plots[9], layout = (3,3),size=(1500,1500))
savefig(combined, "nine_plots_$repetition.pdf")
end

=#

############################################################################
#plot with two simulations

#=

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


=#