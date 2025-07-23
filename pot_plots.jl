using Plots
using Random, Distributions

n=100
vals=exp.(rand(Normal(0, 1), n) )
q=quantile(vals,0.95)
scalefontsizes(1*1.3) 
scatter(vals, title="Time series data",label="" ,
    marker=:circle,
    markercolor=:white,           # non-filled
    markerstrokecolor=:black,     # circle outline color
    markerstrokewidth=0.5,        # thin outline
    markersize=4,
    legend=:topright)

vals_above_q = findall(vals .> q)
scatter!(vals_above_q ,vals[vals_above_q], label="Threshold exceedances",
    marker=:circle,
    markercolor=:red,           # non-filled
    markerstrokecolor=:red,     # circle outline color
    markerstrokewidth=0.5,        # thin outline
    markersize=4,
    legend=:bottomright)
hline!([q], color=:red, linestyle=:dash, label="0.95 quantile", linewidth=1.5,legend=:topright)


