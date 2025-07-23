include("functions.jl")

param=[1.0,1.0]
gridsize=500
#get fine coordinates as vector, just ignore coarse ones
coord_fine=Create_Grid_and_Observation(gridsize,5)[2]

#fbm as matrix and heatmap as plot


#scalefontsizes(1.5) 
plots = Vector{}(undef, 4)


gaussian_sample_tmp=FBM_simu_fast(param,gridsize,1)[1]
heatmap(coord_fine[1:gridsize,1],coord_fine[1:gridsize,1],gaussian_sample_tmp)

for i in 1:4
plots[i]=heatmap(coord_fine[1:gridsize,1],coord_fine[1:gridsize,1],gaussian_sample_tmp[i])
end
plot(plots[1],plots[3],plots[2],plots[4], layout = (2,2),size=(1000,900))

gaussian_sample[4]=gaussian_sample_tmp

for (risk,str) in [(mean,"mean(X)="),(maximum,"sup(X)="),(x->x[50,50],"X(0.5,0.5)=")] #,(minimum,"inf()=")
#create empty plot list
plots = Vector{}(undef, 4)
#plot each realisation
for i in 1:4
plots[i]=heatmap(coord_fine[1:gridsize,1],coord_fine[1:gridsize,1],gaussian_sample[i]; title="$str$(round(risk(gaussian_sample[i]),digits=3))")
end


# Overlay red rectangle
plot!(x_frame, y_frame, linecolor=:red, linewidth=3, label="", line=:solid)

combined = plot(plots[1],plots[3],plots[2],plots[4], layout = (2,2),size=(1000,900))
savefig(combined, "$str.png")
end
combined = plot(plots[1],plots[3],plots[2],plots[4], layout = (2,2),size=(1000,1000),framestyle = :box)
#combined = plot(plots[1],plots[2],plots[3],layout = (3,1),size=(1500,500))

mean(gaussian_sample)
maximum(gaussian_sample)
gaussian_sample[50,50]


using Plots

# Sample plot
plot(rand(10), rand(10), label="Data", legend=false)

# Get plot limits
xli = xlims()
yli = ylims()

# Coordinates of the frame (corners of the rectangle)
x_frame = [xli[1], xli[2], xli[2], xli[1], xli[1]]
y_frame = [yli[1], yli[1], yli[2], yli[2], yli[1]]

# Overlay red rectangle
plot!(x_frame, y_frame, linecolor=:red, linewidth=3, label="", line=:solid, framestyle = :box )



#fbm as vector and surface plot
gaussian_sample_vec=FBM_simu_fast_vec(param,gridsize,1)[1]
surface(coord_fine[:,1],coord_fine[:,2],gaussian_sample_vec)


