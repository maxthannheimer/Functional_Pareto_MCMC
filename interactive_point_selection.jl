using GLMakie, CSV, DataFrames, Plots

# Prepare the figure and axis
f = Figure(resolution = (600, 600))
ax = Axis(f[1, 1];
    title = "Click to select 10 points",
    limits = (0, 11, 0, 11),
    aspect = DataAspect())

# Plot the 10×10 grid in light gray
xs = repeat(1:10, inner=10)
ys = repeat(1:10, outer=10)
scatter!(ax, xs, ys, color = :lightgray, markersize = 8)

selected = Tuple{Int, Int}[]

# Callback on mouse‐click
function mouse_click_callback(ev)
    if ev.button == Mouse.left && ev.action == Mouse.press
        # Get the mouse position using the mouseposition function
        p = mouseposition(ax)
        # Snap to the nearest integer grid point
        x = clamp(round(Int, p[1]), 1, 10)
        y = clamp(round(Int, p[2]), 1, 10)
        push!(selected, (x, y))
        # Mark the selected point in red
        scatter!(ax, [x], [y], color = :red, markersize = 12)
        
        # Refresh the figure by displaying it
        display(f)  # display function automatically handles the render
        
        # Once we have 10 points, save & notify
        if length(selected) == 10
            df = DataFrame(x = Int[], y = Int[])
            for pt in selected
                push!(df, (pt[1], pt[2]))
            end
            CSV.write("selected_points.csv", df)
            println("✅ 10 points selected and saved to selected_points.csv")
        end
    end
end

# Attach the callback function to mouse click event
on(ax.scene.events.mousebutton) do ev
    mouse_click_callback(ev)
end

display(f)
