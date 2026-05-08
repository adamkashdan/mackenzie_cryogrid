using Plots

# Define the spatial grid
x = range(-50, 150, length=200) # Distance to Shoreline [km]
y = range(-1000, 0, length=200) # Depth [m]

# Synthetic Temperature Field Generator
function temp_field(x, y)
    depth = -y
    
    if x < 0
        # Onshore: surface -10, 0 degree at 400m
        T = -10.0 + (10.0 / 400.0) * depth
    else
        # Offshore: Abrupt step at shoreline
        T_surf = -1.5
        core_depth = 150.0
        T_core = -3.5 * exp(-x/60.0) - 1.0
        
        if depth <= core_depth
            T = T_surf + (T_core - T_surf) * (depth / core_depth)
        elseif depth <= 400.0
            T = T_core + (0.0 - T_core) * ((depth - core_depth) / (400.0 - core_depth))
        else
            T = 0.0 + 0.03 * (depth - 400.0)
        end
    end
    
    # Add a bit of realistic noise
    T += 0.1 * randn()
    
    return T
end

# Synthetic Ice Saturation Generator
function ice_sat(T)
    if T < -2.0
        return 0.9 + 0.1 * rand()
    elseif T < 0.0
        return 0.9 * (abs(T) / 2.0)
    else
        return 0.0
    end
end

# Generate Data Grids
T_grid = [temp_field(xi, yi) for yi in y, xi in x]
Ice_grid = [ice_sat(T_grid[j, i]) for j in 1:length(y), i in 1:length(x)]

# Define Colormaps
cmap_temp = cgrad(:RdYlBu, rev=true)
cmap_ice = cgrad(:RdYlBu, rev=true)

# Plot Temperature Panel
p1 = contourf(x, y, T_grid, 
              levels=30, 
              color=cmap_temp, 
              clims=(-10, 20),
              colorbar=:right,
              colorbar_title="Температура грунта [°C]",
              ylabel="Глубина [м]",
              title="Шельф Маккензи (Mackenzie)",
              linewidth=0,
              right_margin=25Plots.mm)

# Overlay solid contour lines (-1, -2, -5)
contour!(p1, x, y, T_grid, levels=[-1, -2, -5], color=:black, linewidth=1.5)

# Overlay dashed 0 contour line
contour!(p1, x, y, T_grid, levels=[0], color=:black, linewidth=1.5, linestyle=:dash)

# Add annotations for the contour lines
annotate!(p1, -25, -250, text("-5", 10, :black, :center))
annotate!(p1, 50, -150, text("-2", 10, :black, :center))
annotate!(p1, 100, -50, text("-1", 10, :black, :center))
annotate!(p1, 100, -420, text("0", 10, :black, :center))

# Plot Ice Saturation Panel
p2 = contourf(x, y, Ice_grid, 
              levels=30, 
              color=cmap_ice, 
              clims=(0, 1),
              colorbar=:right,
              colorbar_title="Льдистость",
              xlabel="Расстояние от берега [км]",
              ylabel="Глубина [м]",
              linewidth=0,
              right_margin=25Plots.mm)

# Overlay 0.5 contour line for Ice Saturation
contour!(p2, x, y, Ice_grid, levels=[0.5], color=:black, linewidth=2.0)

# Combine Layout
l = @layout [a; b]
p_final = plot(p1, p2, layout=l, size=(800, 600), fontfamily="Helvetica")

# Save the plot
savefig(p_final, "mackenzie_reproduction.png")
println("Plot saved as mackenzie_reproduction.png")
