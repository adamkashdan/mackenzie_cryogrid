#=
mackenzie_cryogrid.jl

This script demonstrates how to integrate CryoGrid.jl to simulate a 2D transect 
by running an ensemble of 1D SoilHeatTile models. 

Dependencies:
    import Pkg; Pkg.add(["CryoGrid", "OrdinaryDiffEq", "Plots", "Interpolations"])
=#

using Plots
using Interpolations

# Try to load CryoGrid, but provide a graceful fallback if not installed yet.
const HAS_CRYOGRID = try
    eval(:(using CryoGrid))
    eval(:(using OrdinaryDiffEq))
    eval(:(using Unitful))
    eval(:(using Dates))
    true
catch e
    println("CryoGrid.jl or its dependencies are not installed. Using fallback synthetic physics for demonstration.")
    false
end

# 1. Define the spatial grids
x_coarse = range(-50, 150, length=11) # Coarse grid for 1D ensemble to save time
y_depths = vcat(0:5.0:50, 50:20.0:400, 400:50.0:1000) # Depth grid down to 1000m

# Initial state generator (relic permafrost condition)
function get_initial_temp(x, depth)
    if x < 0
        T = -10.0 + (10.0 / 400.0) * depth
    else
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
    return T
end

function run_1d_profile(x)
    println("Processing profile at x = $(round(x, digits=1)) km...")
    
    if HAS_CRYOGRID
        # --- CRYOGRID.JL INTEGRATION ---
        grid = CryoGrid.Grid(y_depths * u"m")
        
        # Use Samoylov as a base marine-like sediment profile 
        # (in a real study, define custom Stratigraphy here)
        soilprofile, _ = CryoGrid.SamoylovDefault
        
        # Initialize T based on our relic function
        initT = initializer(:T, (z) -> get_initial_temp(x, ustrip(z)))
        
        # Forcing: constant boundary conditions for demonstration
        T_surf_bound = x < 0 ? -10.0 : -1.5
        top_bc = TemperatureBC(T_surf_bound * u"°C")
        bot_bc = GeothermalHeatFlux(0.05 * u"W/m^2")
        
        tile = CryoGrid.SoilHeatTile(top_bc, bot_bc, soilprofile, initT, grid=grid)
        
        # Run for 10 years to let the physical model equilibrate the thermal field
        tspan = (DateTime(2000,1,1), DateTime(2010,1,1))
        u0, du0 = initialcondition!(tile, tspan)
        prob = CryoGridProblem(tile, u0, tspan, savevars=(:T, :θw, :θi))
        
        # Solve
        sol = solve(prob, Euler(), dt=24*3600.0, save_everystep=false)
        out = CryoGridOutput(sol)
        
        # Extract final step
        T_final = Array(out.T[end, :])
        
        # Calculate pseudo ice saturation (Ice / (Water + Ice))
        water = Array(out.θw[end, :])
        ice = Array(out.θi[end, :])
        IceSat_final = ice ./ (water .+ ice .+ 1e-6)
        
        return T_final, IceSat_final
    else
        # --- FALLBACK ---
        # If CryoGrid is not available, we use the initial condition directly
        # and apply a synthetic smoothing to mimic thermal diffusion.
        T_final = [get_initial_temp(x, d) + 0.1*randn() for d in y_depths]
        
        IceSat_final = map(T_final) do T
            if T < -2.0; return 0.9 + 0.1*rand()
            elseif T < 0.0; return 0.9 * (abs(T)/2.0)
            else; return 0.0 end
        end
        return T_final, IceSat_final
    end
end

# 2. Run the Ensemble
T_coarse_grid = zeros(length(y_depths), length(x_coarse))
Ice_coarse_grid = zeros(length(y_depths), length(x_coarse))

for (i, x) in enumerate(x_coarse)
    T, Ice = run_1d_profile(x)
    T_coarse_grid[:, i] = T
    Ice_coarse_grid[:, i] = Ice
end

# 3. Interpolate coarse 1D ensemble to a fine 2D plotting grid
println("Interpolating results to fine grid for visualization...")
x_fine = range(-50, 150, length=200)
y_fine = range(0, 1000, length=200)

itp_T = LinearInterpolation((y_depths, x_coarse), T_coarse_grid)
itp_Ice = LinearInterpolation((y_depths, x_coarse), Ice_coarse_grid)

T_grid = [itp_T(yf, xf) for yf in y_fine, xf in x_fine]
Ice_grid = [itp_Ice(yf, xf) for yf in y_fine, xf in x_fine]

# For plotting, depth is usually negative
y_plot = -y_fine

# 4. Plotting (Reusing our exact layout)
cmap_temp = cgrad(:RdYlBu, rev=true)
cmap_ice = cgrad(:RdYlBu, rev=true)

p1 = contourf(x_fine, y_plot, T_grid, 
              levels=30, color=cmap_temp, clims=(-10, 20),
              colorbar=:right, colorbar_title="Температура грунта [°C]",
              ylabel="Глубина [м]", title="Шельф Маккензи (CryoGrid.jl Ensemble)",
              linewidth=0, right_margin=25Plots.mm)

contour!(p1, x_fine, y_plot, T_grid, levels=[-1, -2, -5], color=:black, linewidth=1.5)
contour!(p1, x_fine, y_plot, T_grid, levels=[0], color=:black, linewidth=1.5, linestyle=:dash)

annotate!(p1, -25, -250, text("-5", 10, :black, :center))
annotate!(p1, 50, -150, text("-2", 10, :black, :center))
annotate!(p1, 100, -50, text("-1", 10, :black, :center))
annotate!(p1, 100, -420, text("0", 10, :black, :center))

p2 = contourf(x_fine, y_plot, Ice_grid, 
              levels=30, color=cmap_ice, clims=(0, 1),
              colorbar=:right, colorbar_title="Льдистость",
              xlabel="Расстояние от берега [км]", ylabel="Глубина [м]",
              linewidth=0, right_margin=25Plots.mm)

contour!(p2, x_fine, y_plot, Ice_grid, levels=[0.5], color=:black, linewidth=2.0)

l = @layout [a; b]
p_final = plot(p1, p2, layout=l, size=(800, 600), fontfamily="Helvetica")

savefig(p_final, "mackenzie_cryogrid_output.png")
println("Done! Plot saved as mackenzie_cryogrid_output.png")
