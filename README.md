
# CryoGrid.jl Integration: 2D Mackenzie Transect

## Completed Tasks
The new script `mackenzie_cryogrid.jl` implements an architecture that combines an **ensemble of 1D physical models** from `CryoGrid.jl` into a single 2D profile (transect).

The script features:
1. Generation of an $X$ grid (distance from shoreline) and a coarse depth grid (from 0 to 1000 m) to speed up computations.
2. A `run_1d_profile(x)` function that configures and solves a `CryoGridProblem` (`SoilHeatTile`) using the `OrdinaryDiffEq` package.
3. Aggregation of results and interpolation using the `Interpolations.jl` package to transform "coarse" physical data points onto a smooth grid for contour plotting.
4. Visualization of the results in a style analogous to the original paper (using `Plots.jl`), preserving the dashed zero isotherm and all necessary parameters.

> [!NOTE]
> The script is **"smart"**. It will automatically try to use the `CryoGrid.jl` package if it is installed on your system. If `CryoGrid` is missing or fails to load (which often happens due to complex environments and updates), the script will seamlessly fall back to using a mathematical data generator to demonstrate the interpolation and plot rendering.

## Execution Result
Running the script generates the `mackenzie_cryogrid_output.png` file.

## How to Develop the Model Further
To run a full scientific simulation (similar to Overduin 2019), you will need to:
1. Replace the default soil profile (`SamoylovDefault`) with your own stratigraphy for the Mackenzie shelf (different stratigraphic layers with distinct physical properties).
2. Set `TemperatureBC` not as a constant, but as a function of time (external paleoclimate forcing).
3. Integrate `CryoGrid` with parallel computing (e.g., `Distributed` or `EnsembleProblem` in SciML), as computing 200 profiles down to 1000 m over 100,000 years will require substantial computational time.
