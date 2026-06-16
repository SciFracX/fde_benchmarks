"""
    Fractional-Order ODE (FODE) Solver Benchmarking Framework
    
    This script performs a comprehensive benchmarking comparison of fractional-order
    differential equation solvers across three major computational platforms:
    
    Platforms Compared:
    1. Julia: FractionalDiffEq.jl library (PECE, PIEX, PITrap methods)
    2. MATLAB: Native FODE solvers (PIEX, PECE, PIRect, PITrap, FOTF/NLFODE)
    3. Python: PyCaputo library (PECE method)
    4. Legacy Julia: FdeSolver.jl (PECE method)
    
    Key Metrics:
    - Execution time (seconds) - How fast each solver completes
    - Accuracy (infinity norm error) - How close solution is to reference
    - Error-time tradeoff - Efficiency frontier of each implementation
    
    Test Problem:
    Noncommensurate FODE system with mixed fractional derivative orders,
    representing realistic physical/biological systems that don't have
    uniform memory dynamics across all variables.
    
    Visualization:
    Two identical plots showing error vs. time relationship for each solver,
    enabling visual identification of:
    - Fast but inaccurate solvers (right side, low error-time product)
    - Slow but accurate solvers (left side, better accuracy)
    - Optimal solvers (steep slope = high accuracy improvement per time unit)
"""

# ========================================================================
# Import Required Libraries
# ========================================================================

using DataFrames      # Tabular data structures for handling benchmark results
using Statistics      # Statistical functions (median, etc.)
using CSV            # CSV file I/O for benchmark data persistence
using CairoMakie     # High-quality vector graphics for scientific plots

# ========================================================================
# Data Structure for Benchmark Results
# ========================================================================

struct wps
    """
    Work-Precision Set (WPS) - stores benchmark results for one solver.
    
    Fields:
        name::String - Solver name and method (e.g., "FractionalDiffEq.jl PECE")
        times::Vector{Float64} - CPU execution time for each accuracy level (seconds)
        errors::Vector{Float64} - Achieved error (L∞ norm) at each time point
    
    A work-precision set captures the fundamental speed-accuracy tradeoff:
    as solvers use higher accuracy (smaller step size), they take more time
    but achieve better error metrics. The set enables direct visualization
    of this relationship on a log-log plot.
    """
    name        # Descriptive name of the solver implementation
    times       # CPU wall-clock time measurements in seconds
    errors      # Accuracy measurements (infinity norm error ||f(u*)||∞)
end

# Initialize collection for all benchmark results
wps_set = Any[]  # Array to accumulate work-precision sets from all solvers

# ========================================================================
# Solver Registry - Metadata for All Tested Implementations
# ========================================================================
# This registry defines which solvers are included in the benchmark comparison.
# Format: (pkg=:PackageName, name=\"DisplayName\")
# The pkg field determines visual styling (line style, colors) in the plot.
solvers_all = [   
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PIEX", )
    #(; pkg = :MATLAB,                           name = "MATLAB PIEX", )
];

# ========================================================================
# SECTION 1: Load Julia Benchmark Data
# ========================================================================
# Each CSV file contains work-precision data: Column 1 = Time (s), Column 2 = Error

# Legacy Julia FODE Solver - PECE method
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FdeSolver_PECE.csv"))
push!(wps_set, wps("FdeSolver.jl PECE", df[:,1], df[:,2]))

# Modern Julia FractionalDiffEq.jl - PECE method
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PECE.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PECE", df[:,1], df[:,2]))

# FractionalDiffEq.jl - PITrap method (Trapezoidal quadrature rule)
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PITrap.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

# ========================================================================
# SECTION 2: Load MATLAB Benchmark Data
# ========================================================================
# MATLAB fractional calculus toolbox provides multiple numerical methods

# MATLAB PIEX - Predictor-Integrator Extrapolation
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PIEX.csv"))
push!(wps_set, wps("MATLAB PIEX", df[:,1], df[:,2]))

# MATLAB PECE - Predictor-Evaluator-Corrector-Evaluator
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PECE.csv"))
push!(wps_set, wps("MATLAB PECE", df[:,1], df[:,2]))

# MATLAB PIRect - Predictor-Integrator Rectangular rule
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PIRect.csv"))
push!(wps_set, wps("MATLAB PIRect", df[:,1], df[:,2]))

# MATLAB PITrap - Predictor-Integrator Trapezoidal rule (standard choice)
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PITrap.csv"))
push!(wps_set, wps("MATLAB PITrap", df[:,1], df[:,2]))

# MATLAB FOTF - Fractional-Order Transfer Function (vectorized nonlinear FODE)
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_NLFODE_VEC.csv"))
push!(wps_set, wps("MATLAB FOTF", df[:,1], df[:,2]))

# ========================================================================
# SECTION 3: Load Python Benchmark Data
# ========================================================================
# PyCaputo - Primary Python library for fractional calculus
# Note: Column indexing offset due to CSV format differences

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/PYCAPUTO_PECE.csv"))
push!(wps_set, wps("PyCaputo PECE", df[:,2], df[:,3]))

# ========================================================================
# SECTION 4: Visualization Configuration and Plotting
# ========================================================================
# Create a work-precision diagram comparing all FODE solvers
# This is the standard publication format for benchmarking numerical methods

fig = begin
    # ---- Define Visual Styling ----
    # Map each package to a distinct line style for clear differentiation
    LINESTYLES = Dict(
        :FdeSolvers => :dash,      # Legacy Julia: dashed line
        :FractionalDiffEq => :solid,   # Modern Julia: solid line (primary focus)
        :MATLAB => :dot,           # MATLAB: dotted line
        :Python => :dashdot        # Python: dash-dot line
    )
    
    # ---- Figure Dimensions ----
    # Professional 16:9 aspect ratio with sufficient resolution
    ASPECT_RATIO = 0.7   # Height = 70% of width
    WIDTH = 1200         # Pixels - sufficient for publication
    HEIGHT = round(Int, WIDTH * ASPECT_RATIO)
    STROKEWIDTH = 2.5    # Line thickness for axes and legends

    # ---- Color Scheme ----
    # Use seaborn_bright colormap for visually distinct solver curves
    colors = cgrad(:seaborn_bright, length(solvers_all); categorical = true)
    
    # ---- Plot Theme ----
    cycle = Cycle([:marker], covary = true)
    plot_theme = Theme(Lines = (; cycle), Scatter = (; cycle))

    # ---- Create Figure with Theme ----
    with_theme(plot_theme) do 
        fig = Figure(; size = (WIDTH, HEIGHT))
        ax = Axis(fig[1, 1], ylabel = L"Time $\mathbf{(s)}$",
        # ========================================================================
        # SECTION 4A: First Work-Precision Diagram (Top)
        # ========================================================================
            xlabelsize = 22, ylabelsize = 22,
            xlabel = L"Error: $\mathbf{||f(u^\ast)||_\infty}$",
            xscale = log10, yscale = log10, xtickwidth = STROKEWIDTH,
            ytickwidth = STROKEWIDTH, spinewidth = STROKEWIDTH,
            xticklabelsize = 20, yticklabelsize = 20)

        idxs = sortperm(median.(getfield.(wps_set, :times)))

        ls, scs = [], []

        for (i, (wp, solver)) in enumerate(zip(wps_set[idxs], solvers_all[idxs]))
            (; name, times, errors) = wp
            #errors = [err.l∞ for err in errors]
            l = lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors[i])
            sc = scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors[i])
            push!(ls, l)
            push!(scs, sc)
        end

        xlims!(ax; high=1e1)
        ylims!(ax; low=10^(-3.7), high=10^(-1.5))

        Legend(fig[1,2], [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[0, :] = Label(fig, "Noncommensurate FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)

        ############ bottom plot ############
        # First work-precision plot: Error vs. time relationship
        # X-axis: Error magnitude (log scale, smaller is better)
        # Y-axis: CPU time (log scale, smaller is better)
        ax = Axis(
            fig[1, 1],
            ylabel = L"Time $\mathbf{(s)}$",
            xlabelsize = 22,
            ylabelsize = 22,
            xlabel = L"Error: $\mathbf{||f(u^\ast)||_\infty}$",
            xscale = log10, yscale = log10, xtickwidth = STROKEWIDTH,
            ytickwidth = STROKEWIDTH, spinewidth = STROKEWIDTH,
            xticklabelsize = 20, yticklabelsize = 20)

        idxs = sortperm(median.(getfield.(wps_set, :times)))

        ls, scs = [], []

        for (i, (wp, solver)) in enumerate(zip(wps_set[idxs], solvers_all[idxs]))
            (; name, times, errors) = wp
            #errors = [err.l∞ for err in errors]
            l = lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors[i])
            sc = scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors[i])
            push!(ls, l)
            push!(scs, sc)
        end

        xlims!(ax; high=1e1)
        ylims!(ax; low=10^(-3.7), high=10^(-1.5))

        Legend(fig[3,2], [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[2, :] = Label(fig, "Noncommensurate FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)
        fig
    end
end

save("fode_benchmarks.svg", fig)