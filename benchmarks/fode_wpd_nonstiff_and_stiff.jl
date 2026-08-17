using DataFrames
using Statistics
using CSV
using CairoMakie

# Lightweight container for one solver curve (name + time/error samples).
struct wps
    name
    times
    errors
end

# Non-stiff benchmark dataset collection.
wps_set1 = Any[]


# Solver metadata for the non-stiff panel.
# `pkg` controls line style; `name` is shown in the legend.
solvers_all_nonstiff = [
    (; pkg = :FdeSolvers,                       name = "FdeSolver.jl PECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PITrap", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PIRect", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl KernelCompression", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl SoE", )    
    (; pkg = :MATLAB,                           name = "MATLAB PIEX", )
    (; pkg = :MATLAB,                           name = "MATLAB PECE", )
    (; pkg = :MATLAB,                           name = "MATLAB PIRect", )
    (; pkg = :MATLAB,                           name = "MATLAB PITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB FOTF", )
    (; pkg = :Python,                           name = "pycaputo PECE", )
];

##### Julia #####
# Load Julia non-stiff benchmark outputs (columns: time, error).
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FdeSolver_PECE.csv"))
push!(wps_set1, wps("FdeSolver.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PECE.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PITrap.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PIRect.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl PIRect", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_KernelCompression.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl KernelCompression", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_SoE.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl SoE", df[:,1], df[:,2]))

##### MATLAB #####
# Load MATLAB non-stiff benchmark outputs.
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PIEX.csv"))
push!(wps_set1, wps("MATLAB PIEX", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PECE.csv"))
push!(wps_set1, wps("MATLAB PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PIRect.csv"))
push!(wps_set1, wps("MATLAB PIRect", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_PITrap.csv"))
push!(wps_set1, wps("MATLAB PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_NLFODE_VEC.csv"))
push!(wps_set1, wps("MATLAB FOTF", df[:,1], df[:,2]))

##### Python #####
# Python CSV includes index column first, so use columns 2 and 3.
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/PYCAPUTO_PECE.csv"))
push!(wps_set1, wps("PyCaputo PECE", df[:,2], df[:,3]))



################################
################################
################################

# Stiff benchmark dataset collection.
wps_set2 = Any[]


# Solver metadata for the stiff panel.
solvers_all_stiff = [   
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl BDF", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl NewtonGregory", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl Trapezoid", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PITrap", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PIRect", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl KernelCompression", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl SoE", )
    (; pkg = :MATLAB,                           name = "MATLAB BDF", )
    (; pkg = :MATLAB,                           name = "MATLAB NewtonGregory", )
    (; pkg = :MATLAB,                           name = "MATLAB Trapezoid", )
    (; pkg = :MATLAB,                           name = "MATLAB PITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB PIRect", )
];

##### Julia #####
# Load Julia stiff benchmark outputs.
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_BDF.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl BDF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_NewtonGregory.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl NewtonGregory", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_Trapzoid.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl Trapezoid", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PITrap.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PIRect.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl PIRect", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_KernelCompression.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl KernelCompression", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_SoE.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl SoE", df[:,1], df[:,2]))

##### MATLAB #####
# Load MATLAB stiff benchmark outputs.
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_BDF.csv"))
push!(wps_set2, wps("MATLAB BDF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_NewtonGregory.csv"))
push!(wps_set2, wps("MATLAB NewtonGregory", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_Trapzoid.csv"))
push!(wps_set2, wps("MATLAB Trapezoid", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_PITrap.csv"))
push!(wps_set2, wps("MATLAB PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_PIRect.csv"))
push!(wps_set2, wps("MATLAB PIRect", df[:,1], df[:,2]))


# Build a two-panel figure: top = non-stiff, bottom = stiff.

    # Map package families to line styles for visual consistency.
fig = begin
    LINESTYLES = Dict(:FdeSolvers => :dash, :FractionalDiffEq => :solid, :MATLAB => :dot, :Python => :dashdot)
    PKG_ORDER = Dict(:FractionalDiffEq => 1, :MATLAB => 2, :FdeSolvers => 3, :Python => 4)
    legend_order(solvers) = sortperm(eachindex(solvers); by = i -> (get(PKG_ORDER, solvers[i].pkg, typemax(Int)), i))
    ASPECT_RATIO = 0.7
    WIDTH = 1200
    HEIGHT = round(Int, WIDTH * ASPECT_RATIO)
    STROKEWIDTH = 2.5
    # Use categorical palettes with one color per solver curve.

    colors_nonstiff = cgrad(:seaborn_bright, length(solvers_all_nonstiff); categorical = true)

    # Marker cycle so each curve has a distinct marker in addition to color/style.
    colors_stiff = cgrad(:seaborn_bright, length(solvers_all_stiff); categorical = true)
    cycle = Cycle([:marker], covary = true)
    plot_theme = Theme(Lines = (; cycle), Scatter = (; cycle))

    with_theme(plot_theme) do 
        # Top axis: non-stiff benchmark (log-log time vs error profile).
        fig = Figure(; size = (WIDTH, HEIGHT))
        ax = Axis(fig[1, 1], ylabel = L"Time $\mathbf{(s)}$",
            xlabelsize = 22, ylabelsize = 22,
            xlabel = L"Error: $\mathbf{||u-u^\ast||_\infty}$",
            xscale = log10, yscale = log10, xtickwidth = STROKEWIDTH,
            ytickwidth = STROKEWIDTH, spinewidth = STROKEWIDTH,
            xticklabelsize = 20, yticklabelsize = 20)
        # Sort by median runtime so legend/curves follow speed ranking.

        idxs = sortperm(median.(getfield.(wps_set1, :times)))

        ls, scs = [], []
        # Plot all non-stiff solver curves as line + marker overlays.

        for (i, (wp, solver)) in enumerate(zip(wps_set1[idxs], solvers_all_nonstiff[idxs]))
            (; name, times, errors) = wp
            #errors = [err.l∞ for err in errors]
            l = lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors_nonstiff[i])
            sc = CairoMakie.scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors_nonstiff[i])
            push!(ls, l)
            push!(scs, sc)
        end

        # Axis ranges tuned for non-stiff benchmark data spread.
        CairoMakie.xlims!(ax; high=1e2)
        CairoMakie.ylims!(ax; low=10^(-4.3), high=10^(-1.5))

        # Right-side legend for non-stiff panel.
        plotted_solvers = solvers_all_nonstiff[idxs]
        legend_idxs = legend_order(plotted_solvers)
        Legend(fig[1,2], [[ls[i], scs[i]] for i in legend_idxs],
            [plotted_solvers[i].name for i in legend_idxs], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[0, 1] = Label(fig, "Non-stiff FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)

        ############ bottom plot ############
        # Bottom axis: stiff benchmark (same visual grammar, different ranges).
        ax = Axis(fig[3, 1], ylabel = L"Time $\mathbf{(s)}$",
            xlabelsize = 22, ylabelsize = 22,
            xlabel = L"Error: $\mathbf{||u-u^\ast||_\infty}$",
            xscale = log10, yscale = log10, xtickwidth = STROKEWIDTH,
            ytickwidth = STROKEWIDTH, spinewidth = STROKEWIDTH,
            xticklabelsize = 20, yticklabelsize = 20)

        # Sort stiff solvers by median runtime for comparable ordering.
        idxs = sortperm(median.(getfield.(wps_set2, :times)))

        ls, scs = [], []

        # Plot stiff solver curves.
        for (i, (wp, solver)) in enumerate(zip(wps_set2[idxs], solvers_all_stiff[idxs]))
            (; name, times, errors) = wp
            l = lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors_stiff[i])
            sc = CairoMakie.scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors_stiff[i])
            push!(ls, l)
            push!(scs, sc)
        end

        # Axis ranges tuned for stiff benchmark data spread.
        CairoMakie.xlims!(ax; high=10^(-0.5))
        CairoMakie.ylims!(ax; low=10^(-4.4), high=10^(-2.0))

        # Right-side legend for stiff panel.
        plotted_solvers = solvers_all_stiff[idxs]
        legend_idxs = legend_order(plotted_solvers)
        Legend(fig[3,2], [[ls[i], scs[i]] for i in legend_idxs],
            [plotted_solvers[i].name for i in legend_idxs], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[2, 1] = Label(fig, "Stiff FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)
        fig
    end
end

# Export combined non-stiff/stiff performance profile figure.
save("fode_general_benchmarks.svg", fig)
