using DataFrames
using Statistics
using CSV
using CairoMakie

# Container for one solver benchmark curve.
# `times` and `errors` are vectors read from benchmark CSV files.
struct wps
    name
    times
    errors
end

# Collection of all curves used in this figure.
wps_set = Any[];


# Solver metadata used for legend text and linestyle grouping.
# The order here should match the order of pushed datasets in `wps_set`.
solvers_all = [
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl MTPECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl MTPIEX", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl MTPITrap", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl MTPIRect", )
    (; pkg = :MATLAB,                           name = "MATLAB MTPECE", )
    (; pkg = :MATLAB,                           name = "MATLAB MTPIEX", )
    (; pkg = :MATLAB,                           name = "MATLAB MTPITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB MTPIRect", )
    (; pkg = :MATLAB,                           name = "MATLAB FOTF", )
    (; pkg = :MATLAB,                           name = "MATLAB Matrix Discretization", )
];

##### Julia #####
# Load Julia multi-term benchmark outputs (columns: time, error).
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPECE.csv"))
push!(wps_set, wps("FractionalDiffEq.jl MTPECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIEX.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPITrap.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIRect.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PIRect", df[:,1], df[:,2]))

##### MATLAB #####
# Load MATLAB multi-term benchmark outputs.
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_MTPECE.csv"))
push!(wps_set, wps("PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_MTPIEX.csv"))
push!(wps_set, wps("PIEX", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_MTPITrap.csv"))
push!(wps_set, wps("PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_MTPIRect.csv"))
push!(wps_set, wps("PIRect", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_MTCAPUTO9.csv"))
push!(wps_set, wps("FOTF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/MATLAB_Matrix.csv"))
push!(wps_set, wps("Matrix Discretization", df[:,1], df[:,2]))

# Build a single-panel performance profile (time vs error in log-log scale).
fig = begin
    # Package family -> line style mapping.
    LINESTYLES = Dict(:FdeSolvers => :dash, :FractionalDiffEq => :solid, :MATLAB => :dot, :Python => :dashdot)
    ASPECT_RATIO = 0.7
    WIDTH = 1200
    HEIGHT = round(Int, WIDTH * ASPECT_RATIO)
    STROKEWIDTH = 2.5

    # Categorical palette gives each solver a distinct color.
    colors = cgrad(:seaborn_bright, length(solvers_all); categorical = true)

    # Add marker variation in addition to color/linestyle.
    cycle = Cycle([:marker], covary = true)
    plot_theme = Theme(Lines = (; cycle), Scatter = (; cycle))

    with_theme(plot_theme) do 
        fig = Figure(; size = (WIDTH, HEIGHT))
        # Axis uses log-log scaling to emphasize efficiency/accuracy trade-offs.
        ax = Axis(fig[1, 1], ylabel = L"Time $\mathbf{(s)}$",
            xlabelsize = 22, ylabelsize = 22,
            xlabel = L"Error: $\mathbf{||u-u^\ast||_2}$",
            xscale = log10, yscale = log10, xtickwidth = STROKEWIDTH,
            ytickwidth = STROKEWIDTH, spinewidth = STROKEWIDTH,
            xticklabelsize = 20, yticklabelsize = 20)

        # Rank curves by median runtime so faster methods appear first.
        idxs = sortperm(median.(getfield.(wps_set, :times)))

        # Find MATLAB MTPIEX (label contains PIEX but not FractionalDiffEq).
        # This method is manually placed last for a preferred visual/legend order.

        piex_idx = findfirst(i ->
            occursin("PIEX", wps_set[i].name) &&
            !(occursin("FractionalDiffEq", wps_set[i].name)),
            eachindex(wps_set)
        )

        # Move MATLAB PIEX to the end while preserving relative order of others.

        idxs = vcat(filter(!=(piex_idx), idxs), piex_idx)

        ls, scs = [], []

        # Draw each solver as line + scatter points.
        for (i, (wp, solver)) in enumerate(zip(wps_set[idxs], solvers_all[idxs]))
            (; name, times, errors) = wp
            #errors = [err.l∞ for err in errors]
            l = CairoMakie.lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors[i])
            sc = CairoMakie.scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors[i])
            push!(ls, l)
            push!(scs, sc)
        end

        # Axis limits tuned for this multi-term benchmark dataset.
        CairoMakie.xlims!(ax; high=10^(1.2))
        CairoMakie.ylims!(ax; low=10^(-3), high=10^(1.8))

        # Legend combines corresponding line and marker entries per solver.
        axislegend(ax, [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[0, :] = Label(fig, "Linear Multi-terms FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)
        fig
    end
end

# Save final figure for manuscript/report use.
save("linear_fode_benchmarks.svg", fig)