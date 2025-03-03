using DataFrames
using Statistics
using CSV
using CairoMakie

struct wps
    name
    times
    errors
end
wps_set = Any[];


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
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPECE.csv"))
push!(wps_set, wps("FdeSolver.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIEX.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPITrap.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIRect.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PIRect", df[:,1], df[:,2]))

##### MATLAB #####
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

fig = begin
    LINESTYLES = Dict(:FdeSolvers => :dash, :FractionalDiffEq => :solid, :MATLAB => :dot, :Python => :dashdot)
    ASPECT_RATIO = 0.7
    WIDTH = 1200
    HEIGHT = round(Int, WIDTH * ASPECT_RATIO)
    STROKEWIDTH = 2.5

    colors = cgrad(:seaborn_bright, length(solvers_all); categorical = true)
    cycle = Cycle([:marker], covary = true)
    plot_theme = Theme(Lines = (; cycle), Scatter = (; cycle))

    with_theme(plot_theme) do 
        fig = Figure(; size = (WIDTH, HEIGHT))
        ax = Axis(fig[1, 1], ylabel = L"Time $\mathbf{(s)}$",
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

        xlims!(ax; high=10^(1.2))
        ylims!(ax; low=10^(-3), high=10^(1.8))

        axislegend(ax, [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[0, :] = Label(fig, "Linear Multi-terms FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)
        fig
    end
end

save("linear_fode_benchmarks.svg", fig)