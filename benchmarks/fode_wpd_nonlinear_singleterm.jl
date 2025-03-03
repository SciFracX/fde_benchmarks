using DataFrames
using Statistics
using CSV
using CairoMakie

struct wps
    name
    times
    errors
end
wps_set = Any[]


solvers_all = [
    (; pkg = :FdeSolvers,                       name = "FdeSolver.jl PECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PECE", )  
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PIEX", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PITrap", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PIRect", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl BDF", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl Trapezoid", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl NewtonGregory", )
    (; pkg = :MATLAB,                           name = "MATLAB PIEX", )
    (; pkg = :MATLAB,                           name = "MATLAB PECE", )
    (; pkg = :MATLAB,                           name = "MATLAB PIRect", )
    (; pkg = :MATLAB,                           name = "MATLAB PITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB FOTF", )
    (; pkg = :MATLAB,                           name = "MATLAB BDF", )
    (; pkg = :MATLAB,                           name = "MATLAB Trapezoid", )
    (; pkg = :MATLAB,                           name = "MATLAB NewtonGregory", )
    (; pkg = :Python,                           name = "PyCaputo PECE", )
];

##### Julia #####
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FdeSolver_PECE.csv"))
push!(wps_set, wps("FdeSolver.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_PECE.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_PIEX.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PIEX", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_PITrap.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_PIRect.csv"))
push!(wps_set, wps("FractionalDiffEq.jl PIRect", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_BDF.csv"))
push!(wps_set, wps("FractionalDiffEq.jl BDF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_Trapzoid.csv"))
push!(wps_set, wps("FractionalDiffEq.jl Trapzoid", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Single_FractionalDiffEq_NewtonGregory.csv"))
push!(wps_set, wps("FractionalDiffEq.jl NewtonGregory", df[:,1], df[:,2]))

##### MATLAB #####
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_PIEX.csv"))
push!(wps_set, wps("MATLAB PIEX", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_PECE.csv"))
push!(wps_set, wps("MATLAB PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_PIRect.csv"))
push!(wps_set, wps("MATLAB PIRect", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_PITrap.csv"))
push!(wps_set, wps("MATLAB PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_NLFODE_VEC.csv"))
push!(wps_set, wps("MATLAB FOTF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_BDF.csv"))
push!(wps_set, wps("MATLAB BDF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_Trapzoid.csv"))
push!(wps_set, wps("MATLAB Trapzoid", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_MATLAB_NewtonGregory.csv"))
push!(wps_set, wps("MATLAB NewtonGregory", df[:,1], df[:,2]))

##### Python #####
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Singleterm_PYCAPUTO_PECE.csv"))
push!(wps_set, wps("PyCaputo PECE", df[:,2], df[:,3]))

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
            xlabel = L"Error: $\mathbf{||u-u^\ast||^2}$",
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
            sc = CairoMakie.scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors[i])
            push!(ls, l)
            push!(scs, sc)
        end

        CairoMakie.xlims!(ax; high=1e0)
        CairoMakie.ylims!(ax; low=10^(-5), high=10^(-2))

        Legend(fig[1,2], [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[0, :] = Label(fig, "Single Term FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)
        fig
    end
end

save("singleterm_fode_benchmarks.svg", fig)