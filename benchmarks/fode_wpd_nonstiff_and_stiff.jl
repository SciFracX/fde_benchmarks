using DataFrames
using Statistics
using CSV
using CairoMakie

struct wps
    name
    times
    errors
end
wps_set1 = Any[]


solvers_all_nonstiff = [
    (; pkg = :FdeSolvers,                       name = "FdeSolver.jl PECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PECE", )    
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB PIEX", )
    (; pkg = :MATLAB,                           name = "MATLAB PECE", )
    (; pkg = :MATLAB,                           name = "MATLAB PIRect", )
    (; pkg = :MATLAB,                           name = "MATLAB PITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB FOTF", )
    (; pkg = :Python,                           name = "PyCaputo PECE", )
];

##### Julia #####
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FdeSolver_PECE.csv"))
push!(wps_set1, wps("FdeSolver.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PECE.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl PECE", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PITrap.csv"))
push!(wps_set1, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

##### MATLAB #####
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
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/PYCAPUTO_PECE.csv"))
push!(wps_set1, wps("PyCaputo PECE", df[:,2], df[:,3]))



################################
################################
################################

wps_set2 = Any[]


solvers_all_stiff = [   
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl BDF", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl NewtonGregory", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl Trapzoid", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PITrap", )
    (; pkg = :FractionalDiffEq,                 name = "FractionalDiffEq.jl PIRect", )
    (; pkg = :MATLAB,                           name = "MATLAB BDF", )
    (; pkg = :MATLAB,                           name = "MATLAB NewtonGregory", )
    (; pkg = :MATLAB,                           name = "MATLAB Trapzoid", )
    (; pkg = :MATLAB,                           name = "MATLAB PITrap", )
    (; pkg = :MATLAB,                           name = "MATLAB PIRect", )
];

##### Julia #####
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_BDF.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl BDF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_NewtonGregory.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl NewtonGregory", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_Trapzoid.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl Trapzoid", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PITrap.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PIRect.csv"))
push!(wps_set2, wps("FractionalDiffEq.jl PIRect", df[:,1], df[:,2]))

##### MATLAB #####
df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_BDF.csv"))
push!(wps_set2, wps("MATLAB BDF", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_NewtonGregory.csv"))
push!(wps_set2, wps("MATLAB NewtonGregory", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_Trapzoid.csv"))
push!(wps_set2, wps("MATLAB Trapzoid", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_PITrap.csv"))
push!(wps_set2, wps("MATLAB PITrap", df[:,1], df[:,2]))

df = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_MATLAB_PIRect.csv"))
push!(wps_set2, wps("MATLAB PIRect", df[:,1], df[:,2]))



fig = begin
    LINESTYLES = Dict(:FdeSolvers => :dash, :FractionalDiffEq => :solid, :MATLAB => :dot, :Python => :dashdot)
    ASPECT_RATIO = 0.7
    WIDTH = 1200
    HEIGHT = round(Int, WIDTH * ASPECT_RATIO)
    STROKEWIDTH = 2.5

    colors_nonstiff = cgrad(:seaborn_bright, length(solvers_all_nonstiff); categorical = true)
    colors_stiff = cgrad(:seaborn_bright, length(solvers_all_stiff); categorical = true)
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

        idxs = sortperm(median.(getfield.(wps_set1, :times)))

        ls, scs = [], []

        for (i, (wp, solver)) in enumerate(zip(wps_set1[idxs], solvers_all_nonstiff[idxs]))
            (; name, times, errors) = wp
            #errors = [err.l∞ for err in errors]
            l = lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors_nonstiff[i])
            sc = scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors_nonstiff[i])
            push!(ls, l)
            push!(scs, sc)
        end

        xlims!(ax; high=1e1)
        ylims!(ax; low=10^(-3.7), high=10^(-1.5))

        Legend(fig[1,2], [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all_nonstiff[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[0, 1] = Label(fig, "Non-stiff FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)

        ############ bottom plot ############
        ax = Axis(fig[3, 1], ylabel = L"Time $\mathbf{(s)}$",
            xlabelsize = 22, ylabelsize = 22,
            xlabel = L"Error: $\mathbf{||u-u^\ast||^2}$",
            xscale = log10, yscale = log10, xtickwidth = STROKEWIDTH,
            ytickwidth = STROKEWIDTH, spinewidth = STROKEWIDTH,
            xticklabelsize = 20, yticklabelsize = 20)

        idxs = sortperm(median.(getfield.(wps_set2, :times)))

        ls, scs = [], []

        for (i, (wp, solver)) in enumerate(zip(wps_set2[idxs], solvers_all_stiff[idxs]))
            (; name, times, errors) = wp
            l = lines!(ax, errors, times; linestyle = LINESTYLES[solver.pkg], label = name,
                linewidth = 5, color = colors_stiff[i])
            sc = scatter!(ax, errors, times; label = name, markersize = 16, strokewidth = 2,
                color = colors_stiff[i])
            push!(ls, l)
            push!(scs, sc)
        end

        xlims!(ax; high=10^(-0.5))
        ylims!(ax; low=10^(-4.4), high=10^(-2.3))

        Legend(fig[3,2], [[l, sc] for (l, sc) in zip(ls, scs)],
            [solver.name for solver in solvers_all_stiff[idxs]], "FODE Solvers";
            framevisible=true, framewidth = STROKEWIDTH, position = :rb,
            titlesize = 20, labelsize = 16, patchsize = (40.0f0, 20.0f0))

        fig[2, 1] = Label(fig, "Stiff FODE Benchmark",
            fontsize = 24, tellwidth = false, font = :bold)
        fig
    end
end

save("fode_general_benchmarks.svg", fig)