using FractionalDiffEq
using GlobalSensitivity
using QuasiMonteCarlo
using Statistics
using DataFrames
using CairoMakie

function viscoelastic_qoi(p)
    m, c, α, k, F0, ω = p

    T = 30.0
    tspan = (0.0, T)
    x0 = 0.0
    v0 = 0.0

    # Example schematic form:
    # m*x'' + c*D^α x + k*x = F0*sin(ω*t)
    #
    # Replace this part with the exact FractionalDiffEq.jl problem constructor
    # used by your package for multi-term FODEs.

    prob = MultiTermsFODEProblem([m, c, k],
        [2.0, α, 0.0],
        (u, p, t) -> F0 * sin(ω * t),
        [x0, v0],
        tspan,
        [0]
    )

    sol = solve(prob, MTPITrap(), dt = 0.01)

    x = [u[1] for u in sol.u]

    q_peak = maximum(abs.(x))

    q_rms = sqrt(mean(abs2, x))

    q_final = abs(x[end])

    return [q_peak, q_rms, q_final]
end

lb = [0.8, 0.5, 1.1, 8.0, 0.8, 1.0]

ub = [1.2, 8.0, 1.9, 12.0, 1.2, 1.5]

bounds = [[lb[i], ub[i]] for i in eachindex(lb)]

sobol_res = gsa(

    viscoelastic_qoi,

    Sobol(),

    bounds;

    samples = 5000

)

param_names = [L"m", L"c", L"\alpha", L"k", L"F_0", L"\omega"]

function orient_sensitivity_matrix(matrix, nparams)

    size(matrix, 1) == nparams && return matrix

    size(matrix, 2) == nparams && return permutedims(matrix)

    error("Unexpected sensitivity matrix size: $(size(matrix))")

end

S1 = orient_sensitivity_matrix(sobol_res.S1, length(param_names))

ST = orient_sensitivity_matrix(sobol_res.ST, length(param_names))

# Use RMS response instead of peak displacement for a smoother QoI.

df = DataFrame(

    parameter = param_names,

    S1 = S1[:, 2],

    ST = ST[:, 2],

)

sort!(df, :ST, rev = true)

fig = Figure(size = (720, 420))

ax = Axis(

    fig[1, 1],

    xlabel = "Parameter",

    ylabel = "Sobol sensitivity index",

    title = "Global sensitivity analysis of RMS displacement"

)

x = 1:nrow(df)

barplot!(

    ax,

    x .- 0.18,

    df.S1,

    width = 0.35,

    label = L"First-order index $S_1$"

)

barplot!(

    ax,

    x .+ 0.18,

    df.ST,

    width = 0.35,

    label = L"Total-order index $S_t$"

)

ax.xticks = (x, df.parameter)

axislegend(ax, position = :rt)

save("sobol_sensitivity_indices_rms.pdf", fig)

save("sobol_sensitivity_indices_rms.png", fig)