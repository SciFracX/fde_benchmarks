using FractionalDiffEq
using GlobalSensitivity
using QuasiMonteCarlo
using Statistics
using DataFrames
using CairoMakie

# -----------------------------------------------------------------------------
# Quantity of interest (QoI)
# -----------------------------------------------------------------------------
# This script performs a Sobol global sensitivity study for a viscoelastic
# oscillator modeled with a multi-term fractional differential equation. The
# central task is to define a quantity of interest that can be evaluated many
# times for different parameter combinations.
#
# The parameter vector `p` is interpreted as:
#   m   -> mass
#   c   -> fractional damping coefficient
#   α   -> fractional derivative order
#   k   -> linear stiffness
#   F0  -> forcing amplitude
#   ω   -> forcing frequency
#
# For each parameter sample, the model is integrated over a fixed time horizon
# and three scalar response metrics are returned:
#   1. Peak absolute displacement
#   2. Root-mean-square displacement
#   3. Final absolute displacement
#
# In the post-processing section below, the RMS response is used as the main
# QoI because it is generally smoother than the peak value and therefore tends
# to produce more stable sensitivity rankings.
function viscoelastic_qoi(p)
    m, c, α, k, F0, ω = p

    # Simulation horizon. A long enough interval is used to capture both the
    # transient response and the later-time oscillatory behavior.
    T = 30.0
    tspan = (0.0, T)

    # Initial displacement and velocity.
    x0 = 0.0
    v0 = 0.0

    # Schematic model form:
    #   m * x'' + c * D^α x + k * x = F0 * sin(ω * t)
    #
    # The exact constructor below should match the multi-term FODE interface
    # provided by the version of FractionalDiffEq.jl used in this project.
    prob = MultiTermsFODEProblem([
            m,
            c,
            k,
        ],
        [2.0, α, 0.0],
        (u, p, t) -> F0 * sin(ω * t),
        [x0, v0],
        tspan,
        [0]
    )

    # Numerical integration of the fractional system.
    sol = solve(prob, MTPITrap(), dt = 0.01)

    # Extract the displacement history from the state vector.
    x = [u[1] for u in sol.u]

    # Maximum excursion of the displacement trajectory.
    q_peak = maximum(abs.(x))

    # RMS displacement, which captures the overall response amplitude across the
    # full time interval.
    q_rms = sqrt(mean(abs2, x))

    # Terminal displacement magnitude at the final simulation time.
    q_final = abs(x[end])

    return [q_peak, q_rms, q_final]
end

# -----------------------------------------------------------------------------
# Parameter uncertainty bounds
# -----------------------------------------------------------------------------
# Sobol analysis assumes a rectangular parameter domain. The lower and upper
# bounds below define a moderate uncertainty range for each physical parameter.
lb = [0.8, 0.5, 1.1, 8.0, 0.8, 1.0]
ub = [1.2, 8.0, 1.9, 12.0, 1.2, 1.5]

# Convert paired lower and upper limits into the format expected by `gsa`.
bounds = [[lb[i], ub[i]] for i in eachindex(lb)]

# -----------------------------------------------------------------------------
# Global sensitivity analysis
# -----------------------------------------------------------------------------
# Use a fairly large Monte Carlo budget so the estimated Sobol indices are
# stable enough for interpretation and comparison.
sobol_res = gsa(
    viscoelastic_qoi,
    Sobol(),
    bounds;
    samples = 5000,
)

# Human-readable labels for the parameter axis and the reporting table.
param_names = [L"m", L"c", L"\alpha", L"k", L"F_0", L"\omega"]

# Depending on the result type, `GlobalSensitivity.jl` may store the Sobol
# matrices in either parameter-major or QoI-major orientation. This helper
# normalizes the layout so that rows always correspond to parameters.
function orient_sensitivity_matrix(matrix, nparams)
    size(matrix, 1) == nparams && return matrix
    size(matrix, 2) == nparams && return permutedims(matrix)

    error("Unexpected sensitivity matrix size: $(size(matrix))")
end

# First-order indices measure the direct effect of each parameter alone.
S1 = orient_sensitivity_matrix(sobol_res.S1, length(param_names))

# Total-order indices measure the direct effect plus all interaction effects.
ST = orient_sensitivity_matrix(sobol_res.ST, length(param_names))

# Assemble a compact table for the RMS response.
# Column 2 is selected because the QoI vector is ordered as:
#   [peak displacement, RMS displacement, final displacement]
df = DataFrame(
    parameter = param_names,
    S1 = S1[:, 2],
    ST = ST[:, 2],
)

# Sort parameters by total-order importance so the strongest contributors appear
# first in the final plot.
sort!(df, :ST, rev = true)

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------
# The final figure compares first-order and total-order Sobol indices. The gap
# between the two bars indicates how much interaction with other parameters is
# contributing to the total influence of each variable.
fig = Figure(size = (720, 420))

ax = Axis(
    fig[1, 1],
    xlabel = "Parameter",
    ylabel = "Sobol sensitivity index",
    title = "Global sensitivity analysis of RMS displacement",
)

# Use the sorted ranking to place parameters on the x-axis.
x = 1:nrow(df)

barplot!(
    ax,
    x .- 0.18,
    df.S1,
    width = 0.35,
    label = L"First-order index $S_1$",
)

# Overlay the total-order contribution so interaction effects are visible.
barplot!(
    ax,
    x .+ 0.18,
    df.ST,
    width = 0.35,
    label = L"Total-order index $S_t$",
)

# Replace numeric ticks with parameter labels for readability.
ax.xticks = (x, df.parameter)

axislegend(ax, position = :rt)

# Save both vector and raster outputs
save("sobol_sensitivity_indices_rms.pdf", fig)
save("sobol_sensitivity_indices_rms.png", fig)