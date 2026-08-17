"""
    Global sensitivity analysis of steady oscillation amplitudes

The uncertain forcing frequency makes terminal displacement and pointwise
trajectory summaries strongly phase dependent. This script therefore uses
phase-insensitive quantities of interest computed over the final complete
forcing cycles of each trajectory:

1. multi-cycle RMS displacement;
2. half peak-to-peak amplitude.

The Sobol figure reports the first quantity. Set `FODE_SOBOL_SAMPLES` to a
smaller value only for smoke tests; the publication default remains 5000.
"""

using CairoMakie
using CSV
using DataFrames
using FractionalDiffEq
using GlobalSensitivity
using QuasiMonteCarlo
using Random
using Statistics

const OUTPUT_DIR = @__DIR__
const SOLVER_DT = parse(Float64, get(ENV, "FODE_DT", "0.01"))
const SENSITIVITY_TMAX = parse(Float64, get(ENV, "FODE_SENSITIVITY_TMAX", "60.0"))
const STEADY_CYCLES = parse(Int, get(ENV, "FODE_STEADY_CYCLES", "5"))
const SOBOL_SAMPLES = parse(Int, get(ENV, "FODE_SOBOL_SAMPLES", "5000"))
const SOBOL_SEED = parse(Int, get(ENV, "FODE_SOBOL_SEED", "20250017"))
const OUTPUT_STEM = get(ENV, "FODE_SENSITIVITY_OUTPUT_STEM", "sobol_sensitivity_indices_rms")

"""Extract the final `cycles` complete forcing periods from a trajectory."""
function final_cycle_indices(t, ω; cycles = STEADY_CYCLES)
    period = 2π / ω
    window_start = last(t) - cycles * period
    window_start >= first(t) || error(
        "The simulation horizon must contain at least $(cycles) forcing periods",
    )
    return findall(>=(window_start), t)
end

"""
Evaluate phase-insensitive steady-response quantities of interest.

The parameter vector is `[m, c, α, k, F₀, ω]`. The returned quantities are
`[steady multi-cycle RMS, steady half peak-to-peak amplitude]`.
"""
function viscoelastic_qoi(p)
    m, c, α, k, F0, ω = p
    tspan = (0.0, SENSITIVITY_TMAX)

    prob = MultiTermsFODEProblem(
        [m, c, k],
        [2.0, α, 0.0],
        (u, p, t) -> F0 * sin(ω * t),
        [0.0, 0.0],
        tspan,
    )

    sol = solve(prob, MTPITrap(), dt = SOLVER_DT)
    x = first(sol.u) isa Number ? collect(sol.u) : [u[1] for u in sol.u]
    idx = final_cycle_indices(sol.t, ω)
    x_steady = @view x[idx]

    q_rms_steady = sqrt(mean(abs2, x_steady))
    q_amplitude_steady = (maximum(x_steady) - minimum(x_steady)) / 2

    return [q_rms_steady, q_amplitude_steady]
end

const LOWER_BOUNDS = [0.8, 0.5, 1.1, 8.0, 0.8, 1.0]
const UPPER_BOUNDS = [1.2, 8.0, 1.9, 12.0, 1.2, 1.5]

"""Normalize parameter-major and QoI-major Sobol matrix layouts."""
function orient_sensitivity_matrix(matrix, nparams)
    size(matrix, 1) == nparams && return matrix
    size(matrix, 2) == nparams && return permutedims(matrix)
    error("Unexpected sensitivity matrix size: $(size(matrix))")
end

function main()
    minimum_periods = SENSITIVITY_TMAX / (2π / LOWER_BOUNDS[6])
    minimum_periods >= STEADY_CYCLES || error(
        "FODE_SENSITIVITY_TMAX is too short for FODE_STEADY_CYCLES",
    )
    SOBOL_SAMPLES >= 8 || error("FODE_SOBOL_SAMPLES must be at least 8")

    println(
        "Computing Sobol indices from $(SOBOL_SAMPLES) samples; " *
        "QoIs use the final $(STEADY_CYCLES) complete forcing periods ...",
    )
    # Construct two independently shifted Sobol design matrices. Calling the
    # range-based convenience API would split one deterministic high-dimensional
    # Sobol set into A and B, which does not preserve the intended QMC design.
    sampler = SobolSample(R = Shift(rng = MersenneTwister(SOBOL_SEED)))
    design_a, design_b = QuasiMonteCarlo.generate_design_matrices(
        SOBOL_SAMPLES,
        LOWER_BOUNDS,
        UPPER_BOUNDS,
        sampler,
        2,
    )
    sobol_result = gsa(
        viscoelastic_qoi,
        Sobol(),
        design_a,
        design_b,
    )

    parameter_names = ["m", "c", "alpha", "k", "F0", "omega"]
    parameter_labels = [L"m", L"c", L"\alpha", L"k", L"F_0", L"\omega"]
    first_order = orient_sensitivity_matrix(sobol_result.S1, length(parameter_labels))
    total_order = orient_sensitivity_matrix(sobol_result.ST, length(parameter_labels))

    # Column 1 corresponds to the steady multi-cycle RMS QoI. Column 2 is the
    # steady half peak-to-peak amplitude and is retained for cross-checking.
    results = DataFrame(
        parameter = parameter_names,
        plot_label = parameter_labels,
        S1 = first_order[:, 1],
        ST = total_order[:, 1],
    )
    sort!(results, :ST, rev = true)

    fig = Figure(size = (760, 440), fontsize = 16)
    ax = Axis(
        fig[1, 1],
        xlabel = "Parameter",
        ylabel = "Sobol sensitivity index",
        title = "Global sensitivity of steady multi-cycle RMS displacement",
    )

    positions = 1:nrow(results)
    barplot!(
        ax,
        positions .- 0.18,
        results.S1,
        width = 0.35,
        color = :steelblue,
        label = L"First-order index $S_1$",
    )
    barplot!(
        ax,
        positions .+ 0.18,
        results.ST,
        width = 0.35,
        color = :darkorange,
        label = L"Total-order index $S_T$",
    )
    ax.xticks = (positions, results.plot_label)
    axislegend(ax, position = :rt, framevisible = true)

    output_pdf = joinpath(OUTPUT_DIR, "$(OUTPUT_STEM).pdf")
    output_png = joinpath(OUTPUT_DIR, "$(OUTPUT_STEM).png")
    output_csv = joinpath(OUTPUT_DIR, "$(OUTPUT_STEM).csv")
    save(output_pdf, fig)
    save(output_png, fig, px_per_unit = 2)
    CSV.write(output_csv, select(results, :parameter, :S1, :ST))

    println(results)
    println("Saved $(output_pdf)")
    println("Saved $(output_png)")
    println("Saved $(output_csv)")

    return results, fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
