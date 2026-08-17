"""
    Long-horizon uncertainty analysis for the fractional Bagley-Torvik oscillator

This script addresses two features that are hidden by pointwise summaries of the
raw displacement trajectories:

1. Under persistent harmonic forcing, the response approaches a nonzero periodic
   state. Damping is therefore assessed from the transient residual after the
   fitted steady-periodic response is removed.
2. When the forcing frequency is uncertain, the trajectories lose phase
   coherence. Pointwise means and medians of displacement then approach zero by
   cancellation and are not representative trajectories. The uncertainty plot
   therefore uses a phase-insensitive, one-cycle moving RMS amplitude.

The uncertain parameters are treated as mutually independent uniform variables
over the ranges specified below. Environment variables may be used for quick
convergence studies without editing the publication defaults; for example,
`FODE_UQ_SAMPLES=32 FODE_OUTPUT_STEM=fractional_uncertainty_smoke`.
"""

using CairoMakie
using CSV
using DataFrames
using FractionalDiffEq
using QuasiMonteCarlo
using Statistics

const OUTPUT_DIR = @__DIR__
const SOLVER_DT = parse(Float64, get(ENV, "FODE_DT", "0.01"))
const OUTPUT_DT = parse(Float64, get(ENV, "FODE_OUTPUT_DT", "0.05"))
const ALPHA_TMAX = parse(Float64, get(ENV, "FODE_ALPHA_TMAX", "80.0"))
const UQ_TMAX = parse(Float64, get(ENV, "FODE_UQ_TMAX", "40.0"))
const N_UQ = parse(Int, get(ENV, "FODE_UQ_SAMPLES", "512"))
const OUTPUT_STEM = get(ENV, "FODE_OUTPUT_STEM", "fractional_uncertainty")
const FIT_CYCLES = 5

const LOWER_BOUNDS = [0.8, 0.5, 1.1, 8.0, 0.8, 1.0]
const UPPER_BOUNDS = [1.2, 8.0, 1.9, 12.0, 1.2, 1.5]

"""Linearly resample a solution onto `saveat`."""
function linear_resample(t, u, saveat)
    y = similar(saveat, Float64)
    k = 1

    for (j, tj) in pairs(saveat)
        while k < length(t) - 1 && t[k + 1] < tj
            k += 1
        end

        if tj <= t[1]
            y[j] = u[1]
        elseif tj >= t[end]
            y[j] = u[end]
        else
            θ = (tj - t[k]) / (t[k + 1] - t[k])
            y[j] = (1 - θ) * u[k] + θ * u[k + 1]
        end
    end

    return y
end

"""
    simulate_fractional_oscillator(p, saveat; dt=SOLVER_DT)

Solve

    m*x'' + c*D^α*x + k*x = F₀*sin(ω*t),  x(0)=x'(0)=0,

on the interval defined by `saveat` and return the displacement on that grid.
"""
function simulate_fractional_oscillator(p, saveat; dt = SOLVER_DT)
    m, c, α, k, F0, ω = p
    tspan = (first(saveat), last(saveat))

    prob = MultiTermsFODEProblem(
        [m, c, k],
        [2.0, α, 0.0],
        (u, p, t) -> F0 * sin(ω * t),
        [0.0, 0.0],
        tspan,
    )

    sol = solve(prob, MTPITrap(), dt = dt)
    x_raw = first(sol.u) isa Number ? collect(sol.u) : [u[1] for u in sol.u]
    return linear_resample(sol.t, x_raw, saveat)
end

"""
Fit `a*sin(ωt) + b*cos(ωt)` over the last `cycles` forcing periods.

The fitted periodic response is used only to separate the persistent forced
response from the decaying transient; it is not interpreted as a mean curve.
"""
function fit_steady_periodic_response(t, x, ω; cycles = FIT_CYCLES)
    period = 2π / ω
    fit_start = max(first(t), last(t) - cycles * period)
    idx = findall(>=(fit_start), t)
    length(idx) >= 3 || error("Not enough points to fit the periodic response")

    design = hcat(sin.(ω .* t[idx]), cos.(ω .* t[idx]))
    coefficients = design \ x[idx]
    fitted = coefficients[1] .* sin.(ω .* t) .+ coefficients[2] .* cos.(ω .* t)
    amplitude = hypot(coefficients[1], coefficients[2])
    phase = atan(coefficients[2], coefficients[1])

    return fitted, amplitude, phase
end

"""Return local maxima of `abs.(residual)` for a transient-envelope plot."""
function transient_peak_envelope(t, residual; threshold = 1.0e-10)
    magnitude = abs.(residual)
    peak_indices = [
        i for i in 2:(length(magnitude) - 1) if
        magnitude[i] >= magnitude[i - 1] &&
        magnitude[i] >= magnitude[i + 1] &&
        magnitude[i] > threshold
    ]
    return t[peak_indices], magnitude[peak_indices]
end

"""
Compute a trailing one-forcing-period RMS amplitude.

The output is `NaN` until a complete forcing period is available. Because each
sample uses its own period `2π/ω`, the statistic is insensitive to phase and
comparable across the uncertain forcing frequencies.
"""
function one_cycle_rms(x, t, ω)
    length(t) == length(x) || throw(DimensionMismatch("t and x must have equal length"))
    length(t) >= 2 || error("At least two time points are required")

    output_dt = t[2] - t[1]
    window = max(2, round(Int, (2π / ω) / output_dt) + 1)
    result = fill(NaN, length(x))
    squared = abs2.(x)
    cumulative = cumsum(squared)

    for j in window:length(x)
        first_index = j - window + 1
        window_sum = cumulative[j] - (first_index > 1 ? cumulative[first_index - 1] : 0.0)
        result[j] = sqrt(window_sum / window)
    end

    return result
end

"""Compute pointwise empirical quantiles, ignoring leading `NaN` values."""
function pointwise_quantiles(values, probabilities)
    summaries = fill(NaN, length(probabilities), size(values, 1))

    for j in axes(values, 1)
        valid = filter(isfinite, @view values[j, :])
        isempty(valid) && continue
        summaries[:, j] .= quantile(valid, probabilities)
    end

    return summaries
end

function main()
    ALPHA_TMAX > 2π / 1.2 || error("FODE_ALPHA_TMAX must cover at least one forcing period")
    UQ_TMAX > 2π / LOWER_BOUNDS[6] ||
        error("FODE_UQ_TMAX must cover at least one period of every sampled frequency")
    N_UQ >= 8 || error("FODE_UQ_SAMPLES must be at least 8")

    alpha_saveat = collect(0.0:OUTPUT_DT:ALPHA_TMAX)
    uq_saveat = collect(0.0:OUTPUT_DT:UQ_TMAX)

    p_nominal = [1.0, 3.0, 1.5, 10.0, 1.0, 1.2]
    α_values = [1.1, 1.3, 1.5, 1.7, 1.9]

    alpha_responses = Vector{Vector{Float64}}(undef, length(α_values))
    transient_peaks = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(α_values))
    steady_amplitudes = zeros(length(α_values))

    println("Computing long-horizon fractional-order sweep to t = $(ALPHA_TMAX) ...")
    for (i, α) in pairs(α_values)
        p = copy(p_nominal)
        p[3] = α
        x = simulate_fractional_oscillator(p, alpha_saveat)
        steady, amplitude, phase = fit_steady_periodic_response(alpha_saveat, x, p[6])

        alpha_responses[i] = x
        transient_peaks[i] = transient_peak_envelope(alpha_saveat, x .- steady)
        steady_amplitudes[i] = amplitude
        println("  α = $(α): fitted steady amplitude = $(round(amplitude, digits = 6)), " *
                "phase = $(round(phase, digits = 6)) rad")
    end

    println("Computing $(N_UQ) Sobol quasi-Monte Carlo trajectories to t = $(UQ_TMAX) ...")
    unit_samples = QuasiMonteCarlo.sample(N_UQ, 6, SobolSample())
    samples = [
        LOWER_BOUNDS[j] + unit_samples[j, i] * (UPPER_BOUNDS[j] - LOWER_BOUNDS[j])
        for j in 1:6, i in 1:N_UQ
    ]

    cycle_rms = fill(NaN, length(uq_saveat), N_UQ)
    for i in 1:N_UQ
        x = simulate_fractional_oscillator(@view(samples[:, i]), uq_saveat)
        cycle_rms[:, i] .= one_cycle_rms(x, uq_saveat, samples[6, i])
        if i % max(1, N_UQ ÷ 8) == 0 || i == N_UQ
            println("  completed $(i)/$(N_UQ) trajectories")
        end
    end

    probabilities = [0.05, 0.25, 0.50, 0.75, 0.95]
    q05, q25, q50, q75, q95 = eachrow(pointwise_quantiles(cycle_rms, probabilities))
    width50 = q75 .- q25
    width90 = q95 .- q05

    # Start the ensemble summary only when every sample has completed one full
    # forcing period. This prevents the early curve from changing population.
    summary_start = 2π / minimum(samples[6, :])
    summary_indices = findall(>=(summary_start), uq_saveat)

    # Use four equally sized panels so that the diagnostic panels are not
    # visually subordinate to the raw response and uncertainty summaries.
    fig = Figure(size = (1700, 1000), fontsize = 17)
    left_layout = GridLayout()
    right_layout = GridLayout()
    fig[1, 1] = left_layout
    fig[1, 2] = right_layout

    ax_response = Axis(
        left_layout[1, 1],
        xlabel = "Time",
        ylabel = "Displacement",
        title = "(a) Forced responses over an extended time horizon",
    )
    ax_transient = Axis(
        left_layout[2, 1],
        xlabel = "Time",
        ylabel = "Transient peak amplitude",
        title = "(c) Decay of the transient residual",
        yscale = log10,
    )
    ax_rms = Axis(
        right_layout[1, 1],
        xlabel = "Time",
        ylabel = "One-cycle RMS displacement",
        title = "(b) Phase-insensitive propagation of parameter uncertainty",
    )
    ax_width = Axis(
        right_layout[2, 1],
        xlabel = "Time",
        ylabel = "Interval width",
        title = "(d) Empirical quantile-band widths",
    )

    colors = Makie.wong_colors()
    for (i, α) in pairs(α_values)
        color = colors[i]
        lines!(
            ax_response,
            alpha_saveat,
            alpha_responses[i],
            color = color,
            linewidth = 2,
            label = L"\alpha = %$α",
        )
        peak_t, peak_magnitude = transient_peaks[i]
        lines!(
            ax_transient,
            peak_t,
            peak_magnitude,
            color = color,
            linewidth = 2,
            label = L"\alpha = %$α",
        )
        scatter!(ax_transient, peak_t, peak_magnitude, color = color, markersize = 4)
    end
    axislegend(ax_response, position = :rt, framevisible = true)
    axislegend(ax_transient, position = :rt, framevisible = true, labelsize = 14)

    t_summary = uq_saveat[summary_indices]
    band!(
        ax_rms,
        t_summary,
        q05[summary_indices],
        q95[summary_indices],
        color = (:steelblue, 0.14),
        label = "90% empirical uncertainty interval",
    )
    band!(
        ax_rms,
        t_summary,
        q25[summary_indices],
        q75[summary_indices],
        color = (:orange, 0.24),
        label = "50% interquartile interval",
    )
    # Thin boundary curves keep the interval limits readable while allowing the
    # trajectories and grid to remain visible through the lighter fills.
    lines!(
        ax_rms,
        t_summary,
        q05[summary_indices],
        color = (:steelblue4, 0.60),
        linewidth = 1.0,
    )
    lines!(
        ax_rms,
        t_summary,
        q95[summary_indices],
        color = (:steelblue4, 0.60),
        linewidth = 1.0,
    )
    lines!(
        ax_rms,
        t_summary,
        q25[summary_indices],
        color = (:darkorange3, 0.65),
        linewidth = 1.0,
    )
    lines!(
        ax_rms,
        t_summary,
        q75[summary_indices],
        color = (:darkorange3, 0.65),
        linewidth = 1.0,
    )
    lines!(
        ax_rms,
        t_summary,
        q50[summary_indices],
        color = :black,
        linewidth = 2.8,
        label = "Median one-cycle RMS",
    )
    axislegend(ax_rms, position = :rt, framevisible = true, labelsize = 14)
    # RMS is nonnegative. Including zero avoids visually exaggerating the band
    # through a truncated vertical axis while preserving its numerical width.
    ylims!(ax_rms, 0.0, 0.4)

    lines!(
        ax_width,
        t_summary,
        width90[summary_indices],
        color = :steelblue4,
        linewidth = 2.5,
        label = L"Q_{0.95}-Q_{0.05}",
    )
    lines!(
        ax_width,
        t_summary,
        width50[summary_indices],
        color = :darkorange3,
        linewidth = 2.5,
        label = L"Q_{0.75}-Q_{0.25}",
    )
    axislegend(ax_width, position = :rt, framevisible = true, labelsize = 14)
    ylims!(ax_width, 0.0, 1.08 * maximum(width90[summary_indices]))

    rowsize!(left_layout, 1, Relative(0.50))
    rowsize!(left_layout, 2, Relative(0.50))
    rowsize!(right_layout, 1, Relative(0.50))
    rowsize!(right_layout, 2, Relative(0.50))
    rowgap!(left_layout, 12)
    rowgap!(right_layout, 12)
    colgap!(fig.layout, 24)

    output_pdf = joinpath(OUTPUT_DIR, "$(OUTPUT_STEM).pdf")
    output_png = joinpath(OUTPUT_DIR, "$(OUTPUT_STEM).png")
    amplitude_csv = joinpath(OUTPUT_DIR, "fractional_order_steady_amplitudes.csv")
    quantile_csv = joinpath(OUTPUT_DIR, "fractional_cycle_rms_quantiles.csv")
    save(output_pdf, fig)
    save(output_png, fig, px_per_unit = 2)
    CSV.write(
        amplitude_csv,
        DataFrame(alpha = α_values, steady_amplitude = steady_amplitudes),
    )
    CSV.write(
        quantile_csv,
        DataFrame(
            time = t_summary,
            q05 = q05[summary_indices],
            q25 = q25[summary_indices],
            q50 = q50[summary_indices],
            q75 = q75[summary_indices],
            q95 = q95[summary_indices],
            width50 = width50[summary_indices],
            width90 = width90[summary_indices],
        ),
    )

    println("Saved $(output_pdf)")
    println("Saved $(output_png)")
    println("Saved $(amplitude_csv)")
    println("Saved $(quantile_csv)")
    println("Steady amplitudes: " * join(round.(steady_amplitudes, digits = 6), ", "))

    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
