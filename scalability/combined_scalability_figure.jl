#!/usr/bin/env julia

"""
Combine the time-step, dimension, and interval scalability results into one
publication-style 3×2 figure.

Run from the project root:

    julia --project=. scalability/combined_scalability_figure.jl

This script only reads existing CSV files. Generate them first with
`memory_fode/memory_figure.jl` and `scalability/scalability_figure.jl`.
"""

using CairoMakie
using CSV
using DataFrames
using Statistics

const SCRIPT_DIR = @__DIR__
const PROJECT_DIR = normpath(joinpath(SCRIPT_DIR, ".."))
const TIME_STEP_RESULTS = joinpath(PROJECT_DIR, "memory_fode", "memory_results.csv")
const SCALABILITY_RESULTS = joinpath(SCRIPT_DIR, "scalability_results.csv")
const FIGURE_PDF = joinpath(SCRIPT_DIR, "combined_scalability.pdf")
const FIGURE_SVG = joinpath(SCRIPT_DIR, "combined_scalability.svg")
const FIGURE_PNG = joinpath(SCRIPT_DIR, "combined_scalability.png")
const TIME_STEP_PROBLEM_ID = "coupled_linear_d10_alpha0.8_t20"
const IMPLEMENTATIONS = ["Julia", "MATLAB", "Python"]

const IMPLEMENTATION_LABELS = Dict(
    "Julia" => "FractionalDiffEq.jl (PECE)",
    "MATLAB" => "MATLAB fde_pi12_pc",
    "Python" => "pycaputo (PECE)",
)
const COLORS = Dict(
    "Julia" => "#4063D8",
    "MATLAB" => "#D95319",
    "Python" => "#2E8B57",
)
const MARKERS = Dict(
    "Julia" => :circle,
    "MATLAB" => :rect,
    "Python" => :utriangle,
)

function require_file(path::String, generating_command::String)
    isfile(path) || error(
        "missing $(path). Generate it first with:\n$(generating_command)"
    )
end

function load_time_step_results()
    require_file(
        TIME_STEP_RESULTS,
        "julia --project=. memory_fode/memory_figure.jl",
    )
    data = DataFrame(CSV.File(TIME_STEP_RESULTS))
    required = [
        :problem_id,
        :implementation,
        :measurement,
        :trial,
        :N,
        :runtime_s,
        :peak_rss_mib,
    ]
    missing_columns = setdiff(required, propertynames(data))
    isempty(missing_columns) || error(
        "$(TIME_STEP_RESULTS) is missing columns: $(join(missing_columns, ", "))"
    )
    data = data[data.problem_id .== TIME_STEP_PROBLEM_ID, :]
    nrow(data) > 0 || error(
        "$(TIME_STEP_RESULTS) has no rows for $(TIME_STEP_PROBLEM_ID)"
    )
    return data
end

function load_scalability_results()
    require_file(
        SCALABILITY_RESULTS,
        "julia --project=. scalability/scalability_figure.jl",
    )
    data = DataFrame(CSV.File(SCALABILITY_RESULTS))
    required = [
        :implementation,
        :measurement,
        :experiment,
        :trial,
        :dimension,
        :final_time,
        :runtime_s,
        :peak_rss_mib,
    ]
    missing_columns = setdiff(required, propertynames(data))
    isempty(missing_columns) || error(
        "$(SCALABILITY_RESULTS) is missing columns: $(join(missing_columns, ", "))"
    )
    return data
end

function summarize_series(
    data::DataFrame,
    implementation::String,
    x_column::Symbol;
    experiment::Union{Nothing, String} = nothing,
)
    selected = (data.implementation .== implementation) .&
               (data.measurement .== "solve")
    if !isnothing(experiment)
        selected .&= data.experiment .== experiment
    end
    rows = data[selected, :]
    x = sort(unique(Float64.(rows[!, x_column])))
    runtime_median = Float64[]
    runtime_q25 = Float64[]
    runtime_q75 = Float64[]
    rss_median = Float64[]
    rss_q25 = Float64[]
    rss_q75 = Float64[]
    for value in x
        group = Float64.(rows[!, x_column]) .== value
        runtimes = Float64.(rows.runtime_s[group])
        rss = Float64.(rows.peak_rss_mib[group])
        push!(runtime_median, median(runtimes))
        push!(runtime_q25, quantile(runtimes, 0.25))
        push!(runtime_q75, quantile(runtimes, 0.75))
        push!(rss_median, median(rss))
        push!(rss_q25, quantile(rss, 0.25))
        push!(rss_q75, quantile(rss, 0.75))
    end
    return (;
        x,
        runtime_median,
        runtime_q25,
        runtime_q75,
        rss_median,
        rss_q25,
        rss_q75,
    )
end

function draw_series!(
    runtime_axis,
    memory_axis,
    summary,
    color,
    marker;
    show_band::Bool = true,
)
    if show_band
        band!(
            runtime_axis,
            summary.x,
            summary.runtime_q25,
            summary.runtime_q75;
            color = (color, 0.13),
        )
    end
    runtime_line = lines!(
        runtime_axis,
        summary.x,
        summary.runtime_median;
        color,
        linewidth = 3.0,
    )
    runtime_points = scatter!(
        runtime_axis,
        summary.x,
        summary.runtime_median;
        color,
        marker,
        markersize = 12,
        strokecolor = :white,
        strokewidth = 1,
    )
    if show_band
        band!(
            memory_axis,
            summary.x,
            summary.rss_q25,
            summary.rss_q75;
            color = (color, 0.13),
        )
    end
    lines!(
        memory_axis,
        summary.x,
        summary.rss_median;
        color,
        linewidth = 3.0,
    )
    scatter!(
        memory_axis,
        summary.x,
        summary.rss_median;
        color,
        marker,
        markersize = 12,
        strokecolor = :white,
        strokewidth = 1,
    )
    return runtime_line, runtime_points
end

function draw_baseline!(
    axis,
    data::DataFrame,
    implementation::String,
    x,
    color;
    show_band::Bool = true,
)
    rows = data[
        (data.implementation .== implementation) .&
        (data.measurement .== "baseline"),
        :,
    ]
    nrow(rows) == 0 && return false
    baseline_median = median(Float64.(rows.peak_rss_mib))
    baseline_q25 = quantile(Float64.(rows.peak_rss_mib), 0.25)
    baseline_q75 = quantile(Float64.(rows.peak_rss_mib), 0.75)
    baseline_x = [first(x), last(x)]
    if show_band
        band!(
            axis,
            baseline_x,
            fill(baseline_q25, 2),
            fill(baseline_q75, 2);
            color = (color, 0.055),
        )
    end
    hlines!(
        axis,
        [baseline_median];
        color = (color, 0.72),
        linestyle = :dash,
        linewidth = 1.8,
    )
    return true
end

function format_n_tick(value::Real)
    value < 1000 && return string(round(Int, value))
    value < 10000 && return string(round(value / 1000; digits = 2), "k")
    return string(round(value / 1000; digits = 1), "k")
end

function make_axis(
    figure,
    row::Int,
    column::Int,
    title::String,
    xlabel::String,
    ylabel::String,
    ticks;
    logarithmic_y::Bool,
)
    grid_color = (:gray35, 0.15)
    return Axis(
        figure[row, column];
        title,
        xlabel,
        ylabel,
        xscale = log10,
        yscale = logarithmic_y ? log10 : identity,
        xticks = ticks,
        xgridvisible = false,
        ygridvisible = true,
        ygridcolor = grid_color,
        ygridwidth = 0.8,
        spinewidth = 1.1,
    )
end

function choose_ticks(values; maximum_ticks::Int = 5)
    length(values) <= maximum_ticks && return values
    indices = unique(round.(Int, range(1, length(values); length = maximum_ticks)))
    return values[indices]
end

function make_figure(time_step_data::DataFrame, scalability_data::DataFrame)
    n_values = sort(unique(time_step_data.N[time_step_data.measurement .== "solve"]))
    dimension_values = sort(unique(scalability_data.dimension[
        (scalability_data.measurement .== "solve") .&
        (scalability_data.experiment .== "dimension")
    ]))
    interval_values = sort(unique(scalability_data.final_time[
        (scalability_data.measurement .== "solve") .&
        (scalability_data.experiment .== "interval")
    ]))
    isempty(n_values) && error("time-step results contain no solve measurements")
    isempty(dimension_values) && error("dimension results contain no solve measurements")
    isempty(interval_values) && error("interval results contain no solve measurements")

    figure = Figure(size = (1800, 860), fontsize = 16)
    displayed_n_values = choose_ticks(n_values)
    n_ticks = (Float64.(displayed_n_values), format_n_tick.(displayed_n_values))
    displayed_dimension_values = choose_ticks(dimension_values)
    displayed_interval_values = choose_ticks(interval_values)
    dimension_ticks = (
        Float64.(displayed_dimension_values),
        string.(displayed_dimension_values),
    )
    interval_ticks = (
        Float64.(displayed_interval_values),
        [string(round(Int, value)) for value in displayed_interval_values],
    )

    n_runtime = make_axis(
        figure,
        1,
        1,
        "(a) Runtime vs number of time steps",
        "Number of time steps, N",
        "Runtime (s)",
        n_ticks,
        logarithmic_y = true,
    )
    d_runtime = make_axis(
        figure,
        1,
        2,
        "(b) Runtime vs system dimension",
        "System dimension, d",
        "",
        dimension_ticks,
        logarithmic_y = true,
    )
    t_runtime = make_axis(
        figure,
        1,
        3,
        "(c) Runtime vs integration interval",
        "Final time, T",
        "",
        interval_ticks,
        logarithmic_y = true,
    )
    n_memory = make_axis(
        figure,
        2,
        1,
        "(d) Peak RSS vs number of time steps",
        "Number of time steps, N",
        "Peak RSS (MiB)",
        n_ticks,
        logarithmic_y = false,
    )
    d_memory = make_axis(
        figure,
        2,
        2,
        "(e) Peak RSS vs system dimension",
        "System dimension, d",
        "",
        dimension_ticks,
        logarithmic_y = false,
    )
    t_memory = make_axis(
        figure,
        2,
        3,
        "(f) Peak RSS vs integration interval",
        "Final time, T",
        "",
        interval_ticks,
        logarithmic_y = false,
    )

    legend_entries = Any[]
    legend_labels = String[]
    has_baseline = false
    for implementation in IMPLEMENTATIONS
        color = COLORS[implementation]
        marker = MARKERS[implementation]
        show_band = implementation != "MATLAB"
        n_summary = summarize_series(time_step_data, implementation, :N)
        d_summary = summarize_series(
            scalability_data,
            implementation,
            :dimension;
            experiment = "dimension",
        )
        t_summary = summarize_series(
            scalability_data,
            implementation,
            :final_time;
            experiment = "interval",
        )
        isempty(n_summary.x) && continue
        isempty(d_summary.x) && continue
        isempty(t_summary.x) && continue

        legend_line, legend_points = draw_series!(
            n_runtime,
            n_memory,
            n_summary,
            color,
            marker,
            show_band = show_band,
        )
        draw_series!(
            d_runtime,
            d_memory,
            d_summary,
            color,
            marker;
            show_band,
        )
        draw_series!(
            t_runtime,
            t_memory,
            t_summary,
            color,
            marker;
            show_band,
        )
        has_baseline |= draw_baseline!(
            n_memory,
            time_step_data,
            implementation,
            n_summary.x,
            color,
            show_band = show_band,
        )
        has_baseline |= draw_baseline!(
            d_memory,
            scalability_data,
            implementation,
            d_summary.x,
            color,
            show_band = show_band,
        )
        has_baseline |= draw_baseline!(
            t_memory,
            scalability_data,
            implementation,
            t_summary.x,
            color,
            show_band = show_band,
        )
        push!(legend_entries, [legend_line, legend_points])
        push!(legend_labels, IMPLEMENTATION_LABELS[implementation])
    end

    if has_baseline
        push!(legend_entries, LineElement(
            color = (:gray30, 0.8),
            linestyle = :dash,
            linewidth = 2.0,
        ))
        push!(legend_labels, "Initialized-process RSS baseline")
    end
    Legend(
        figure[0, 1:3],
        legend_entries,
        legend_labels;
        orientation = :horizontal,
        framevisible = false,
        tellheight = true,
        tellwidth = false,
        patchsize = (38.0f0, 18.0f0),
        labelsize = 14,
    )

    linkyaxes!(n_runtime, d_runtime, t_runtime)
    linkyaxes!(n_memory, d_memory, t_memory)
    ylims!(n_memory; low = 0)
    hideydecorations!(d_runtime; grid = false)
    hideydecorations!(t_runtime; grid = false)
    hideydecorations!(d_memory; grid = false)
    hideydecorations!(t_memory; grid = false)

    colgap!(figure.layout, 30)
    rowgap!(figure.layout, 22)
    save(FIGURE_PDF, figure)
    save(FIGURE_SVG, figure)
    save(FIGURE_PNG, figure; px_per_unit = 2)
    return figure
end

function main()
    time_step_data = load_time_step_results()
    scalability_data = load_scalability_results()
    make_figure(time_step_data, scalability_data)
    println("Figures: $(FIGURE_PDF), $(FIGURE_SVG), $(FIGURE_PNG)")
end

main()
