#!/usr/bin/env julia

"""
Run and plot independent state-dimension and time-interval scalability tests.

Default protocol:
  dimension sweep: d = 5,10,15,20,25,30,35,40; T = 20; N = 5120
  interval sweep:  T = 5,10,15,20,25,30,35,40; d = 10; h = 1/256

Examples (run from the project root):
    julia --project=. scalability/scalability_figure.jl
    julia --project=. scalability/scalability_figure.jl --plot-only
    julia --project=. scalability/scalability_figure.jl --languages=Julia,Python
    julia --project=. scalability/scalability_figure.jl --repetitions=10

Environment overrides:
    SCALABILITY_PYTHON=/path/to/python-with-pycaputo
    SCALABILITY_MATLAB=/path/to/matlab
    SCALABILITY_OUTPUT_DIR=/path/for/csv-and-figures
"""

using CairoMakie
using CSV
using DataFrames
using Printf
using Statistics

const SCALABILITY_DIR = @__DIR__
const PROJECT_DIR = normpath(joinpath(SCALABILITY_DIR, ".."))
const OUTPUT_DIR = abspath(get(ENV, "SCALABILITY_OUTPUT_DIR", SCALABILITY_DIR))
const RESULTS_FILE = joinpath(OUTPUT_DIR, "scalability_results.csv")
const FIGURE_PDF = joinpath(OUTPUT_DIR, "scalability.pdf")
const FIGURE_SVG = joinpath(OUTPUT_DIR, "scalability.svg")
const FIGURE_PNG = joinpath(OUTPUT_DIR, "scalability.png")

const IMPLEMENTATIONS = ["Julia", "MATLAB", "Python"]
const DEFAULT_DIMENSIONS = [5, 10, 15, 20, 25, 30, 35, 40]
const DEFAULT_INTERVALS = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
const DEFAULT_DIMENSION_T = 20.0
const DEFAULT_DIMENSION_N = 5120
const DEFAULT_INTERVAL_DIMENSION = 10
const DEFAULT_STEPS_PER_UNIT = 256
const DEFAULT_REPETITIONS = 5
const DEFAULT_MAX_ATTEMPTS = 3

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

Base.@kwdef mutable struct Options
    run_experiment::Bool = true
    rerun::Bool = false
    implementations::Vector{String} = copy(IMPLEMENTATIONS)
    dimensions::Vector{Int} = copy(DEFAULT_DIMENSIONS)
    intervals::Vector{Float64} = copy(DEFAULT_INTERVALS)
    dimension_t::Float64 = DEFAULT_DIMENSION_T
    dimension_n::Int = DEFAULT_DIMENSION_N
    interval_dimension::Int = DEFAULT_INTERVAL_DIMENSION
    steps_per_unit::Int = DEFAULT_STEPS_PER_UNIT
    repetitions::Int = DEFAULT_REPETITIONS
    max_attempts::Int = DEFAULT_MAX_ATTEMPTS
end

function usage()
    println("""
    usage: julia --project=. scalability/scalability_figure.jl [options]

      --languages=Julia,MATLAB,Python
      --dimensions=D1,D2,...          dimension sweep values
      --intervals=T1,T2,...           final times for interval sweep
      --dimension-t=T                 fixed T for dimension sweep
      --dimension-n=N                 fixed N for dimension sweep
      --interval-dimension=D          fixed d for interval sweep
      --steps-per-unit=K              interval sweep uses h=1/K and N=K*T
      --repetitions=R                 fresh processes per point
      --max-attempts=A                attempts after worker failure
      --plot-only                     plot existing scalability_results.csv
      --rerun                         replace selected measurements
      --help                          show this message
    """)
end

function canonical_implementation(name::AbstractString)
    index = findfirst(==(lowercase(strip(name))), lowercase.(IMPLEMENTATIONS))
    isnothing(index) && error(
        "unknown implementation $(repr(name)); choose from $(join(IMPLEMENTATIONS, ", "))"
    )
    return IMPLEMENTATIONS[index]
end

function parse_int_list(raw::AbstractString, name::String)
    values = parse.(Int, filter(!isempty, strip.(split(raw, ','))))
    isempty(values) && error("$(name) must contain at least one integer")
    all(>(0), values) || error("all $(name) values must be positive")
    length(unique(values)) == length(values) || error("$(name) contains duplicates")
    return sort(values)
end

function parse_float_list(raw::AbstractString, name::String)
    values = parse.(Float64, filter(!isempty, strip.(split(raw, ','))))
    isempty(values) && error("$(name) must contain at least one number")
    all(value -> isfinite(value) && value > 0, values) ||
        error("all $(name) values must be finite and positive")
    length(unique(values)) == length(values) || error("$(name) contains duplicates")
    return sort(values)
end

function parse_options(args)
    options = Options()
    for arg in args
        if arg == "--plot-only"
            options.run_experiment = false
        elseif arg == "--rerun"
            options.rerun = true
        elseif arg == "--help" || arg == "-h"
            usage()
            exit(0)
        elseif startswith(arg, "--languages=")
            raw = split(arg, "="; limit = 2)[2]
            names = filter(!isempty, strip.(split(raw, ',')))
            isempty(names) && error("--languages must contain at least one name")
            options.implementations = unique(canonical_implementation.(names))
        elseif startswith(arg, "--dimensions=")
            options.dimensions = parse_int_list(
                split(arg, "="; limit = 2)[2], "--dimensions"
            )
        elseif startswith(arg, "--intervals=")
            options.intervals = parse_float_list(
                split(arg, "="; limit = 2)[2], "--intervals"
            )
        elseif startswith(arg, "--dimension-t=")
            options.dimension_t = parse(Float64, split(arg, "="; limit = 2)[2])
            isfinite(options.dimension_t) && options.dimension_t > 0 ||
                error("--dimension-t must be finite and positive")
        elseif startswith(arg, "--dimension-n=")
            options.dimension_n = parse(Int, split(arg, "="; limit = 2)[2])
            options.dimension_n > 0 || error("--dimension-n must be positive")
        elseif startswith(arg, "--interval-dimension=")
            options.interval_dimension = parse(Int, split(arg, "="; limit = 2)[2])
            options.interval_dimension > 0 ||
                error("--interval-dimension must be positive")
        elseif startswith(arg, "--steps-per-unit=")
            options.steps_per_unit = parse(Int, split(arg, "="; limit = 2)[2])
            options.steps_per_unit > 0 || error("--steps-per-unit must be positive")
        elseif startswith(arg, "--repetitions=")
            options.repetitions = parse(Int, split(arg, "="; limit = 2)[2])
            options.repetitions > 0 || error("--repetitions must be positive")
        elseif startswith(arg, "--max-attempts=")
            options.max_attempts = parse(Int, split(arg, "="; limit = 2)[2])
            options.max_attempts > 0 || error("--max-attempts must be positive")
        else
            error("unknown option $(repr(arg)); use --help for usage")
        end
    end
    options.rerun && !options.run_experiment &&
        error("--rerun cannot be combined with --plot-only")
    for final_time in options.intervals
        raw_steps = final_time * options.steps_per_unit
        isapprox(raw_steps, round(raw_steps); rtol = 0, atol = 64eps(raw_steps)) ||
            error("T=$(final_time) does not give an integer N with the selected steps per unit")
    end
    return options
end

function experiment_configurations(options::Options)
    configurations = NamedTuple[]
    for dimension in options.dimensions
        push!(configurations, (
            experiment = "dimension",
            dimension = dimension,
            final_time = options.dimension_t,
            N = options.dimension_n,
        ))
    end
    for final_time in options.intervals
        push!(configurations, (
            experiment = "interval",
            dimension = options.interval_dimension,
            final_time = final_time,
            N = round(Int, final_time * options.steps_per_unit),
        ))
    end
    return configurations
end

empty_results() = DataFrame(
    implementation = String[],
    measurement = String[],
    experiment = String[],
    trial = Int[],
    dimension = Int[],
    final_time = Float64[],
    N = Int[],
    h = Float64[],
    runtime_s = Float64[],
    peak_rss_mib = Float64[],
)

function load_results()
    isfile(RESULTS_FILE) || return empty_results()
    data = DataFrame(CSV.File(RESULTS_FILE))
    required = [
        :implementation,
        :measurement,
        :experiment,
        :trial,
        :dimension,
        :final_time,
        :N,
        :h,
        :runtime_s,
        :peak_rss_mib,
    ]
    missing_columns = setdiff(required, propertynames(data))
    isempty(missing_columns) || error(
        "$(RESULTS_FILE) is missing columns: $(join(missing_columns, ", "))"
    )
    data.implementation = String.(data.implementation)
    data.measurement = String.(data.measurement)
    data.experiment = String.(data.experiment)
    data.trial = Int.(data.trial)
    data.dimension = Int.(data.dimension)
    data.final_time = Float64.(data.final_time)
    data.N = Int.(data.N)
    data.h = Float64.(data.h)
    data.runtime_s = Float64.(data.runtime_s)
    data.peak_rss_mib = Float64.(data.peak_rss_mib)
    return data[:, required]
end

function write_results(data::DataFrame)
    mkpath(OUTPUT_DIR)
    sort!(data, [
        :implementation,
        :measurement,
        :experiment,
        :dimension,
        :final_time,
        :N,
        :trial,
    ])
    temporary_file = RESULTS_FILE * ".tmp"
    CSV.write(temporary_file, data)
    mv(temporary_file, RESULTS_FILE; force = true)
end

function find_executable(environment_variable::String, fallback::String)
    requested = get(ENV, environment_variable, fallback)
    executable = if occursin('/', requested) || occursin('\\', requested)
        isfile(requested) ? abspath(requested) : nothing
    else
        Sys.which(requested)
    end
    if isnothing(executable) && fallback == "matlab" && Sys.isapple()
        applications = isdir("/Applications") ? readdir("/Applications"; join = true) : String[]
        matlab_apps = sort(filter(
            path -> startswith(basename(path), "MATLAB") && endswith(path, ".app"),
            applications,
        ))
        if !isempty(matlab_apps)
            app_executable = joinpath(last(matlab_apps), "bin", "matlab")
            executable = isfile(app_executable) ? app_executable : nothing
        end
    end
    isnothing(executable) && error(
        "cannot find $(fallback); set $(environment_variable) to its executable path"
    )
    return executable
end

function python_executable()
    if haskey(ENV, "SCALABILITY_PYTHON")
        return find_executable("SCALABILITY_PYTHON", "python3")
    end
    candidates = filter(!isnothing, Any[
        Sys.which("python3"),
        joinpath(homedir(), "miniconda3", "envs", "fde", "bin", "python"),
        joinpath(homedir(), "miniconda3", "bin", "python"),
    ])
    for candidate in unique(String.(candidates))
        isfile(candidate) || continue
        import_check = pipeline(
            ignorestatus(`$(candidate) -c "import pycaputo"`);
            stdout = devnull,
            stderr = devnull,
        )
        success(run(import_check)) && return candidate
    end
    error("cannot find Python with pycaputo; set SCALABILITY_PYTHON to its path")
end

function worker_command(configuration; baseline::Bool = false)
    if baseline
        arguments = ["--baseline"]
    else
        arguments = [
            string(configuration.dimension),
            @sprintf("%.17g", configuration.final_time),
            string(configuration.N),
        ]
    end
    return arguments
end

function worker_command(implementation::String, configuration; baseline::Bool = false)
    arguments = worker_command(configuration; baseline)
    if implementation == "Julia"
        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        script = joinpath(SCALABILITY_DIR, "fode_scalability.jl")
        return `$(julia) --project=$(PROJECT_DIR) --startup-file=no --threads=1 $(script) $(arguments)`
    elseif implementation == "Python"
        python = python_executable()
        script = joinpath(SCALABILITY_DIR, "fode_scalability.py")
        return `$(python) -u $(script) $(arguments)`
    elseif implementation == "MATLAB"
        matlab = find_executable("SCALABILITY_MATLAB", "matlab")
        matlab_dir = replace(SCALABILITY_DIR, "'" => "''")
        call = if baseline
            "fode_scalability('baseline')"
        else
            @sprintf(
                "fode_scalability(%d,%.17g,%d)",
                configuration.dimension,
                configuration.final_time,
                configuration.N,
            )
        end
        return `$(matlab) -singleCompThread -batch $("addpath('$(matlab_dir)'); $(call);")`
    end
    error("unsupported implementation $(repr(implementation))")
end

function parse_peak_rss(stderr_text::String)
    if Sys.isapple()
        matched = match(r"(?m)^\s*(\d+)\s+maximum resident set size\s*$", stderr_text)
        isnothing(matched) && error("could not parse macOS peak RSS from /usr/bin/time")
        return parse(Float64, matched.captures[1]) / 1024.0^2
    elseif Sys.islinux()
        matched = match(
            r"(?im)^\s*Maximum resident set size \(kbytes\):\s*(\d+)\s*$",
            stderr_text,
        )
        isnothing(matched) && error("could not parse Linux peak RSS from /usr/bin/time -v")
        return parse(Float64, matched.captures[1]) / 1024.0
    end
    error("OS-level peak RSS collection is supported on macOS and Linux")
end

function run_one(implementation::String, configuration; baseline::Bool = false)
    isfile("/usr/bin/time") || error("/usr/bin/time is required for RSS measurement")
    worker = worker_command(implementation, configuration; baseline)
    measured = if Sys.isapple()
        `/usr/bin/time -l $(worker)`
    elseif Sys.islinux()
        `/usr/bin/time -v $(worker)`
    else
        error("OS-level peak RSS collection is supported on macOS and Linux")
    end
    measured = addenv(
        measured,
        "OMP_NUM_THREADS" => "1",
        "OPENBLAS_NUM_THREADS" => "1",
        "MKL_NUM_THREADS" => "1",
    )

    stdout_buffer = IOBuffer()
    stderr_buffer = IOBuffer()
    process = run(
        pipeline(measured; stdout = stdout_buffer, stderr = stderr_buffer);
        wait = false,
    )
    wait(process)
    stdout_text = String(take!(stdout_buffer))
    stderr_text = String(take!(stderr_buffer))
    if !success(process)
        label = baseline ? "baseline" :
            "$(configuration.experiment), d=$(configuration.dimension), " *
            "T=$(configuration.final_time), N=$(configuration.N)"
        error(
            "$(implementation), $(label) failed.\n" *
            "stdout:\n$(stdout_text)\nstderr:\n$(stderr_text)"
        )
    end

    matches = collect(eachmatch(
        r"SCALABILITY_RESULT,(\d+),([^,\s]+),(\d+),([^,\s]+),([^,\s]+)",
        stdout_text,
    ))
    length(matches) == 1 || error(
        "expected one result line from $(implementation); received $(length(matches))"
    )
    captures = only(matches).captures
    returned_dimension = parse(Int, captures[1])
    returned_time = parse(Float64, captures[2])
    returned_n = parse(Int, captures[3])
    h = parse(Float64, captures[4])
    runtime_seconds = parse(Float64, captures[5])

    if baseline
        returned_dimension == 0 && returned_time == 0 && returned_n == 0 &&
            h == 0 && runtime_seconds == 0 || error("invalid baseline result")
    else
        returned_dimension == configuration.dimension || error("worker returned wrong d")
        isapprox(returned_time, configuration.final_time; rtol = 32eps()) ||
            error("worker returned wrong T")
        returned_n == configuration.N || error("worker returned wrong N")
        expected_h = configuration.final_time / configuration.N
        isapprox(h, expected_h; rtol = 32eps()) || error("worker returned wrong h")
        runtime_seconds > 0 || error("worker returned non-positive runtime")
    end

    peak_rss_mib = parse_peak_rss(stderr_text)
    peak_rss_mib > 0 || error("/usr/bin/time returned non-positive peak RSS")
    return (; h, runtime_seconds, peak_rss_mib)
end

function measure_with_retries(
    implementation::String,
    configuration,
    max_attempts::Int;
    baseline::Bool = false,
)
    for attempt in 1:max_attempts
        try
            return run_one(implementation, configuration; baseline)
        catch exception
            attempt == max_attempts && rethrow()
            label = baseline ? "baseline" :
                "$(configuration.experiment), d=$(configuration.dimension), " *
                "T=$(configuration.final_time), N=$(configuration.N)"
            @warn "Worker failed; retrying in a fresh process" implementation label attempt max_attempts exception = sprint(showerror, exception)
        end
    end
    error("unreachable retry state")
end

function matches_configuration(row, configuration)
    return row.experiment == configuration.experiment &&
           row.dimension == configuration.dimension &&
           row.final_time == configuration.final_time &&
           row.N == configuration.N
end

function selected_row(row, options::Options, configurations)
    row.implementation in options.implementations || return false
    row.measurement == "baseline" && return true
    return any(configuration -> matches_configuration(row, configuration), configurations)
end

function run_experiment!(data::DataFrame, options::Options)
    configurations = experiment_configurations(options)
    if options.rerun && nrow(data) > 0
        keep = [!selected_row(row, options, configurations) for row in eachrow(data)]
        data = data[keep, :]
        write_results(data)
    end

    total = length(options.implementations) * (length(configurations) + 1) *
            options.repetitions
    completed = 0
    baseline_configuration = (
        experiment = "baseline",
        dimension = 0,
        final_time = 0.0,
        N = 0,
    )

    for implementation in options.implementations
        all_configurations = [(; baseline_configuration..., baseline = true)]
        append!(all_configurations, [(; configuration..., baseline = false)
                                     for configuration in configurations])
        for configuration in all_configurations
            for trial in 1:options.repetitions
                completed += 1
                already_measured = any(
                    (data.implementation .== implementation) .&
                    (data.measurement .== (configuration.baseline ? "baseline" : "solve")) .&
                    (data.experiment .== configuration.experiment) .&
                    (data.trial .== trial) .&
                    (data.dimension .== configuration.dimension) .&
                    (data.final_time .== configuration.final_time) .&
                    (data.N .== configuration.N)
                )
                if already_measured
                    @info "Skipping existing measurement" implementation experiment = configuration.experiment dimension = configuration.dimension final_time = configuration.final_time N = configuration.N trial completed total
                    continue
                end

                @info "Measuring scalability" implementation experiment = configuration.experiment dimension = configuration.dimension final_time = configuration.final_time N = configuration.N trial completed total
                result = measure_with_retries(
                    implementation,
                    configuration,
                    options.max_attempts;
                    baseline = configuration.baseline,
                )
                push!(data, (
                    implementation = implementation,
                    measurement = configuration.baseline ? "baseline" : "solve",
                    experiment = configuration.experiment,
                    trial = trial,
                    dimension = configuration.dimension,
                    final_time = configuration.final_time,
                    N = configuration.N,
                    h = result.h,
                    runtime_s = result.runtime_seconds,
                    peak_rss_mib = result.peak_rss_mib,
                ))
                write_results(data)
                if configuration.baseline
                    @printf(
                        "%s baseline trial %d: peak RSS %.3f MiB\n",
                        implementation,
                        trial,
                        result.peak_rss_mib,
                    )
                else
                    @printf(
                        "%s %s d=%d T=%g N=%d trial %d: runtime %.6g s, peak RSS %.3f MiB\n",
                        implementation,
                        configuration.experiment,
                        configuration.dimension,
                        configuration.final_time,
                        configuration.N,
                        trial,
                        result.runtime_seconds,
                        result.peak_rss_mib,
                    )
                end
            end
        end
    end
    return data
end

function active_data(data::DataFrame, options::Options)
    configurations = experiment_configurations(options)
    keep = [selected_row(row, options, configurations) for row in eachrow(data)]
    return data[keep, :]
end

function summarize(data::DataFrame, implementation::String, experiment::String)
    rows = data[
        (data.implementation .== implementation) .&
        (data.measurement .== "solve") .&
        (data.experiment .== experiment),
        :,
    ]
    x = experiment == "dimension" ? Float64.(sort(unique(rows.dimension))) :
        sort(unique(rows.final_time))
    runtime_median = Float64[]
    runtime_q25 = Float64[]
    runtime_q75 = Float64[]
    rss_median = Float64[]
    rss_q25 = Float64[]
    rss_q75 = Float64[]
    for value in x
        selected = experiment == "dimension" ? rows.dimension .== Int(value) :
            rows.final_time .== value
        runtime_values = rows.runtime_s[selected]
        rss_values = rows.peak_rss_mib[selected]
        push!(runtime_median, median(runtime_values))
        push!(runtime_q25, quantile(runtime_values, 0.25))
        push!(runtime_q75, quantile(runtime_values, 0.75))
        push!(rss_median, median(rss_values))
        push!(rss_q25, quantile(rss_values, 0.25))
        push!(rss_q75, quantile(rss_values, 0.75))
    end
    return (; x, runtime_median, runtime_q25, runtime_q75,
            rss_median, rss_q25, rss_q75)
end

function draw_series!(runtime_axis, memory_axis, summary, color, marker)
    band!(
        runtime_axis,
        summary.x,
        summary.runtime_q25,
        summary.runtime_q75;
        color = (color, 0.13),
    )
    line = lines!(
        runtime_axis,
        summary.x,
        summary.runtime_median;
        color,
        linewidth = 3.0,
    )
    points = scatter!(
        runtime_axis,
        summary.x,
        summary.runtime_median;
        color,
        marker,
        markersize = 13,
        strokecolor = :white,
        strokewidth = 1,
    )
    band!(
        memory_axis,
        summary.x,
        summary.rss_q25,
        summary.rss_q75;
        color = (color, 0.13),
    )
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
        markersize = 13,
        strokecolor = :white,
        strokewidth = 1,
    )
    return line, points
end

function draw_baseline!(axis, data, implementation::String, x, color)
    rows = data[
        (data.implementation .== implementation) .&
        (data.measurement .== "baseline"),
        :,
    ]
    nrow(rows) == 0 && return false
    baseline_median = median(rows.peak_rss_mib)
    baseline_q25 = quantile(rows.peak_rss_mib, 0.25)
    baseline_q75 = quantile(rows.peak_rss_mib, 0.75)
    baseline_x = [first(x), last(x)]
    band!(
        axis,
        baseline_x,
        fill(baseline_q25, 2),
        fill(baseline_q75, 2);
        color = (color, 0.06),
    )
    hlines!(
        axis,
        [baseline_median];
        color = (color, 0.75),
        linestyle = :dash,
        linewidth = 2.0,
    )
    return true
end

function make_figure(data::DataFrame, options::Options)
    data = active_data(data, options)
    solve_data = data[data.measurement .== "solve", :]
    nrow(solve_data) > 0 || error("no active scalability measurements to plot")
    all(solve_data.runtime_s .> 0) || error("runtime must be positive")
    all(solve_data.peak_rss_mib .> 0) || error("peak RSS must be positive")

    dimension_ticks = Float64.(options.dimensions)
    interval_ticks = options.intervals
    grid_color = (:gray35, 0.16)
    figure = Figure(size = (1380, 900), fontsize = 17)

    common_axis = (
        xscale = log10,
        xgridvisible = false,
        ygridvisible = true,
        ygridcolor = grid_color,
        ygridwidth = 0.8,
        spinewidth = 1.2,
    )
    dimension_runtime_axis = Axis(
        figure[1, 1];
        title = "(a) Runtime vs system dimension",
        xlabel = "System dimension, d",
        ylabel = "Runtime (s)",
        yscale = log10,
        xticks = (dimension_ticks, string.(options.dimensions)),
        common_axis...,
    )
    dimension_memory_axis = Axis(
        figure[1, 2];
        title = "(b) Memory vs system dimension",
        xlabel = "System dimension, d",
        ylabel = "Peak RSS (MiB)",
        xticks = (dimension_ticks, string.(options.dimensions)),
        common_axis...,
    )
    interval_runtime_axis = Axis(
        figure[2, 1];
        title = "(c) Runtime vs integration interval",
        xlabel = "Final time, T",
        ylabel = "Runtime (s)",
        yscale = log10,
        xticks = (interval_ticks, [@sprintf("%g", value) for value in interval_ticks]),
        common_axis...,
    )
    interval_memory_axis = Axis(
        figure[2, 2];
        title = "(d) Memory vs integration interval",
        xlabel = "Final time, T",
        ylabel = "Peak RSS (MiB)",
        xticks = (interval_ticks, [@sprintf("%g", value) for value in interval_ticks]),
        common_axis...,
    )

    legend_entries = Any[]
    legend_labels = String[]
    has_baseline = false
    for implementation in options.implementations
        color = COLORS[implementation]
        marker = MARKERS[implementation]
        dimension_summary = summarize(data, implementation, "dimension")
        interval_summary = summarize(data, implementation, "interval")

        legend_line = nothing
        legend_points = nothing
        if !isempty(dimension_summary.x)
            legend_line, legend_points = draw_series!(
                dimension_runtime_axis,
                dimension_memory_axis,
                dimension_summary,
                color,
                marker,
            )
            has_baseline |= draw_baseline!(
                dimension_memory_axis,
                data,
                implementation,
                dimension_summary.x,
                color,
            )
        end
        if !isempty(interval_summary.x)
            interval_line, interval_points = draw_series!(
                interval_runtime_axis,
                interval_memory_axis,
                interval_summary,
                color,
                marker,
            )
            isnothing(legend_line) && (legend_line = interval_line)
            isnothing(legend_points) && (legend_points = interval_points)
            has_baseline |= draw_baseline!(
                interval_memory_axis,
                data,
                implementation,
                interval_summary.x,
                color,
            )
        end
        if !isnothing(legend_line)
            push!(legend_entries, [legend_line, legend_points])
            push!(legend_labels, IMPLEMENTATION_LABELS[implementation])
        end
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
        figure[0, 1:2],
        legend_entries,
        legend_labels;
        orientation = :horizontal,
        framevisible = false,
        tellheight = true,
        tellwidth = false,
        patchsize = (38.0f0, 18.0f0),
        labelsize = 15,
    )
    Label(
        figure[3, 1:2],
        "Dimension sweep: T=$(options.dimension_t), N=$(options.dimension_n). " *
        "Interval sweep: d=$(options.interval_dimension), h=1/$(options.steps_per_unit) " *
        "and N=$(options.steps_per_unit)T. Lines/points: median; bands: IQR; " *
        "dashed lines: initialized-process RSS.";
        fontsize = 13,
        color = :gray35,
        tellwidth = false,
    )

    ylims!(dimension_memory_axis; low = 0)
    ylims!(interval_memory_axis; low = 0)
    colgap!(figure.layout, 46)
    rowgap!(figure.layout, 24)
    mkpath(OUTPUT_DIR)
    save(FIGURE_PDF, figure)
    save(FIGURE_SVG, figure)
    save(FIGURE_PNG, figure; px_per_unit = 2)
    return figure
end

function main(args)
    options = parse_options(args)
    mkpath(OUTPUT_DIR)
    data = load_results()
    if options.run_experiment
        data = run_experiment!(data, options)
    end
    make_figure(data, options)
    println("Results: $(RESULTS_FILE)")
    println("Figures: $(FIGURE_PDF), $(FIGURE_SVG), $(FIGURE_PNG)")
end

main(ARGS)
