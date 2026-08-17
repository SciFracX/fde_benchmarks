#!/usr/bin/env julia

"""
Run and plot the cross-language FODE runtime/peak-RSS experiment.

Examples (run from the project root):
    julia --project=. memory_fode/memory_figure.jl
    julia --project=. memory_fode/memory_figure.jl --ns=320,640,1280 --repetitions=5
    julia --project=. memory_fode/memory_figure.jl --plot-only
    julia --project=. memory_fode/memory_figure.jl --rerun

Environment overrides:
    MEMORY_PYTHON=/path/to/python-with-pycaputo
    MEMORY_MATLAB=/path/to/matlab
    MEMORY_OUTPUT_DIR=/path/for/csv-and-figures

Each (implementation, N) pair runs in a fresh process.  The worker reports
solver-only wall time, while `/usr/bin/time` reports peak RSS for the whole
process, including the language runtime and loaded libraries.
"""

using CairoMakie
using CSV
using DataFrames
using Printf
using Statistics

const MEMORY_DIR = @__DIR__
const PROJECT_DIR = normpath(joinpath(MEMORY_DIR, ".."))
const OUTPUT_DIR = abspath(get(ENV, "MEMORY_OUTPUT_DIR", MEMORY_DIR))
const RESULTS_FILE = joinpath(OUTPUT_DIR, "memory_results.csv")
const FIGURE_PDF = joinpath(OUTPUT_DIR, "memory_consumption.pdf")
const FIGURE_SVG = joinpath(OUTPUT_DIR, "memory_consumption.svg")
const FIGURE_PNG = joinpath(OUTPUT_DIR, "memory_consumption.png")
const DEFAULT_NS = [320, 640, 1280, 2560, 5120, 10240, 20480, 40960]
const DEFAULT_REPETITIONS = 5
const DEFAULT_MAX_ATTEMPTS = 3
const IMPLEMENTATIONS = ["Julia", "MATLAB", "Python"]
const PROBLEM_ID = "coupled_linear_d10_alpha0.8_t20"
const LEGACY_PROBLEM_ID = "scalar_linear_d1_alpha0.8_t5"
const INTERVAL_LENGTH = 20.0

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
    ns::Vector{Int} = copy(DEFAULT_NS)
    implementations::Vector{String} = copy(IMPLEMENTATIONS)
    repetitions::Int = DEFAULT_REPETITIONS
    max_attempts::Int = DEFAULT_MAX_ATTEMPTS
end

function usage()
    println("""
    usage: julia --project=. memory_fode/memory_figure.jl [options]

      --ns=N1,N2,...                 time-step counts (default: $(join(DEFAULT_NS, ',')))
      --languages=Julia,MATLAB,Python
                                      implementations to run
      --repetitions=R                 independent processes per point (default: $(DEFAULT_REPETITIONS))
      --max-attempts=A                attempts after a worker startup failure (default: $(DEFAULT_MAX_ATTEMPTS))
      --plot-only                     plot the existing memory_results.csv
      --rerun                         replace selected existing measurements
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
        elseif startswith(arg, "--ns=") || startswith(arg, "--n=")
            raw = split(arg, "="; limit = 2)[2]
            options.ns = parse.(Int, split(raw, ','))
            isempty(options.ns) && error("--ns must contain at least one integer")
            all(>(0), options.ns) || error("all N values must be positive")
            length(unique(options.ns)) == length(options.ns) ||
                error("--ns contains duplicate values")
            sort!(options.ns)
        elseif startswith(arg, "--languages=")
            raw = split(arg, "="; limit = 2)[2]
            names = filter(!isempty, strip.(split(raw, ',')))
            isempty(names) && error("--languages must contain at least one name")
            options.implementations = unique(canonical_implementation.(names))
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
    return options
end

empty_results() = DataFrame(
    problem_id = String[],
    implementation = String[],
    measurement = String[],
    trial = Int[],
    N = Int[],
    h = Float64[],
    runtime_s = Float64[],
    peak_rss_mib = Float64[],
)

function load_results()
    isfile(RESULTS_FILE) || return empty_results()
    data = DataFrame(CSV.File(RESULTS_FILE))
    legacy_required = [:implementation, :N, :h, :runtime_s, :peak_rss_mib]
    missing_columns = setdiff(legacy_required, propertynames(data))
    isempty(missing_columns) || error(
        "$(RESULTS_FILE) is missing columns: $(join(missing_columns, ", "))"
    )
    if :measurement ∉ propertynames(data)
        insertcols!(data, 2, :measurement => fill("solve", nrow(data)))
    end
    if :trial ∉ propertynames(data)
        trials = zeros(Int, nrow(data))
        counts = Dict{Tuple{String, String, Int}, Int}()
        for row_index in 1:nrow(data)
            key = (
                String(data.implementation[row_index]),
                String(data.measurement[row_index]),
                Int(data.N[row_index]),
            )
            trials[row_index] = get(counts, key, 0) + 1
            counts[key] = trials[row_index]
        end
        insertcols!(data, 3, :trial => trials)
    end
    if :problem_id ∉ propertynames(data)
        insertcols!(data, 1, :problem_id => fill(LEGACY_PROBLEM_ID, nrow(data)))
    end

    required = [
        :problem_id,
        :implementation,
        :measurement,
        :trial,
        :N,
        :h,
        :runtime_s,
        :peak_rss_mib,
    ]
    data.problem_id = String.(data.problem_id)
    data.implementation = String.(data.implementation)
    data.measurement = String.(data.measurement)
    data.trial = Int.(data.trial)
    data.N = Int.(data.N)
    data.h = Float64.(data.h)
    data.runtime_s = Float64.(data.runtime_s)
    data.peak_rss_mib = Float64.(data.peak_rss_mib)
    return data[:, required]
end

function write_results(data::DataFrame)
    mkpath(OUTPUT_DIR)
    sort!(data, [:problem_id, :implementation, :measurement, :N, :trial])
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
    if haskey(ENV, "MEMORY_PYTHON")
        return find_executable("MEMORY_PYTHON", "python3")
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
    error("cannot find a Python installation with pycaputo; set MEMORY_PYTHON to its path")
end

function worker_command(
    implementation::String,
    number_of_steps::Int;
    baseline::Bool = false,
)
    worker_argument = baseline ? "--baseline" : string(number_of_steps)
    if implementation == "Julia"
        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        return `$(julia) --project=$(PROJECT_DIR) --startup-file=no --threads=1 $(joinpath(MEMORY_DIR, "fode_memory.jl")) $(worker_argument)`
    elseif implementation == "Python"
        python = python_executable()
        return `$(python) -u $(joinpath(MEMORY_DIR, "fode_memory.py")) $(worker_argument)`
    elseif implementation == "MATLAB"
        matlab = find_executable("MEMORY_MATLAB", "matlab")
        matlab_dir = replace(MEMORY_DIR, "'" => "''")
        call = baseline ? "fode_memory('baseline')" : "fode_memory($(number_of_steps))"
        expression = "addpath('$(matlab_dir)'); $(call);"
        return `$(matlab) -singleCompThread -batch $(expression)`
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
    error("OS-level peak RSS collection is currently supported on macOS and Linux")
end

function run_one(
    implementation::String,
    number_of_steps::Int;
    baseline::Bool = false,
)
    isfile("/usr/bin/time") || error("/usr/bin/time is required for peak RSS measurement")
    worker = worker_command(implementation, number_of_steps; baseline)
    measured = if Sys.isapple()
        `/usr/bin/time -l $(worker)`
    elseif Sys.islinux()
        `/usr/bin/time -v $(worker)`
    else
        error("OS-level peak RSS collection is currently supported on macOS and Linux")
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
        configuration = baseline ? "baseline" : "N=$(number_of_steps)"
        error(
            "$(implementation), $(configuration) failed.\n" *
            "stdout:\n$(stdout_text)\nstderr:\n$(stderr_text)"
        )
    end

    result_matches = collect(eachmatch(
        r"FODE_MEMORY_RESULT,(\d+),([^,\s]+),([^,\s]+)",
        stdout_text,
    ))
    length(result_matches) == 1 || error(
        "expected exactly one result line from $(implementation), N=$(number_of_steps); " *
        "received $(length(result_matches))"
    )
    captures = only(result_matches).captures
    returned_n = parse(Int, captures[1])
    h = parse(Float64, captures[2])
    runtime_seconds = parse(Float64, captures[3])
    expected_n = baseline ? 0 : number_of_steps
    returned_n == expected_n || error("worker returned N=$(returned_n), expected $(expected_n)")
    if baseline
        h == 0 && runtime_seconds == 0 ||
            error("baseline worker must return h=0 and runtime=0")
    else
        isapprox(h, INTERVAL_LENGTH / number_of_steps; rtol = 32eps()) ||
            error("worker returned an inconsistent step size h=$(h)")
        runtime_seconds > 0 || error("worker returned a non-positive runtime")
    end

    peak_rss_mib = parse_peak_rss(stderr_text)
    peak_rss_mib > 0 || error("/usr/bin/time returned a non-positive peak RSS")
    return (; h, runtime_seconds, peak_rss_mib)
end

function measure_with_retries(
    implementation::String,
    number_of_steps::Int,
    max_attempts::Int;
    baseline::Bool = false,
)
    for attempt in 1:max_attempts
        try
            return run_one(implementation, number_of_steps; baseline)
        catch exception
            attempt == max_attempts && rethrow()
            configuration = baseline ? "baseline" : "N=$(number_of_steps)"
            @warn "Worker failed; retrying in a fresh process" implementation configuration attempt max_attempts exception = sprint(showerror, exception)
        end
    end
    error("unreachable retry state")
end

function selected_pair(row, options::Options)
    selected_problem = row.problem_id == PROBLEM_ID
    selected_implementation = row.implementation in options.implementations
    selected_measurement = row.measurement == "baseline" || row.N in options.ns
    return selected_problem && selected_implementation && selected_measurement
end

function run_experiment!(data::DataFrame, options::Options)
    if options.rerun && nrow(data) > 0
        keep = [!selected_pair(row, options) for row in eachrow(data)]
        data = data[keep, :]
        write_results(data)
    end

    total = length(options.implementations) * (length(options.ns) + 1) * options.repetitions
    completed = 0
    for implementation in options.implementations
        configurations = [(; measurement = "baseline", N = 0, baseline = true)]
        append!(configurations, [
            (; measurement = "solve", N = number_of_steps, baseline = false)
            for number_of_steps in options.ns
        ])

        for configuration in configurations
            for trial in 1:options.repetitions
                completed += 1
                already_measured = any(
                    (data.problem_id .== PROBLEM_ID) .&
                    (data.implementation .== implementation) .&
                    (data.measurement .== configuration.measurement) .&
                    (data.N .== configuration.N) .&
                    (data.trial .== trial)
                )
                if already_measured
                    @info "Skipping existing measurement" implementation measurement = configuration.measurement number_of_steps = configuration.N trial completed total
                    continue
                end

                @info "Measuring FODE process" implementation measurement = configuration.measurement number_of_steps = configuration.N trial completed total
                result = measure_with_retries(
                    implementation,
                    configuration.N,
                    options.max_attempts;
                    baseline = configuration.baseline,
                )
                push!(data, (
                    problem_id = PROBLEM_ID,
                    implementation = implementation,
                    measurement = configuration.measurement,
                    trial = trial,
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
                        "%s N=%d trial %d: runtime %.6g s, peak RSS %.3f MiB\n",
                        implementation,
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

function summarize_solver(data::DataFrame, implementation::String)
    rows = data[
        (data.implementation .== implementation) .&
        (data.measurement .== "solve"),
        :,
    ]
    number_of_steps = sort(unique(rows.N))
    runtime_median = Float64[]
    runtime_q25 = Float64[]
    runtime_q75 = Float64[]
    rss_median = Float64[]
    rss_q25 = Float64[]
    rss_q75 = Float64[]

    for n in number_of_steps
        selected = rows.N .== n
        runtime_values = rows.runtime_s[selected]
        rss_values = rows.peak_rss_mib[selected]
        push!(runtime_median, median(runtime_values))
        push!(runtime_q25, quantile(runtime_values, 0.25))
        push!(runtime_q75, quantile(runtime_values, 0.75))
        push!(rss_median, median(rss_values))
        push!(rss_q25, quantile(rss_values, 0.25))
        push!(rss_q75, quantile(rss_values, 0.75))
    end

    return (;
        number_of_steps,
        runtime_median,
        runtime_q25,
        runtime_q75,
        rss_median,
        rss_q25,
        rss_q75,
    )
end

function empirical_exponent(number_of_steps, runtime_median)
    count = min(length(number_of_steps), 4)
    count >= 2 || return NaN
    indices = (length(number_of_steps) - count + 1):length(number_of_steps)
    x = log.(Float64.(number_of_steps[indices]))
    y = log.(runtime_median[indices])
    centered_x = x .- mean(x)
    denominator = sum(abs2, centered_x)
    denominator > 0 || return NaN
    return sum(centered_x .* (y .- mean(y))) / denominator
end

function format_n_tick(number_of_steps::Int)
    number_of_steps < 1000 && return string(number_of_steps)
    number_of_steps < 10000 && return @sprintf("%.2fk", number_of_steps / 1000)
    return @sprintf("%.1fk", number_of_steps / 1000)
end

function choose_ticks(number_of_steps)
    length(number_of_steps) <= 8 && return number_of_steps
    indices = unique(round.(Int, range(1, length(number_of_steps); length = 8)))
    return number_of_steps[indices]
end

function make_figure(data::DataFrame)
    data = data[data.problem_id .== PROBLEM_ID, :]
    solve_data = data[data.measurement .== "solve", :]
    nrow(solve_data) > 0 || error(
        "no measurements are available for $(PROBLEM_ID); run the experiment first"
    )
    all(solve_data.N .> 0) || error("solver N must be positive for logarithmic axes")
    all(solve_data.runtime_s .> 0) ||
        error("solver runtime must be positive for a logarithmic axis")

    all_number_of_steps = sort(unique(solve_data.N))
    tick_positions = choose_ticks(all_number_of_steps)
    tick_labels = format_n_tick.(tick_positions)
    grid_color = (:gray35, 0.16)

    figure = Figure(size = (1320, 540), fontsize = 18)
    runtime_axis = Axis(
        figure[1, 1];
        title = "(a) Runtime scaling",
        xlabel = "Number of time steps, N",
        ylabel = "Runtime (s)",
        xscale = log10,
        yscale = log10,
        xticks = (tick_positions, tick_labels),
        xgridvisible = false,
        ygridvisible = true,
        ygridcolor = grid_color,
        ygridwidth = 0.8,
        spinewidth = 1.2,
    )
    memory_axis = Axis(
        figure[1, 2];
        title = "(b) Full-process memory footprint",
        xlabel = "Number of time steps, N",
        ylabel = "Peak RSS (MiB)",
        xscale = log10,
        xticks = (tick_positions, tick_labels),
        xgridvisible = false,
        ygridvisible = true,
        ygridcolor = grid_color,
        ygridwidth = 0.8,
        spinewidth = 1.2,
    )

    legend_entries = Any[]
    legend_labels = String[]
    has_baseline = false
    for implementation in IMPLEMENTATIONS
        summary = summarize_solver(data, implementation)
        isempty(summary.number_of_steps) && continue
        color = COLORS[implementation]
        marker = MARKERS[implementation]
        show_iqr_band = implementation != "MATLAB"

        if show_iqr_band
            band!(
                runtime_axis,
                summary.number_of_steps,
                summary.runtime_q25,
                summary.runtime_q75;
                color = (color, 0.15),
            )
        end
        runtime_line = lines!(
            runtime_axis,
            summary.number_of_steps,
            summary.runtime_median;
            color = color,
            linewidth = 3.0,
        )
        runtime_points = scatter!(
            runtime_axis,
            summary.number_of_steps,
            summary.runtime_median;
            color = color,
            marker = marker,
            markersize = 14,
            strokecolor = :white,
            strokewidth = 1,
        )

        if show_iqr_band
            band!(
                memory_axis,
                summary.number_of_steps,
                summary.rss_q25,
                summary.rss_q75;
                color = (color, 0.15),
            )
        end
        lines!(
            memory_axis,
            summary.number_of_steps,
            summary.rss_median;
            color = color,
            linewidth = 3.0,
        )
        scatter!(
            memory_axis,
            summary.number_of_steps,
            summary.rss_median;
            color = color,
            marker = marker,
            markersize = 14,
            strokecolor = :white,
            strokewidth = 1,
        )

        baseline_rows = data[
            (data.implementation .== implementation) .&
            (data.measurement .== "baseline"),
            :,
        ]
        if nrow(baseline_rows) > 0
            has_baseline = true
            baseline_median = median(baseline_rows.peak_rss_mib)
            baseline_q25 = quantile(baseline_rows.peak_rss_mib, 0.25)
            baseline_q75 = quantile(baseline_rows.peak_rss_mib, 0.75)
            baseline_x = [first(summary.number_of_steps), last(summary.number_of_steps)]
            if show_iqr_band
                band!(
                    memory_axis,
                    baseline_x,
                    fill(baseline_q25, 2),
                    fill(baseline_q75, 2);
                    color = (color, 0.07),
                )
            end
            hlines!(
                memory_axis,
                [baseline_median];
                color = (color, 0.75),
                linestyle = :dash,
                linewidth = 2.0,
            )
        end

        exponent = empirical_exponent(
            summary.number_of_steps,
            summary.runtime_median,
        )
        exponent_label = isfinite(exponent) ? @sprintf("p = %.2f", exponent) : "p = n/a"
        push!(legend_entries, [runtime_line, runtime_points])
        push!(
            legend_labels,
            "$(IMPLEMENTATION_LABELS[implementation])  ($(exponent_label))",
        )
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
        figure[2, 1:2],
        "Lines and points: median; shaded bands: interquartile range " *
        "(Julia and Python only). " *
        "The exponent p is fitted from the largest four N values.";
        fontsize = 13,
        color = :gray35,
        tellwidth = false,
    )
    colgap!(figure.layout, 42)
    rowgap!(figure.layout, 8)
    ylims!(memory_axis; low = 0)
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
    make_figure(data)
    println("Results: $(RESULTS_FILE)")
    println("Figures: $(FIGURE_PDF), $(FIGURE_SVG), $(FIGURE_PNG)")
end

main(ARGS)
