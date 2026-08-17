#!/usr/bin/env julia

"""
Run one Julia PECE solve for the dimension/interval scalability benchmark.

Usage:
    julia --project=. scalability/fode_scalability.jl D T N
    julia --project=. scalability/fode_scalability.jl --baseline

The benchmark is the commensurate system

    D_t^0.8 u(t) = A_d u(t),  u(0) = [1.0, ..., 0.1],

where A_d has -1.5 on the diagonal and 0.25 on its first off-diagonals.
The final stdout line is parsed by `scalability_figure.jl`.
"""

using FractionalDiffEq
using LinearAlgebra
using Printf

const ALPHA = 0.8
const T0 = 0.0

function system_matrix(dimension::Int)
    matrix = -1.5 * Matrix{Float64}(I, dimension, dimension)
    for index in 1:(dimension - 1)
        matrix[index, index + 1] = 0.25
        matrix[index + 1, index] = 0.25
    end
    return matrix
end

initial_state(dimension::Int) = collect(range(1.0, 0.1; length = dimension))

function rhs!(du, u, matrix, _t)
    mul!(du, matrix, u)
    return nothing
end

function solve_pece(dimension::Int, final_time::Float64, number_of_steps::Int)
    h = (final_time - T0) / number_of_steps
    orders = fill(ALPHA, dimension)
    problem = FODEProblem(
        ODEFunction{true}(rhs!),
        orders,
        initial_state(dimension),
        (T0, final_time),
        system_matrix(dimension),
    )
    return solve(problem, PECE(); dt = h), h
end

function parse_positive_int(value::String, name::String)
    parsed = try
        parse(Int, value)
    catch
        error("$(name) must be a positive integer; received $(repr(value))")
    end
    parsed > 0 || error("$(name) must be positive")
    return parsed
end

function parse_positive_float(value::String, name::String)
    parsed = try
        parse(Float64, value)
    catch
        error("$(name) must be positive; received $(repr(value))")
    end
    isfinite(parsed) && parsed > T0 || error("$(name) must be finite and positive")
    return parsed
end

function main(args)
    if args == ["--baseline"]
        @printf("SCALABILITY_RESULT,0,0,0,0,0\n")
        return
    end
    length(args) == 3 || error("usage: fode_scalability.jl D T N | --baseline")
    dimension = parse_positive_int(args[1], "D")
    final_time = parse_positive_float(args[2], "T")
    number_of_steps = parse_positive_int(args[3], "N")
    h = (final_time - T0) / number_of_steps

    # Compile the same solver path using the requested h over at most 64 steps.
    warmup_steps = min(number_of_steps, 64)
    warmup_final_time = T0 + warmup_steps * h
    warmup_solution, _ = solve_pece(dimension, warmup_final_time, warmup_steps)
    isempty(warmup_solution.u) && error("Julia PECE warm-up returned no values")
    warmup_solution = nothing
    GC.gc(true)

    solution = nothing
    runtime_seconds = @elapsed begin
        solution, returned_h = solve_pece(dimension, final_time, number_of_steps)
        returned_h == h || error("internal step-size mismatch")
    end

    length(solution.t) == number_of_steps + 1 ||
        error("Julia PECE returned an unexpected number of time points")
    isapprox(last(solution.t), final_time; rtol = 0, atol = 64eps(final_time)) ||
        error("Julia PECE did not reach the requested final time")
    final_state = last(solution.u)
    length(final_state) == dimension || error("Julia PECE returned an invalid state size")
    all(isfinite, final_state) || error("Julia PECE returned a non-finite solution")

    @printf(
        "SCALABILITY_RESULT,%d,%.17g,%d,%.17g,%.17g\n",
        dimension,
        final_time,
        number_of_steps,
        h,
        runtime_seconds,
    )
end

main(ARGS)
