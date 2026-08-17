#!/usr/bin/env julia

"""
Run one Julia PECE solve for the cross-language FODE memory benchmark.

The benchmark is the 10-state system D_t^0.8 u = A*u on [0, 20], where A
has -1.5 on the diagonal and 0.25 on its first off-diagonals.  The initial
state is [1.0, 0.9, ..., 0.1].

Usage:
    julia --project=. memory_fode/fode_memory.jl N
    julia --project=. memory_fode/fode_memory.jl --baseline

The last stdout line is machine-readable by `memory_figure.jl`.  Runtime is
measured around the numerical solve only; peak RSS is deliberately measured
by the parent process so that all three languages use the same OS metric.
"""

using FractionalDiffEq
using LinearAlgebra
using Printf

const ALPHA = 0.8
const T0 = 0.0
const TFINAL = 20.0
const SYSTEM_SIZE = 10
const ORDERS = fill(ALPHA, SYSTEM_SIZE)
const U0 = collect(range(1.0, 0.1; length = SYSTEM_SIZE))
const SYSTEM_MATRIX = let matrix = -1.5 * Matrix{Float64}(I, SYSTEM_SIZE, SYSTEM_SIZE)
    for index in 1:(SYSTEM_SIZE - 1)
        matrix[index, index + 1] = 0.25
        matrix[index + 1, index] = 0.25
    end
    matrix
end

function rhs!(du, u, _p, _t)
    mul!(du, SYSTEM_MATRIX, u)
    return nothing
end

function solve_pece(number_of_steps::Int; final_time::Float64 = TFINAL)
    h = (final_time - T0) / number_of_steps
    fun = ODEFunction{true}(rhs!)
    problem = FODEProblem(fun, ORDERS, U0, (T0, final_time))
    return solve(problem, PECE(); dt = h), h
end

function main(args)
    length(args) == 1 || error("usage: fode_memory.jl N")
    if only(args) == "--baseline"
        # Package loading has already happened at this point.  Do not invoke
        # the numerical method: this captures the initialized Julia process
        # against which full solver-process RSS can be compared.
        @printf("FODE_MEMORY_RESULT,0,0,0\n")
        return
    end
    number_of_steps = try
        parse(Int, only(args))
    catch
        error("N must be a positive integer; received $(repr(only(args)))")
    end
    number_of_steps > 0 || error("N must be positive")

    # Compile the solver path before timing.  The warm-up is intentionally
    # small, while remaining inside this process for full-process RSS.
    h = (TFINAL - T0) / number_of_steps
    warmup_steps = min(number_of_steps, 64)
    warmup_final_time = T0 + warmup_steps * h
    warmup_solution, _ = solve_pece(warmup_steps; final_time = warmup_final_time)
    isempty(warmup_solution.u) && error("Julia PECE warm-up returned no values")
    warmup_solution = nothing
    GC.gc(true)

    solution = nothing
    runtime_seconds = @elapsed begin
        solution, returned_h = solve_pece(number_of_steps)
        returned_h == h || error("internal step-size mismatch")
    end

    isempty(solution.u) && error("Julia PECE returned no values")
    length(solution.t) == number_of_steps + 1 ||
        error("Julia PECE returned an unexpected number of time points")
    isapprox(last(solution.t), TFINAL; rtol = 0, atol = 64eps(TFINAL)) ||
        error("Julia PECE did not reach the requested final time")
    final_state = last(solution.u)
    final_values = final_state isa Number ? (final_state,) : final_state
    all(isfinite, final_values) || error("Julia PECE returned a non-finite solution")

    @printf("FODE_MEMORY_RESULT,%d,%.17g,%.17g\n",
            number_of_steps, h, runtime_seconds)
end

main(ARGS)
