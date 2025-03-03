using FractionalDiffEq
using SpecialFunctions
using FdeSolver
using DataFrames
using BenchmarkTools
using CSV
using LinearAlgebra

alpha = 0.5;
f(u, h, p, t) = 2/gamma(3-alpha)*t^(2-alpha) - 1/gamma(2-alpha)*t^(1-alpha) + 2*tau*t - tau^2 - tau - u + h
h(p, t) = 0.0
tau = 0.5;
u0 = 0.0;
tspan = (0.0, 5.0)
analytical(u, h, p, t) = [t^2-t]
fun = DDEFunction(f)#, analytic = analytical)
prob = FDDEProblem(f, alpha, u0, h, constant_lags = [tau], tspan)
sol = solve(prob, DelayPECE(), dt=0.01)

# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[];
E3 = Float64[];T3 = Float64[];E4 = Float64[];T4 = Float64[];
E5 = Float64[];T5 = Float64[];

for n in range(3, 7)
    println("n: $n")# to print out the current step of runing
    dt = 2.0^-n #stepsize of computting
        #computing the time
    t1= @benchmark solve($(prob), $(DelayPECE()), dt = $(dt)) seconds=1
    t2= @benchmark solve($(prob), $(DelayPIEX()), dt = $(dt)); seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)

    #computing the error
    sol1 = solve(prob, DelayPECE(), dt = dt);
    sol2 = solve(prob, DelayPIEX(), dt = dt);

    ery1=norm(sol1.u - analytical.(0, 0, 0, sol1.t))
    ery2=norm(sol2.u - analytical.(0, 0, 0, sol2.t))

    push!(E1, ery1)
    push!(E2, ery2)
end

#save data
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPECE.csv",  DataFrame(time = eval(Symbol("T1")), error = eval(Symbol("E1"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIEX.csv",  DataFrame(time = eval(Symbol("T3")), error = eval(Symbol("E3"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPITrap.csv",  DataFrame(time = eval(Symbol("T4")), error = eval(Symbol("E4"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIRect.csv",  DataFrame(time = eval(Symbol("T5")), error = eval(Symbol("E5"))))

