using FractionalDiffEq
using SpecialFunctions
using FdeSolver
using DataFrames
using BenchmarkTools
using CSV
using LinearAlgebra

tspan = (0, 100)
rightfun(u, p, t) = 6*cos(t)
prob = MultiTermsFODEProblem([1, 1, 1, 4, 1, 4], [3, 2.5, 2, 1, 0.5, 0], rightfun, [1, 1, -1], tspan)
arightfun(u, p, t) = 6*cos(t)-(-2*t^2-gamma(3)/(2*gamma(2.5))*t^1.5+gamma(2)/(gamma(1.5))*t^0.5+7);
aprob = MultiTermsFODEProblem([1, 1, 1, 4, 1, 4], [3, 2.5, 2, 1, 0.5, 0], arightfun, [0, 0, 0], tspan)
ana(t) = [sqrt(2)*sin(t+pi/4)+t^2/2-t-1];

analytical(t) = [sqrt(2)*sin(t + pi/4)]

# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[];
E3 = Float64[];T3 = Float64[];E4 = Float64[];T4 = Float64[];
E5 = Float64[];T5 = Float64[];
h = Float64[]

for n in range(3, 7)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computing the time
    t1= @benchmark solve($(prob), $(MTPECE()), dt = $(h)) seconds=1
    #t2= @benchmark FDEsolver(F, $(tSpan), $(y0), $(β), $(par), JF = JF, h=$(h)) seconds=1
    t3= @benchmark solve($(prob), $(MTPIEX()), dt = $(h)); seconds=1
    t4 = @benchmark solve($(prob), $(MTPITrap()), dt = $(h)); seconds=1
    t5 = @benchmark solve($(prob), $(MTPIRect()), dt = $(h)); seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    #push!(T2, mean(t2).time / 10^9)
    push!(T3, mean(t3).time / 10^9)
    push!(T4, mean(t4).time / 10^9)
    push!(T5, mean(t5).time / 10^9)

    #computing the error
    sol1 = solve(prob, MTPECE(), dt = h);
    #_, y2 = FDEsolver(F, tSpan, y0, β, par, JF = JF, h=h)
    sol3 =  solve(prob, MTPIEX(), dt = h);
    sol4 =  solve(prob, MTPITrap(), dt = h);
    sol5 =  solve(prob, MTPIRect(), dt = h);

    ery1=norm(sol1.u - analytical.(sol1.t))
    #ery2=sol2.errors[:final]
    ery3=norm(sol3.u - analytical.(sol3.t))
    ery4=norm(sol4.u - analytical.(sol4.t))
    ery5=norm(sol5.u - analytical.(sol5.t))

    push!(E1, ery1)
    push!(E3, ery3)
    push!(E4, ery4)
    push!(E5, ery5)
end

#save data
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPECE.csv",  DataFrame(time = eval(Symbol("T1")), error = eval(Symbol("E1"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIEX.csv",  DataFrame(time = eval(Symbol("T3")), error = eval(Symbol("E3"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPITrap.csv",  DataFrame(time = eval(Symbol("T4")), error = eval(Symbol("E4"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_MTPIRect.csv",  DataFrame(time = eval(Symbol("T5")), error = eval(Symbol("E5"))))

