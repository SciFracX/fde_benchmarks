using FractionalDiffEq, FdeSolver, Plots, CSV, SpecialFunctions, BenchmarkTools
using DataFrames, LinearAlgebra

f!(y, par, t) = -10*y

u0 = 1.0
tspan = (0.0, 5.0)
order = 0.8
analytical_solution(u, p, t) = [mittleff(order, -10*t^order)]
fun = ODEFunction(f!)#, analytic = analytical_solution)
prob = FODEProblem(fun, order, u0, tspan)
sol = solve(prob, PIEX(), dt=0.001)
# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[]
E3 = Float64[];T3 = Float64[];E4 = Float64[];T4 = Float64[];
E5 = Float64[];T5 = Float64[];E6 = Float64[];T6 = Float64[];
E7 = Float64[];T7 = Float64[];E8 = Float64[];T8 = Float64[];
h = Float64[]

for n in range(3, 7)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computing the time
    t1= @benchmark solve($(prob), $(FdeSolverPECE()), dt = $(h)) seconds=1
    t2= @benchmark solve($(prob), $(PIEX()), dt = $(h)); seconds=1
    t3= @benchmark solve($(prob), $(PECE()), dt = $(h)); seconds=1
    t4 = @benchmark solve($(prob), $(PITrap()), dt = $(h)); seconds=1
    t5 = @benchmark solve($(prob), $(PIRect()), dt = $(h)); seconds=1
    t6= @benchmark solve($(prob), $(BDF()), dt = $(h)); seconds=1
    t7 = @benchmark solve($(prob), $(NewtonGregory()), dt = $(h)); seconds=1
    t8 = @benchmark solve($(prob), $(Trapezoid()), dt = $(h)); seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
    push!(T3, mean(t3).time / 10^9)
    push!(T4, mean(t4).time / 10^9)
    push!(T5, mean(t5).time / 10^9)
    push!(T6, mean(t6).time / 10^9)
    push!(T7, mean(t7).time / 10^9)
    push!(T8, mean(t8).time / 10^9)

    #computing the error
    sol1 = solve(prob, FdeSolverPECE(), dt = h);
    sol2 =  solve(prob, PIEX(), dt = h);
    sol3 =  solve(prob, PECE(), dt = h);
    sol4 =  solve(prob, PITrap(), dt = h);
    sol5 =  solve(prob, PIRect(), dt = h);
    sol6 =  solve(prob, BDF(), dt = h);
    sol7 =  solve(prob, NewtonGregory(), dt = h);
    sol8 =  solve(prob, Trapezoid(), dt = h);

    exa = analytical_solution.(0, 0, sol3.t)

    ery1=norm(sol1.u - exa, 2)
    ery3=norm(sol2.u - exa, 2)
    ery3=norm(sol3.u - exa, 2)
    ery4=norm(sol4.u - exa, 2)
    ery5=norm(sol5.u - exa, 2)
    ery6=norm(sol6.u - exa, 2)
    ery7=norm(sol7.u - exa, 2)
    ery8=norm(sol8.u - exa, 2)

    push!(E1, ery1)
    push!(E2, ery1)
    push!(E3, ery3)
    push!(E4, ery4)
    push!(E5, ery5)
    push!(E6, ery6)
    push!(E7, ery7)
    push!(E8, ery8)
end

#save data
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FdeSolver_PECE.csv",  DataFrame(time = eval(Symbol("T1")), error = eval(Symbol("E1"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_PIEX.csv",  DataFrame(time = eval(Symbol("T2")), error = eval(Symbol("E2"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_PECE.csv",  DataFrame(time = eval(Symbol("T3")), error = eval(Symbol("E3"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_PITrap.csv",  DataFrame(time = eval(Symbol("T4")), error = eval(Symbol("E4"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_PIRect.csv",  DataFrame(time = eval(Symbol("T5")), error = eval(Symbol("E5"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_BDF.csv",  DataFrame(time = eval(Symbol("T6")), error = eval(Symbol("E6"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_NewtonGregory.csv",  DataFrame(time = eval(Symbol("T7")), error = eval(Symbol("E7"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Linear_Single_FractionalDiffEq_Trapzoid.csv",  DataFrame(time = eval(Symbol("T8")), error = eval(Symbol("E8"))))

