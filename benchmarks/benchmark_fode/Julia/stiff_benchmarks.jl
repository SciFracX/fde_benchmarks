using FractionalDiffEq
using SpecialFunctions
using FdeSolver
using DataFrames
using BenchmarkTools
using CSV

beta = [0.5, 0.5, 0.5]
β = 0.5
#σ₁ = β; σ₂ = 2β; σ₃ = 1+β; σ₄ = 5β; σ₅ = 2; σ₆ = 2+β
σ = [β, 2β, 1+β, 5β, 2, 2+β]
a1 = 0.5; a2 = 0.8; a3 = 1; a4 = 1; a5 = 1; a6 = 1

A = [-10000 0 1;
     -0.05 -0.08 -0.2;
     1 0 -1]
B = [-0.6 0 0.2;
     -0.1 -0.2 0;
     0 -0.5 -0.8]
u0 = [1, 1, 1]
Γ(k) = gamma(σ[k]+1)/gamma(σ[k]-β+1)
g(t) = [a1*Γ(1)*t^(σ[1]-β) + a2*Γ(2)*t^(σ[2]-β); a3*Γ(3)*t^(σ[3]-β) + a4*Γ(4)*t^(σ[4]-β); a5*Γ(5)*t^(σ[5]-β) + a6*Γ(6)*t^(σ[6]-β)] - (A + B) * [a1*t^σ[1] + a2*t^σ[2] + u0[1]; a3*t^σ[3] + a4*t^σ[4] + u0[2]; a5*t^σ[5] + a6*t^σ[6] + u0[3]]
function f!(du, u, p, t)
    du .= A*u + B*u + g(t)
end
function jac(J, u, p, t)
    J .= A+B
end
tspan = (0.0, 1.0)
analytic_solution(u, p, t) = [a1*t^σ[1] + a2*t^σ[2] + u0[1]; a3*t^σ[3] + a4*t^σ[4] + u0[2]; a5*t^σ[5] + a6*t^σ[6] + u0[3]]
fun1 = ODEFunction(f!, analytic = analytic_solution)
prob = FODEProblem(fun1, beta, u0, tspan)
fun2 = ODEFunction(f!, jac = jac, analytic = analytic_solution)
prob1 = FODEProblem(fun2, beta, u0, tspan)


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
    t2= @benchmark solve($(prob1), $(FdeSolverPECE()), dt = $(h)) seconds=1
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
    push!(T6, mean(t3).time / 10^9)
    push!(T7, mean(t4).time / 10^9)
    push!(T8, mean(t5).time / 10^9)

    #computing the error
    sol1 = solve(prob, FdeSolverPECE(), dt = h);
    sol2 = solve(prob1, FdeSolverPECE(), dt = h)
    sol3 =  solve(prob, PECE(), dt = h);
    sol4 =  solve(prob, PITrap(), dt = h);
    sol5 =  solve(prob, PIRect(), dt = h);
    sol6 =  solve(prob, BDF(), dt = h);
    sol7 =  solve(prob, NewtonGregory(), dt = h);
    sol8 =  solve(prob, Trapezoid(), dt = h);

    ery1=sol1.errors[:final]
    ery2=sol2.errors[:final]
    ery3=sol3.errors[:final]
    ery4=sol4.errors[:final]
    ery5=sol5.errors[:final]
    ery6=sol6.errors[:final]
    ery7=sol7.errors[:final]
    ery8=sol8.errors[:final]

    push!(E1, ery1)
    push!(E2, ery2)
    push!(E3, ery3)
    push!(E4, ery4)
    push!(E5, ery5)
    push!(E6, ery6)
    push!(E7, ery7)
    push!(E8, ery8)
end

#save data
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FdeSolver_PECE.csv",  DataFrame(time = eval(Symbol("T1")), error = eval(Symbol("E1"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FdeSolver_NR.csv",  DataFrame(time = eval(Symbol("T2")), error = eval(Symbol("E2"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PECE.csv",  DataFrame(time = eval(Symbol("T3")), error = eval(Symbol("E3"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PITrap.csv",  DataFrame(time = eval(Symbol("T4")), error = eval(Symbol("E4"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_PIRect.csv",  DataFrame(time = eval(Symbol("T5")), error = eval(Symbol("E5"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_BDF.csv",  DataFrame(time = eval(Symbol("T6")), error = eval(Symbol("E6"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_NewtonGregory.csv",  DataFrame(time = eval(Symbol("T7")), error = eval(Symbol("E7"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/Stiff_FractionalDiffEq_Trapzoid.csv",  DataFrame(time = eval(Symbol("T8")), error = eval(Symbol("E8"))))

