using FractionalDiffEq, FdeSolver, Plots, CSV, SpecialFunctions, BenchmarkTools
using DataFrames
#= MATLAB scripts
alpha = [0.5, 0.2, 0.6];
f_fun = @(t,y) [(((y(2)-0.5).*(y(3)-0.3)).^(1/6) + sqrt(t))/sqrt(pi);
        gamma(2.2)*(y(1)-1);
        gamma(2.8)/gamma(2.2)*(y(2)-0.5)];
J_fun = @(t,y) [0, (y(2)-0.5).^(-5/6).*(y(3)-0.3).^(1/6)/6/sqrt(pi), (y(2)-0.5).^(1/6).*(y(3)-0.3).^(-5/6)/6/sqrt(pi);
        gamma(2.2), 0 , 0;
        0 , gamma(2.8)/gamma(2.2) , 0];
t0 = 0;
T = 5;
y0 = [ 1 ; 0.500000001 ; 0.300000001 ] ;
h = 0.01
[t, y] = FDE_PI1_Im(alpha,f_fun,J_fun,t0,T,y0,h) ;

exa = @(t) [t+1, t^1.2+0.5, t^1.8+0.3];
exact = exa(t1);
=#

function f!(du, u, p, t)
    du[1] = 1/sqrt(pi)*(((u[2]-0.5)*(u[3]-0.3))^(1/6) + sqrt(t))
    du[2] = gamma(2.2)*(u[1]-1)
    du[3] = gamma(2.8)/gamma(2.2)*(u[2]-0.5)
end
analytical_solution(u, p, t) = [t+1, t^1.2+0.5, t^1.8+0.3]
u0 = [1.0, 0.5, 0.3]
tspan = (0.0, 5.0)
order = [0.5, 0.2, 0.6]
fun = ODEFunction(f!, analytic = analytical_solution)
prob = FODEProblem(fun, order, u0, tspan)

# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[]
E3 = Float64[];T3 = Float64[];E4 = Float64[];T4 = Float64[];
E5 = Float64[];T5 = Float64[];E6 = Float64[];T6 = Float64[];
h = Float64[]

for n in range(3, 7)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computing the time
    t1= @benchmark solve($(prob), $(FdeSolverPECE()), dt = $(h)) seconds=1
    #t2= @benchmark FDEsolver(F, $(tSpan), $(y0), $(β), $(par), JF = JF, h=$(h)) seconds=1
    t3= @benchmark solve($(prob), $(PECE()), dt = $(h)); seconds=1
    t4 = @benchmark solve($(prob), $(PITrap()), dt = $(h)); seconds=1
    t5 = @benchmark solve($(prob), $(PIRect()), dt = $(h)); seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    #push!(T2, mean(t2).time / 10^9)
    push!(T3, mean(t3).time / 10^9)
    push!(T4, mean(t4).time / 10^9)
    push!(T5, mean(t5).time / 10^9)

    #computing the error
    sol1 = solve(prob, FdeSolverPECE(), dt = h);
    #_, y2 = FDEsolver(F, tSpan, y0, β, par, JF = JF, h=h)
    sol3 =  solve(prob, PECE(), dt = h);
    sol4 =  solve(prob, PITrap(), dt = h);
    sol5 =  solve(prob, PIRect(), dt = h);

    ery1=sol1.errors[:final]
    #ery2=sol2.errors[:final]
    ery3=sol3.errors[:final]
    ery4=sol4.errors[:final]
    ery5=sol5.errors[:final]

    push!(E1, ery1)
    push!(E3, ery3)
    push!(E4, ery4)
    push!(E5, ery5)
end

#save data
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FdeSolver_PECE.csv",  DataFrame(time = eval(Symbol("T1")), error = eval(Symbol("E1"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PECE.csv",  DataFrame(time = eval(Symbol("T3")), error = eval(Symbol("E3"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PITrap.csv",  DataFrame(time = eval(Symbol("T4")), error = eval(Symbol("E4"))))
CSV.write("/Users/quqingyu/SciFracX/paper/benchmarks/data/FractionalDiffEq_PIRect.csv",  DataFrame(time = eval(Symbol("T5")), error = eval(Symbol("E5"))))

