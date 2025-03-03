using FractionalDiffEq, Plots

using FractionalDiffEq, Plots
function chua!(du, x, p, t)
    a, b, c, m0, m1 = p
    du[1] = a*(x[2]-x[1]-(m1*x[1]+m0*(abs(x[1]+1)-abs(x[1]-1))))
    du[2] = x[1]-x[2]+x[3]
    du[3] = -b*x[2]-c*x[3]
end
α = [0.93, 0.99, 0.92];
x0 = [0.2; -0.1; 0.1];
tspan = (0, 100);
p = [10.725, 10.593, 0.268, -0.1927, -0.7872]
prob = FODEProblem(chua!, α, x0, tspan, p)
sol = solve(prob, PECE(), dt=0.01)
plot(sol, vars=(1,2,3), linewidth=1.3, title="Chua System", xaxis = "u₁", yaxis="u₂", zaxis="u₃", legend=false)
savefig("chua.svg")