using FractionalDiffEq
using Plots

function fun!(du, u, p, t)
    du[1:end-1] .= u[2:end]
    du[end] = -t^0.1*mittleff(1, 1.545, t)/mittleff(1, 1.445, t)*ℯ^t*u[1]u[112] + ℯ^(-2*t) - u[201]^2
end

tspan = (0.0, 1.0)
alpha = 0.005*ones(291)
u0 = zeros(291)
u0[1] = 1.0
u0[201] = -1.0
prob = FODEProblem(fun!, alpha, u0, tspan)
sol = solve(prob, PITrap(), dt=0.001)

plot(sol, vars = (0, 1))
plot!(sol.t, exp.(-sol.t))