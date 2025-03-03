using FractionalDiffEq, Plots
tspan = (0.0, 20.0)
rhs(u, p, t) = ifelse(t > 1, 8, 0)
u0 = [0.0, 0.0]
prob = MultiTermsFODEProblem([1, 1/2, 1/2], [2, 1.5, 0], rhs, u0, tspan)
sol = solve(prob, MTPECE(), dt=0.01)
Plots.plot(sol, label = "u", xlabel = "t", ylabel = "u(t)", title = "Bagley Torvik Equation", lw = 3, legend = :topright)
savefig("bagley_torvik.svg")