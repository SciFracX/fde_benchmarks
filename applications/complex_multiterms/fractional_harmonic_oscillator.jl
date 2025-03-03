using FractionalDiffEq, CairoMakie
tspan = (0.0, 50.0)
α = -1; β = 1.2; γ = 4;
rightfun(u, p, t) = 0
prob1 = MultiTermsFODEProblem([1, α, β, γ], [2, 1, 0.75, 0], rightfun, [-0.5, -0.5], tspan)
prob2 = MultiTermsFODEProblem([1, α, β, γ], [2, 1, 0.78, 0], rightfun, [-0.5, -0.5], tspan)
prob3 = MultiTermsFODEProblem([1, α, β, γ], [2, 1, 0.81695, 0], rightfun, [-0.5, -0.5], tspan)
prob4 = MultiTermsFODEProblem([1, α, β, γ], [2, 1, 0.85, 0], rightfun, [-0.5, -0.5], tspan)
prob5 = MultiTermsFODEProblem([1, α, β, γ], [2, 1, 0.95, 0], rightfun, [-0.5, -0.5], tspan)

sol1 = solve(prob1, MTPECE(), dt=0.01)
sol2 = solve(prob2, MTPECE(), dt=0.01)
sol3 = solve(prob3, MTPECE(), dt=0.01)
sol4 = solve(prob4, MTPECE(), dt=0.01)
sol5 = solve(prob5, MTPECE(), dt=0.01)

sol1t = sol1.t; sol1u = reduce(vcat, sol1.u);
sol2t = sol1.t; sol2u = reduce(vcat, sol2.u);
sol3t = sol1.t; sol3u = reduce(vcat, sol3.u);
sol4t = sol1.t; sol4u = reduce(vcat, sol4.u);
sol5t = sol1.t; sol5u = reduce(vcat, sol5.u);

ASPECT_RATIO = 1.0
WIDTH = 1200
HEIGHT = round(Int, WIDTH * ASPECT_RATIO)
STROKEWIDTH = 2.5
fig = Figure(; size = (WIDTH, HEIGHT))
ax = Axis(fig[1, 1])
lines!(ax, sol1t, sol1u)
Label(fig[0, :], L"$\theta=0.75$", fontsize = 24, tellwidth = false, font = :bold)
ax = Axis(fig[3, 1])
lines!(ax, sol2t, sol2u)
Label(fig[2, :], L"$\theta=0.78$", lineheight = 0.1, fontsize = 24, tellwidth = false, font = :bold)
ax = Axis(fig[5, 1])
lines!(ax, sol3t, sol3u)
Label(fig[4, :], L"$\theta=0.81695$", fontsize = 24, tellwidth = false, font = :bold)
ax = Axis(fig[7, 1])
lines!(ax, sol4t, sol4u)
Label(fig[6, :], L"$\theta=0.9$", fontsize = 24, tellwidth = false, font = :bold)
ax = Axis(fig[9, 1])
lines!(ax, sol5t, sol5u)
Label(fig[8, :], L"$\theta=0.95$", fontsize = 24, tellwidth = false, font = :bold)
fig
save("fractional_harmonic_oscillator.svg", fig)

Plots.plot(sol1, lw=2, legend=:bottomright, label = "u")
Plots.plot(sol2, lw=2, legend=:bottomright, label = "u")
Plots.plot(sol3, lw=2, legend=:bottomright, label = "u")
Plots.plot(sol4, lw=2, legend=:bottomright, label = "u")
Plots.plot(sol5, lw=2, legend=:bottomright, label = "u")