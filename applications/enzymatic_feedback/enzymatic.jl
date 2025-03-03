using FractionalDiffEq, Plots
function EnzymeKinetics!(dy, y, ϕ, p, t)
    dy[1] = 10.5-y[1]/(1+0.0005*ϕ[4]^3)
    dy[2] = y[1]/(1+0.0005*ϕ[4]^3)-y[2]
    dy[3] = y[2]-y[3]
    dy[4] = y[3]-0.5*y[4]
end

u0 = [60, 10, 10, 20]; alpha = [0.95, 0.95, 0.95, 0.95]
beta(t) = 0.97 - 1/(ℯ^t+1)
tspan = (0.0, 150.0); tau = 4
h(p, t) = [60, 10, 10, 20]
prob = FDDEProblem(EnzymeKinetics!, beta, u0, h, constant_lags = [tau], tspan)
sol = solve(prob, DelayPECE(), dt=0.01)
plot(sol)