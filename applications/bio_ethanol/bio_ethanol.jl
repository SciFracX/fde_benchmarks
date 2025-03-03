using FractionalDiffEq
using OrdinaryDiffEq
using CairoMakie
using CSV
using DataFrames

function fde_lotka_volterra(du, u, p, t)
    kc, km, ks, kp = p
    du[1] = kc * u[1] * u[2] - km * u[1]
    du[2] = -ks * u[1] * u[2]
    du[3] = kp * u[1] * u[2]
end

function ode_lotka_volterra(du, u, p, t)
    kc, km, ks, kp = p
    du[1] = kc * u[1] * u[2] - km * u[1]
    du[2] = -ks * u[1] * u[2]
    du[3] = kp * u[1] * u[2]
end
order = [0.7769, 0.8767, 0.998]
fde_p = [0.0041, 2.3875*10^(-14), 0.0585, 0.0156]
ode_p = [0.0041, 2.3875*10^(-14), 0.0585, 0.0156]
#u0 = [0.5, 2.5, 90.0]
u0 = [0.5, 90.0, 2.0]
tspan = (0.0, 50.0)

fde_prob = FODEProblem(fde_lotka_volterra, order, u0, tspan, fde_p)
ode_prob = ODEProblem(ode_lotka_volterra, u0, tspan, fit1.param)
fde_sol = solve(fde_prob, PITrap(), dt=1)
ode_sol = solve(ode_prob, Tsit5(), dt=1)

biomass = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/biomass.csv"))
product = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/product.csv"))
substrate = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/substrate.csv"))

fig = Figure()
axb = Axis(fig[1, 1], xlabel = "Time (h)", ylabel = "Biomass (g/L)")
axs = Axis(fig[1, 2], xlabel = "Time (h)", ylabel = "Substrate (g/L)")
axp = Axis(fig[2, 1], xlabel = "Time (h)", ylabel = "Product (g/L)")
scb = CairoMakie.scatter!(axb, biomass[!, 1], biomass[!, 2], color = :blue, label = "Experimental Biomass")
scp = CairoMakie.scatter!(axp, product[!, 1], product[!, 2], color = :green, label = "Experimental Product")
scs = CairoMakie.scatter!(axs, substrate[!, 1], substrate[!, 2], color = :red, label = "Experimental Substrate")
fde_lineb=lines!(axb, fde_sol.t, fde_sol[1, :], color = :blue, label = "Biomass (model)")
fde_lines=lines!(axs, fde_sol.t, fde_sol[2, :], color = :red, label = "Substrate (model)")
fde_linep=lines!(axp, fde_sol.t, fde_sol[3, :], color = :green, label = "Product (model)")
ode_lineb=lines!(axb, ode_sol.t, ode_sol[1, :], color = :blue, label = "Biomass (model)")
ode_lines=lines!(axs, ode_sol.t, ode_sol[2, :], color = :red, label = "Substrate (model)")
ode_linep=lines!(axp, ode_sol.t, ode_sol[3, :], color = :green, label = "Product (model)")
fig

using Optim
using LsqFit
using OrdinaryDiffEq

function lotka_volterra(du, u, p, t)
    kc, km, ks, kp = p
    du[1] = kc * u[1] * u[2] - km * u[1]
    du[2] = -ks * u[1] * u[2]
    du[3] = kp * u[1] * u[2]
end


data = Matrix(hcat(biomass[!, 2], substrate[!, 2], product[!, 2])')
t = [0,2,4,6,8,10,12,16,28,44,48]

function lsq_ode_estimator(time_array, phi)
    tspan = (0, 50)
    ini_cond = data[:, 1]
    oprob = ODEProblem(lotka_volterra, ini_cond, tspan, phi)
    osol  = solve(oprob, Tsit5(), saveat=time_array)
    estimated = reduce(hcat, osol.u)
    return vec(estimated')
end

function lsq_fde_estimator(time_array, phi)
    tspan = (0, 50)
    ini_cond = data[:, 1]
    oprob = FODEProblem(lotka_volterra, order, ini_cond, tspan, phi)
    osol  = solve(oprob, PITrap(), dt = 1)
    fde_sol = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]
    estimated = reduce(hcat, fde_sol)
    return vec(estimated')
end

p0_ode = [0.5, 0.5, 0.5, 0.5]
p0_fde = [0.05, 0.0, 0.05, 0.05]
lb = [0.0, 0.0, 0.0, 0.0]
ub = [0.1, 0.1, 0.1, 0.1]
fit1 = curve_fit(lsq_ode_estimator, t, vec(data), p0_ode)
fit2 = curve_fit(lsq_fde_estimator, t, vec(data), p0_fde, lower=lb, upper=ub)
using Plots



function floudas_one(dz_dt, z, phi, t)
    r_1 = phi[1]*z[1]
    r_2 = phi[2]*z[2]

    dz_dt[1] = -r_1
    dz_dt[2] = r_1 - r_2
end

ode_fun = floudas_one

data = Float64[1.0 0.57353 0.328937 0.188654 0.108198 0.0620545 0.0355906 0.0204132 0.011708 0.00671499; 
        0.0 0.401566 0.589647 0.659731 0.666112 0.639512 0.597179 0.54867 0.499168 0.451377]
t = Float64[0.0, 0.111111, 0.222222, 0.333333, 0.444444, 0.555556, 0.666667, 0.777778, 0.888889, 1.0]

function lsq_ss_estimator(time_array, phi)
    tspan = (t[1], t[end])
    ini_cond = data[:,1]
    oprob = ODEProblem(ode_fun, ini_cond, tspan, phi)
    osol  = solve(oprob, Tsit5(), saveat=t)
    estimated = reduce(hcat, osol.u)
    return vec(estimated)
end

p0 = [5., 5.]
fit = curve_fit(lsq_ss_estimator, t, vec(data), p0)