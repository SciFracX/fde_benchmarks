# Dataset subset
using CSV
using Optim
using DataFrames
using OrdinaryDiffEq
using FractionalDiffEq
using StatsBase
using CairoMakie

# parameters
kc=0.0041 # Transmission coeﬃcient from infected individuals
km=2.3875*10^(-14) # Relative transmissibility of hospitalized patients
ks=0.0585 # Transmission coeﬃcient due to super-spreaders
kp=0.0156 # Rate at which exposed become infectious
function lotka_volterra(du, u, p, t)
    kc, km, ks, kp = p
    du[1] = kc * u[1] * u[2] - km * u[1]
    du[2] = -ks * u[1] * u[2]
    du[3] = kp * u[1] * u[2]
end
biomass = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/biomass.csv"))
product = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/product.csv"))
substrate = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/substrate.csv"))

X0=[0.5, 90.0, 2.5] # initial values
tspan=(0, 50) # time span [initial time, final time]
data = Matrix(hcat(biomass[!, 2], substrate[!, 2], product[!, 2])')
t = [0,2,4,6,8,10,12,16,28,44,48]
par=[kc, km, ks, kp] # parameters
## optimazation of β for integer order model

function loss_1(b)# loss function
	par=copy(b)
    prob = ODEProblem(lotka_volterra, X0, tspan, par)
	#_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    sol = solve(prob, Tsit5(), saveat=t)
    #appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    appX = reduce(hcat, sol.u)
    rmsd(data, appX; normalize=:true) # Normalized root-mean-square error
end

function loss_2(b)# loss function
	par=b[4:end]
    order = b[1:3]
    prob = FODEProblem(lotka_volterra, order, X0, tspan, par)
	#_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    osol = solve(prob, PITrap(), dt = 1)
    #appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    app = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]
    appX = reduce(hcat, app)
    rmsd(data, appX; normalize=:true) # Normalized root-mean-square error
end

function moser_luong_lotka_volterra(du, u, p, t)
    um, ks, kp, ms, yxs, yps, alpha, beta = p
    n=1;m=9
    mu = um*u[2]^n/(ks+u[2]^n)*(1-(u[3]/kp)^m)/(1+exp(-t+1))
    du[1] = mu*u[1]
    du[2] = -(1/yxs + 1/yps) * mu * u[1] - ms * u[1]
    du[3] = alpha*yps*mu*u[1] + beta*u[1]
end
function loss_3(b)# loss function
	par=copy(b)
    prob = ODEProblem(moser_luong_lotka_volterra, X0, tspan, par)
	#_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    sol = solve(prob, Tsit5(), saveat=t)
    #appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    appX = reduce(hcat, sol.u)
    rmsd(data, appX; normalize=:true) # Normalized root-mean-square error
end

p_lo_1=[0.003, 0.0, 0.05, 0.015] #lower bound for β
p_up_1=[0.005, 0.0005, 0.06, 0.016] # upper bound for β
p_vec_1=[0.004, 0.0001, 0.055, 0.0155] #  initial guess for β
Res1=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p1=vcat(Optim.minimizer(Res1))



p_lo_2=[0.7, 0.8, 0.9, 0.003, 0.0, 0.05, 0.015] #lower bound for β
p_up_2=[0.8, 0.9, 1.0, 0.005, 0.0005, 0.06, 0.016] # upper bound for β
p_vec_2=[0.75, 0.85, 0.95, 0.004, 0.0001, 0.055, 0.0155] #  initial guess for β
Res2=optimize(loss_2,p_lo_2,p_up_2,p_vec_2,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p2=vcat(Optim.minimizer(Res2))





p_lo_3=[0.3, 15, 120, 0.0, 0.0, 0.4, 20, 0.0] #lower bound for β
p_up_3=[0.4, 25, 140, 0.1, 0.1, 0.5, 25, 0.05] # upper bound for β
p_vec_3=[0.35, 19, 125, 0.01, 0.01, 0.45, 24, 0.03] #  initial guess for β
Res3=optimize(loss_3,p_lo_3,p_up_3,p_vec_3,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p3=vcat(Optim.minimizer(Res3))


using Plots
aprob1 = ODEProblem(lotka_volterra, X0, tspan, p1)
ode_sol = solve(aprob1, Tsit5())

aprob2 = FODEProblem(lotka_volterra, p2[1:3], X0, tspan, p2[4:end])
fde_sol = solve(aprob2, PITrap(), dt = 0.01)





fig = Figure(size = (700, 490))
axb = Axis(fig[1, 1], xlabel = "Time (h)", ylabel = "Concentration of biomass (g/L)")

axs = Axis(fig[1, 2], xlabel = "Time (h)", ylabel = "Concentration of substrate (g/L)")

axp = Axis(fig[2, 1], xlabel = "Time (h)", ylabel = "Concentration of product (g/L)")

ax_bar = Axis(fig[2, 2], xticks = (1:2, ["Integer order model", "Fractional order model"]), title = "Normalized RMSE")
ax_bar.xlabelsize=3
scb = CairoMakie.scatter!(axb, biomass[!, 1], biomass[!, 2], color = :blue, label = "Experimental Biomass")
scp = CairoMakie.scatter!(axp, product[!, 1], product[!, 2], color = :green, label = "Experimental Product")
scs = CairoMakie.scatter!(axs, substrate[!, 1], substrate[!, 2], color = :red, label = "Experimental Substrate")
fde_lineb=lines!(axb, fde_sol.t, fde_sol[1, :], color = :blue, label = "Biomass (model)")
fde_lines=lines!(axs, fde_sol.t, fde_sol[2, :], color = :red, label = "Substrate (model)")
fde_linep=lines!(axp, fde_sol.t, fde_sol[3, :], color = :green, label = "Product (model)")
ode_lineb=lines!(axb, ode_sol.t, ode_sol[1, :], linestyle = :dash, color = :blue, label = "Biomass (model)")
ode_lines=lines!(axs, ode_sol.t, ode_sol[2, :], linestyle = :dash, color = :red, label = "Substrate (model)")
ode_linep=lines!(axp, ode_sol.t, ode_sol[3, :], linestyle = :dash, color = :green, label = "Product (model)")
barplot!(ax_bar, [1, 2], [err1, err2])

axislegend(axb, [fde_lineb, ode_lineb], ["Fractional order", "Integer order"], position = :rb, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(axs, [fde_lines, ode_lines], ["Fractional order", "Integer order"], position = :rt, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(axp, [fde_linep, ode_linep], ["Fractional order", "Integer order"], position = :rb, labelsize=10, rowgap = -5, patchsize = (12, 22))
fig
save("bio_ethanol.svg", fig)






















aprob3 = ODEProblem(moser_luong_lotka_volterra, X0, tspan, p3)
sol = solve(aprob1, Tsit5())
plot(sol)
scatter!(t, data[1, :], label="biomass")
scatter!(t, data[2, :], label="substrate")
scatter!(t, data[3, :], label="product")

sol1_err = solve(aprob1, Tsit5(), saveat=t)
osol = solve(aprob2, PITrap(), dt = 1)
sol3_err = solve(aprob3, Tsit5(), saveat=t)

sol2_err = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]

err1 = rmsd(data, reduce(hcat, sol1_err.u); normalize=true)
err2 = rmsd(data, reduce(hcat, sol2_err); normalize=true)
err3 = rmsd(data, reduce(hcat, sol3_err.u); normalize=true)

bar([1, 2, 3], [err1, err2, err3], label="error", title="Normalized RMSE", xticks=([1, 2, 3], ["ODE", "FDE", "Moser-Luong"]))